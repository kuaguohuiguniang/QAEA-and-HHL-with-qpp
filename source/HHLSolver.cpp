#include "HHLSolver.hpp"
#include "cmath"
#include <cassert>
#include <iostream>
#include <string>
#include <tuple>
#include <numeric>

/**
 * @brief Applies the quantum phase estimation stage of HHL.
 *
 * The clock register is placed into superposition, controlled powers of
 * exp(-iAt0) are applied to the target register, and the inverse QFT maps the
 * accumulated phases into computational-basis clock states.
 */
void HHLSolver::apply_QPE(qpp::ket& state, const qpp::cmat& A, 
                   const std::vector<qpp::idx>& clock_regs, 
                   const std::vector<qpp::idx>& target_regs) {
    std::cout << "Applying QPE" << std::endl;
    for (idx i = 0; i < static_cast<idx>(clock_regs.size()); ++i) {
        /// Apply Hadamard gates to the clock qubits.
        state = apply(state, gt.H, {clock_regs[i]});
        std::cout << "H" << clock_regs[i] << " ";
    }
    std::cout << std::endl;
    idx powU = 1;
    cmat U = expm(-1_i * A * t0);
    idx nq_c = static_cast<idx>(clock_regs.size());
    for (idx i = 0; i < static_cast<idx>(clock_regs.size()); ++i) {
        /// Apply controlled-U^(2^i) to the target register.
        state = applyCTRL(state, U, {clock_regs[nq_c - i - 1]}, target_regs);
        U = powm(U, 2);
        powU *= 2;
    }
    state = applyTFQ(state, clock_regs);
}

/**
 * @brief Applies the eigenvalue-dependent ancilla rotation used by HHL.
 *
 * Each possible clock-register value is interpreted as a discretized
 * eigenvalue estimate. The ancilla is then rotated by an angle proportional to
 * C/lambda, encoding the reciprocal eigenvalue information needed for the
 * linear-system solution.
 */
void HHLSolver::apply_rotation(qpp::ket& state, 
                               const std::vector<qpp::idx>& clock_regs, 
                               qpp::idx ancilla_reg) {

    idx nq_c = static_cast<idx>(clock_regs.size());
    idx N_dim = 1LL << nq_c; ///< Number of computational states in the clock register.

    std::cout << "Applying Conditional Rotation (Ancilla: " << ancilla_reg << ")..." << std::endl;

    /// Iterate through every possible basis state k of the clock register.
    for (idx k = 0; k < N_dim; ++k) {
        
        /// Step 1: calculate the eigenvalue lambda corresponding to integer k.
        
        /// Decode negative phases using a two's-complement-style convention.
        long long signed_k = k;
        if (k >= (N_dim / 2)) {
            signed_k = static_cast<long long>(k) - static_cast<long long>(N_dim);
        }

        double lambda = -(2.0 * pi * signed_k) / (t0 * static_cast<double>(N_dim));

        /// Step 2: calculate theta = 2 * arcsin(C / lambda).
        
        double theta = 0.0;
        
        /// Avoid division by zero at lambda = 0.
        if (std::abs(lambda) >= 1e-9) { 
            double ratio = C / lambda;
            
            if (ratio > 1.0) ratio = 1.0;
            else if (ratio < -1.0) ratio = -1.0;
            
            theta = 2.0 * std::asin(ratio);
        }

        if (std::abs(theta) < 1e-9) continue;


        /// Step 3: construct controls so RY(theta) fires only on clock state |k>.
        /// Standard controls activate on |1>, so the zero bits of k are flipped.

        std::vector<idx> zero_bit_qubits;
        
        for (idx bit = 0; bit < nq_c; ++bit) {
            bool is_one = (k >> bit) & 1;

            idx physical_qubit = clock_regs[nq_c - 1 - bit];

            if (!is_one) {
                /// Flip a zero bit to make the multi-control condition active.
                state = apply(state, gt.X, {physical_qubit});
                zero_bit_qubits.push_back(physical_qubit);
            }
        }

        /// Step 4: apply the multi-controlled RY gate to the ancilla.
        state = applyCTRL(state, gt.RY(theta), clock_regs, {ancilla_reg});

        /// Step 5: uncompute the X gates used to activate zero controls.
        for (idx q : zero_bit_qubits) {
            state = apply(state, gt.X, {q});
        }
    }
}

/**
 * @brief Uncomputes the phase-estimation register.
 *
 * This applies the inverse of apply_QPE: forward QFT, inverse controlled
 * Hamiltonian evolutions, and final Hadamard gates on the clock register.
 */
void HHLSolver::apply_inv_QPE(qpp::ket& state, const qpp::cmat& A, 
                              const std::vector<qpp::idx>& clock_regs, 
                              const std::vector<qpp::idx>& target_regs) {

    std::cout << "Applying Inverse QPE (Uncomputation)..." << std::endl;

    /// Step 1: apply forward QFT, reversing the final inverse QFT in QPE.
    state = applyQFT(state, clock_regs);

    /// Step 2: apply inverse controlled unitaries exp(+iAt0).
    
    idx nq_c = static_cast<idx>(clock_regs.size());
    
    cmat U_inv = expm(1_i * A * t0);
    
    long long powU = 1;
    
    for (idx i = 0; i < nq_c; ++i) {
        idx ctrl_qubit = clock_regs[nq_c - i - 1];
        state = applyCTRL(state, U_inv, {ctrl_qubit}, target_regs);
        U_inv = powm(U_inv, 2); 
        powU *= 2;
    }

    /// Step 3: apply Hadamard gates to return the clock register to |0...0>.
    for (idx i = 0; i < nq_c; ++i) {
        state = apply(state, gt.H, {clock_regs[i]});
    }
}

HHLSolver::HHLSolver(int precision, double evolution_time, double rotation_const)
    : n_clock(precision), t0(evolution_time), C(rotation_const) {};

/**
 * @brief Runs the complete HHL circuit simulation.
 *
 * The method allocates clock, target, and ancilla registers, prepares
 * |0...0>|b>|0>, applies the HHL circuit stages, measures the ancilla, and
 * extracts the target-register state after successful postselection.
 */
qpp::ket HHLSolver::solve(const LinearSystem& sys) {
    std::cout << "Starting HHL Solver..." << std::endl;

    idx nq_c = static_cast<idx>(n_clock);
    idx nq_tar = static_cast<idx>(std::ceil(std::log2(sys.N))); 
    idx nq_a = 1;                         
    idx nq = nq_c + nq_tar + nq_a;        

    /// Define register indices.
    std::vector<idx> clock_regs(nq_c);
    std::iota(clock_regs.begin(), clock_regs.end(), 0);

    std::vector<idx> target_regs(nq_tar);
    std::iota(target_regs.begin(), target_regs.end(), nq_c);

    idx ancilla_reg = nq - 1;

    /// Equivalent documented form: |0>^n_clock tensor |b> tensor |0>_ancilla.
    qpp::ket state = kron(qpp::mket(std::vector<idx>(nq_c, 0)), sys.b, qpp::mket({0}));
    std::cout << "Initial state prepared." << std::endl;

    /// Step 1: apply quantum phase estimation.
    apply_QPE(state, sys.A, clock_regs, target_regs);

    /// Step 2: apply the conditional ancilla rotation.
    apply_rotation(state, clock_regs, ancilla_reg);

    /// Step 3: uncompute the phase-estimation register.
    apply_inv_QPE(state, sys.A, clock_regs, target_regs);

    std::cout << "HHL Circuit completed. Measuring Ancilla..." << std::endl;

    /// Step 4: measure the ancilla and extract the target-register solution.

    auto measured = qpp::measure_seq(state, {ancilla_reg});
    idx outcome = std::get<0>(measured)[0];
    auto states = std::get<2>(measured);

    if (outcome == 1) {
        std::cout << "Success! Ancilla measured 1.\n";

        qpp::ket post = states; ///< Collapsed full ket after ancilla measurement.

        std::vector<idx> dims_total(nq-1, 2);
        std::vector<idx> qubits_to_discard = clock_regs; ///< Clock qubits to trace out.
        qpp::cmat rho_target = qpp::ptrace(post, qubits_to_discard, dims_total);
        Eigen::SelfAdjointEigenSolver<qpp::cmat> es(rho_target);
        qpp::ket x = es.eigenvectors().col(rho_target.rows() - 1);

        /// Fix the global phase for more readable output.
        if (std::abs(x(0)) > 1e-12)
            x *= std::conj(x(0)) / std::abs(x(0));

        return x;
    } else {
        std::cout << "Failure. Ancilla measured 0.\n";
        return qpp::ket::Zero(1ULL << nq_tar);
    }
}
