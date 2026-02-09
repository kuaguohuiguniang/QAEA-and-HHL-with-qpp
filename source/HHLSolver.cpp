#include "HHLSolver.hpp"
#include "cmath"
#include <cassert>
#include <iostream>
#include <string>
#include <tuple>
#include <numeric>

void HHLSolver::apply_QPE(qpp::ket& state, const qpp::cmat& A, 
                   const std::vector<qpp::idx>& clock_regs, 
                   const std::vector<qpp::idx>& target_regs) {
    std::cout << "Applying QPE" << std::endl;
    for (idx i = 0; i < static_cast<idx>(clock_regs.size()); ++i) {
        // apply Hadamard on counting qubits
        state = apply(state, gt.H, {clock_regs[i]});
        std::cout << "H" << clock_regs[i] << " ";
    }
    std::cout << std::endl;
    idx powU = 1;
    cmat U = expm(-1_i * A * t0);
    idx nq_c = static_cast<idx>(clock_regs.size());
    for (idx i = 0; i < static_cast<idx>(clock_regs.size()); ++i) {
        // apply controlled-U^(2^i) with clock qubit as control and target_regs as target
        state = applyCTRL(state, U, {clock_regs[nq_c - i - 1]}, target_regs);
        U = powm(U, 2);
        powU *= 2;
    }
    state = applyTFQ(state, clock_regs);
}

void HHLSolver::apply_rotation(qpp::ket& state, 
                               const std::vector<qpp::idx>& clock_regs, 
                               qpp::idx ancilla_reg) {

    idx nq_c = static_cast<idx>(clock_regs.size());
    idx N_dim = 1LL << nq_c; // 2^n_clock

    std::cout << "Applying Conditional Rotation (Ancilla: " << ancilla_reg << ")..." << std::endl;

    // Iterate through every possible basis state 'k' of the clock register
    for (idx k = 0; k < N_dim; ++k) {
        
        // --- 1. Calculate the Eigenvalue lambda corresponding to integer k ---
        
        // Handle 2's complement for negative phases
        // If k >= N/2, it represents a negative integer (k - N)
        long long signed_k = k;
        if (k >= (N_dim / 2)) {
            signed_k = static_cast<long long>(k) - static_cast<long long>(N_dim);
        }

        double lambda = -(2.0 * pi * signed_k) / (t0 * static_cast<double>(N_dim));

        // --- 2. Calculate Rotation Angle theta ---
        // theta = 2 * arcsin(C / lambda)
        
        double theta = 0.0;
        
        // Avoid division by zero (singularity at lambda=0)
        if (std::abs(lambda) >= 1e-9) { 
            double ratio = C / lambda;
            
            if (ratio > 1.0) ratio = 1.0;
            else if (ratio < -1.0) ratio = -1.0;
            
            theta = 2.0 * std::asin(ratio);
        }

        if (std::abs(theta) < 1e-9) continue;


        // --- 3. Construct the Control Logic ---
        // Apply RY(theta) ONLY if the clock register is in state |k>.
        // Since standard controls activate on |1>, we flip the |0> bits of k using X gates.

        std::vector<idx> zero_bit_qubits;
        
        for (idx bit = 0; bit < nq_c; ++bit) {
            bool is_one = (k >> bit) & 1;

            idx physical_qubit = clock_regs[nq_c - 1 - bit];

            if (!is_one) {
                // If the bit in k is 0, apply X to make it 1 for the control
                state = apply(state, gt.X, {physical_qubit});
                zero_bit_qubits.push_back(physical_qubit);
            }
        }

        // --- 4. Apply Multi-Controlled Ry Gate ---
        // Controls: All clock qubits. Target: Ancilla.
        state = applyCTRL(state, gt.RY(theta), clock_regs, {ancilla_reg});

        // --- 5. Restore State (Uncompute X gates) ---
        for (idx q : zero_bit_qubits) {
            state = apply(state, gt.X, {q});
        }
    }
}

void HHLSolver::apply_inv_QPE(qpp::ket& state, const qpp::cmat& A, 
                              const std::vector<qpp::idx>& clock_regs, 
                              const std::vector<qpp::idx>& target_regs) {

    std::cout << "Applying Inverse QPE (Uncomputation)..." << std::endl;

    // 1. Forward QFT (Reverse of the last step of QPE)
    state = applyQFT(state, clock_regs);

    // 2. Inverse Controlled-Unitaries (Reverse of the middle step)
    // Forward used: exp(-i * A * t0)
    // Inverse uses: exp(+i * A * t0)
    
    idx nq_c = static_cast<idx>(clock_regs.size());
    
    cmat U_inv = expm(1_i * A * t0);
    
    long long powU = 1;
    
    for (idx i = 0; i < nq_c; ++i) {
        idx ctrl_qubit = clock_regs[nq_c - i - 1];
        state = applyCTRL(state, U_inv, {ctrl_qubit}, target_regs);
        U_inv = powm(U_inv, 2); 
        powU *= 2;
    }

    // 3. Hadamard Gates (Reverse of the first step)
    for (idx i = 0; i < nq_c; ++i) {
        state = apply(state, gt.H, {clock_regs[i]});
    }
}

HHLSolver::HHLSolver(int precision, double evolution_time, double rotation_const)
    : n_clock(precision), t0(evolution_time), C(rotation_const) {};

qpp::ket HHLSolver::solve(const LinearSystem& sys) {
    std::cout << "Starting HHL Solver..." << std::endl;

    idx nq_c = static_cast<idx>(n_clock);
    idx nq_tar = static_cast<idx>(std::ceil(std::log2(sys.N))); 
    idx nq_a = 1;                         
    idx nq = nq_c + nq_tar + nq_a;        

    // Define register indices
    std::vector<idx> clock_regs(nq_c);
    std::iota(clock_regs.begin(), clock_regs.end(), 0);

    std::vector<idx> target_regs(nq_tar);
    std::iota(target_regs.begin(), target_regs.end(), nq_c);

    idx ancilla_reg = nq - 1;

    // Initialize state |0>^n_clock ⊗ |b> ⊗ |0>_ancilla
    qpp::ket state = kron(qpp::mket(std::vector<idx>(nq_c, 0)), sys.b, qpp::mket({0}));
    std::cout << "Initial state prepared." << std::endl;

    // 1. Apply Quantum Phase Estimation
    apply_QPE(state, sys.A, clock_regs, target_regs);

    // 2. Apply Conditional Rotation
    apply_rotation(state, clock_regs, ancilla_reg);

    // 3. Apply Inverse Quantum Phase Estimation
    apply_inv_QPE(state, sys.A, clock_regs, target_regs);

    std::cout << "HHL Circuit completed. Measuring Ancilla..." << std::endl;

    // 4. Extract Solution Vector

    auto measured = qpp::measure_seq(state, {ancilla_reg});
    idx outcome = std::get<0>(measured)[0];
    auto states = std::get<2>(measured);

    if (outcome == 1) {
        std::cout << "Success! Ancilla measured 1.\n";

        qpp::ket post = states;                 // collapsed full ket

        std::vector<idx> dims_total(nq-1, 2);
        std::vector<idx> qubits_to_discard = clock_regs; // discard clock qubits
        qpp::cmat rho_target = qpp::ptrace(post, qubits_to_discard, dims_total);
        Eigen::SelfAdjointEigenSolver<qpp::cmat> es(rho_target);
        qpp::ket x = es.eigenvectors().col(rho_target.rows() - 1);

        // optional: fix global phase
        if (std::abs(x(0)) > 1e-12)
            x *= std::conj(x(0)) / std::abs(x(0));

        return x;
    } else {
        std::cout << "Failure. Ancilla measured 0.\n";
        return qpp::ket::Zero(1ULL << nq_tar);
    }
}