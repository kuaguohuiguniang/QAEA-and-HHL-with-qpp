/*#pragma once
#include "input.hpp"
#include "qpp/qpp.hpp"
#include <vector>

class FourierSolver {
private:
    // Parameters defining the precision of the approximation
    int J; // Truncation order (e.g., number of terms in the sum)
    
    // Internal struct to hold one term of the LCU decomposition: c_k * exp(-i * A * t_k)
    struct LCUTerm {
        double time;       // t_k
        std::complex<double> coeff; // c_k
    };

    // Helper to generate the specific coefficients (Fourier series / Riemann sum)
    // You will implement the math from the paper here.
    std::vector<LCUTerm> generate_lcu_terms();

    // --- The Two Main Oracles ---

    // 1. Prepare (V)
    // Initializes the control register into the state sum sqrt(c_k) |k>
    void apply_prepare_V(qpp::ket& state, 
                         const std::vector<qpp::idx>& control_regs,
                         const std::vector<LCUTerm>& terms);

    // 2. Select (U)
    // Applies Controlled-Hamiltonian evolution: sum |k><k| tensor exp(-i * A * t_k)
    void apply_select_U(qpp::ket& state, const qpp::cmat& A,
                        const std::vector<qpp::idx>& control_regs, 
                        const std::vector<qpp::idx>& target_regs,
                        const std::vector<LCUTerm>& terms);

    // 3. Inverse Prepare (V_dagger)
    void apply_prepare_V_dagger(qpp::ket& state, 
                                const std::vector<qpp::idx>& control_regs,
                                const std::vector<LCUTerm>& terms);

public:
    FourierSolver(int precision_order);

    qpp::ket solve(const LinearSystem& sys);
};*/

/*
#pragma once
#include "input.hpp"
#include <qpp/qpp.hpp>
#include <vector>
#include <complex>
#include <cmath>

class FourierSolver {
private:
    // --- Derived System Properties ---
    double kappa;       // Condition number of A (lambda_max / lambda_min)
    double max_eig;     // Spectral norm ||A||

    // --- Algorithm Parameters (Auto-optimized) ---
    int J;              // Truncation order for the Fourier series (Frequency bandwidth)
    int K;              // Number of integration steps (Cutoff term)
    double delta_z;     // Fundamental time step (dt) to avoid aliasing (must satisfy max_eig * dt < pi)
    double delta_y;     // Integration step size (dy)
    double epsilon;     // Desired precision (passed in constructor)

    // Internal struct to hold one term of the Linear Combination of Unitaries (LCU)
    // Represents the operation: coeff * exp(-i * A * time)
    struct LCUTerm {
        double time;                // t_m = j * k * delta_z
        std::complex<double> coeff; // c_m (Combined weight of Fourier coeff and integration weight)
    };

    // --- Helper Functions ---

    // 1. Parameter Optimizer
    // Estimates kappa and max_eig from A, then calculates optimal J, K, delta_z, delta_y
    // scaling primarily as O(kappa * log(1/epsilon))
    void optimize_parameters(const qpp::cmat& A);

    // 2. Term Generator
    // Generates the flattened list of (time, coeff) pairs by iterating j and k
    std::vector<LCUTerm> generate_lcu_terms();

    // 3. Prepare Oracle (V)
    // Initializes the control registers into the superposition state: sum_m sqrt(c_m) |m>
    void apply_prepare_V(qpp::ket& state, 
                         const std::vector<qpp::idx>& control_regs,
                         const std::vector<LCUTerm>& terms);

    // 4. Select Oracle (U)
    // Applies Controlled-Hamiltonian evolution: sum_m |m><m| \otimes exp(-i * A * t_m)
    void apply_select_U(qpp::ket& state, const qpp::cmat& A,
                        const std::vector<qpp::idx>& control_regs, 
                        const std::vector<qpp::idx>& target_regs,
                        const std::vector<LCUTerm>& terms);

    // 5. Inverse Prepare (V_dagger)
    // Uncomputes the coefficient state to allow for interference (V^dagger)
    void apply_prepare_V_dagger(qpp::ket& state, 
                                const std::vector<qpp::idx>& control_regs,
                                const std::vector<LCUTerm>& terms);

public:
    // Constructor
    // precision: The desired error tolerance (e.g., 1e-3). Default provided.
    FourierSolver(double precision = 1e-3);

    // The Main Driver
    // 1. Auto-calculates J, K, delta_z based on sys.A
    // 2. Runs the Prepare-Select-Inverse sequence (LCU Lemma)
    // 3. Measures control registers (post-selection on |0>) and returns the solution vector
    qpp::ket solve(const LinearSystem& sys);
};
*/
#pragma once

#include "input.hpp"
#include <qpp/qpp.hpp>
#include <vector>
#include <complex>

class FourierSolver {
public:
    // ---------------- Public Types ----------------

    struct Result {
        qpp::ket solution;     // Postselected solution state |x>
        double success_prob;   // Probability of success (measuring |0> on control)
    };

    // ---------------- Constructor ----------------

    // epsilon: target operator approximation error
    explicit FourierSolver(double epsilon = 1e-3);

    // ---------------- Main API ----------------

    // Solves A x = b using Fourier-based LCU (simulation-level)
    Result solve(const LinearSystem& sys) const;

private:
    // ---------------- Internal Data Models ----------------

    struct FourierParams {
        double spectral_norm;  // ||A||_2
        double kappa;          // condition number
        int J;
        int K;
        double delta_z;
        double delta_y;
    };

    struct LCUTerm {
        double time;
        std::complex<double> coeff;
        double weight;         // |coeff|
    };

    struct LCUExpansion {
        std::vector<LCUTerm> terms;
        double alpha;          // sum |coeff|
    };

    // ---------------- Classical Preprocessing ----------------

    FourierParams compute_fourier_params(const qpp::cmat& A) const;

    LCUExpansion build_lcu_expansion(const FourierParams& params) const;

    // ---------------- Quantum Subroutines ----------------

    qpp::ket prepare_control_state(const LCUExpansion& lcu) const;

    void apply_select_operator(qpp::ket& state,
                               const qpp::cmat& A,
                               const LCUExpansion& lcu,
                               const std::vector<qpp::idx>& control,
                               const std::vector<qpp::idx>& target) const;

    void unprepare_control_state(qpp::ket& state,
                                 const LCUExpansion& lcu) const;

private:
    double epsilon_;  // stored once, immutable
};

