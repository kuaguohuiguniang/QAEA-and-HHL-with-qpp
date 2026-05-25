#pragma once
/**
 * @file HHLSolver.hpp
 * @brief Declares a Quantum++ simulation of the HHL linear-system solver.
 */

#include "HHLinput.hpp"
#include <qpp/qpp.hpp> 
#include <vector>

using namespace qpp; 
using namespace qpp::literals;

/**
 * @brief Simulates the Harrow-Hassidim-Lloyd algorithm for a linear system.
 *
 * The solver expects a validated and normalized LinearSystem. It prepares the
 * input state, applies quantum phase estimation, performs the eigenvalue-
 * dependent ancilla rotation, uncomputes the phase register, and postselects on
 * the ancilla outcome associated with a successful solve.
 */
class HHLSolver { 
    private: 
        int n_clock;    ///< Number of clock qubits used for phase estimation.
        double t0;      ///< Hamiltonian evolution time used in exp(-iAt0).
        double C;       ///< Rotation constant; should not exceed the smallest relevant eigenvalue.
 
        /**
         * @brief Applies quantum phase estimation to encode eigenvalue phases.
         *
         * @param state State vector containing clock, target, and ancilla registers.
         * @param A Hermitian matrix defining the Hamiltonian simulation.
         * @param clock_regs Clock register qubit indices.
         * @param target_regs Target register qubit indices storing |b>.
         */
        void apply_QPE(qpp::ket& state, const qpp::cmat& A, 
                        const std::vector<qpp::idx>& clock_regs, const std::vector<qpp::idx>& target_regs);
        
        /**
         * @brief Rotates the ancilla according to the eigenvalue encoded in the clock register.
         *
         * @param state State vector after phase estimation.
         * @param clock_regs Clock register qubit indices.
         * @param ancilla_reg Ancilla qubit index.
         */
        void apply_rotation(qpp::ket& state, const std::vector<qpp::idx>& clock_regs, qpp::idx ancilla_reg);

        /**
         * @brief Applies inverse quantum phase estimation to uncompute the clock register.
         *
         * @param state State vector after the conditional ancilla rotation.
         * @param A Hermitian matrix defining the Hamiltonian simulation.
         * @param clock_regs Clock register qubit indices.
         * @param target_regs Target register qubit indices.
         */
        void apply_inv_QPE(qpp::ket& state, const qpp::cmat& A, 
                            const std::vector<qpp::idx>& clock_regs, const std::vector<qpp::idx>& target_regs);

    public: 
        /**
         * @brief Constructs an HHL solver with fixed simulation parameters.
         *
         * @param precision Number of clock qubits used for eigenvalue precision.
         * @param evolution_time Time parameter t0 for Hamiltonian evolution.
         * @param rotation_const Ancilla rotation constant C.
         */
        HHLSolver(int precision, double evolution_time, double rotation_const); 
        
        /**
         * @brief Runs the HHL circuit simulation and returns the postselected solution state.
         *
         * @param sys Validated linear system Ax=b.
         * @return Normalized solution state |x> if postselection succeeds; otherwise a zero ket.
         */
        qpp::ket solve(const LinearSystem& sys);
};
