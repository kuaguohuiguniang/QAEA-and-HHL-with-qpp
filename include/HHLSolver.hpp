#pragma once

#include "HHLinput.hpp"
#include <qpp/qpp.hpp> 
#include <vector>

using namespace qpp; 
using namespace qpp::literals;

class HHLSolver { 
    private: 
        int n_clock;    // Number of precision qubits for the eigenvalue estimation 
        double t0;      // The time evolution parameter 
        double C;       // Rotation constant (must be <= smallest eigenvalue to avoid arcsin error)
 
        // 1. Quantum Phase Estimation 
        // Applies Hadamard to clock, then Controlled-Unitary evolutions, then Inverse QFT 
        void apply_QPE(qpp::ket& state, const qpp::cmat& A, 
                        const std::vector<qpp::idx>& clock_regs, const std::vector<qpp::idx>& target_regs);
        
        // 2. Conditional Rotation 
        // Rotates the Ancilla qubit based on the integer state stored in Clock registers 
        void apply_rotation(qpp::ket& state, const std::vector<qpp::idx>& clock_regs, qpp::idx ancilla_reg);

        // 3. Inverse QPE 
        // The exact reverse of step 1 (Uncompute) 
        void apply_inv_QPE(qpp::ket& state, const qpp::cmat& A, 
                            const std::vector<qpp::idx>& clock_regs, const std::vector<qpp::idx>& target_regs);

    public: 
        // Constructor 
        HHLSolver(int precision, double evolution_time, double rotation_const); 
        
        // The Main Driver 
        // Orchestrates the registers and calls the 3 helpers 
        qpp::ket solve(const LinearSystem& sys);
};