#pragma once
/**
 * @file quantum_ode_solver.hpp
 * @brief Declares the high-level ODE solver wrapper for the nozzle simulation.
 */

#include "interval_stepper.hpp"
#include <cmath>

namespace nozzle_qns {

    /**
     * @brief Runtime diagnostics from a nozzle solve.
     */
    struct Diagnostics {
        idx outer_steps_done = 0; ///< Number of completed outer intervals.
        double t = 0.0;           ///< Final reported simulation time.
        double max_abs_U = 0.0;   ///< Maximum absolute conservative component.
    };

    /**
     * @brief Final state and diagnostics from a nozzle solve.
     */
    struct Solution {
        std::vector<UVec> U_final; ///< Conservative state at final time.
        Diagnostics diag{};        ///< Solver diagnostics.
    };

    /**
     * @brief Advances the nozzle state through all configured outer intervals.
     */
    class QuantumODESolver {
        public:
            /**
             * @brief Constructs a solver from a one-interval stepper.
             */
            explicit QuantumODESolver(OneIntervalStepper stepper);

            /**
             * @brief Solves the configured ODE problem from an initial state.
             *
             * @param U0 Initial conservative state.
             * @return Final state and diagnostics.
             */
            Solution solve(std::vector<UVec> U0);

        private:
            OneIntervalStepper stepper_;
    };

} // namespace nozzle_qns
