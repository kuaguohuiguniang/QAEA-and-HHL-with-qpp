#pragma once
#include "interval_stepper.hpp"
#include <cmath>

namespace nozzle_qns {

    struct Diagnostics {
        idx outer_steps_done = 0;
        double t = 0.0;
        double max_abs_U = 0.0;
    };

    struct Solution {
        std::vector<UVec> U_final;
        Diagnostics diag{};
    };

    class QuantumODESolver {
        public:
            explicit QuantumODESolver(OneIntervalStepper stepper);

            Solution solve(std::vector<UVec> U0);

        private:
            OneIntervalStepper stepper_;
    };

} // namespace nozzle_qns
