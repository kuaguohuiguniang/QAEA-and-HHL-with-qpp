#pragma once
#include "mean_estimator.hpp"
#include <algorithm>
#include <cmath>

namespace nozzle_qns {

    struct IntegralEstimate {
        // estimated mean of f over u in [0,1] (discrete knots)
        double mean_f1{}, mean_f2{}, mean_f3{};

        // diagnostic min/max used for normalization
        double f1_min{}, f1_max{};
        double f2_min{}, f2_max{};
        double f3_min{}, f3_max{};

        // backend diagnostics
        MeanEstimate mean_diag{};
    };

    class MeanIntegralComputer {
        public:
            explicit MeanIntegralComputer(IMeanEstimator& backend) : backend_(backend) {}

            // Input: knot samples of f (not g). Output: estimated mean of f.
            IntegralEstimate estimate_mean_f_from_knots(
                const std::vector<double>& f1_knots,
                const std::vector<double>& f2_knots,
                const std::vector<double>& f3_knots,
                double eps_target,
                double delta_target) const;

        private:
            IMeanEstimator& backend_;
    };

} // namespace nozzle_qns
