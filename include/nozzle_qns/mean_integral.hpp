#pragma once
/**
 * @file mean_integral.hpp
 * @brief Declares normalization and mean-integration helpers for knot samples.
 */

#include "mean_estimator.hpp"
#include <algorithm>
#include <cmath>

namespace nozzle_qns {

    /**
     * @brief Mean estimate for the three conservative derivative components.
     */
    struct IntegralEstimate {
        double mean_f1{}, mean_f2{}, mean_f3{}; ///< Estimated means of f over sampled knots.

        double f1_min{}, f1_max{}; ///< Normalization range for the first component.
        double f2_min{}, f2_max{}; ///< Normalization range for the second component.
        double f3_min{}, f3_max{}; ///< Normalization range for the third component.

        MeanEstimate mean_diag{}; ///< Diagnostics returned by the mean-estimation backend.
    };

    /**
     * @brief Converts raw derivative knot samples into backend mean estimates.
     */
    class MeanIntegralComputer {
        public:
            /**
             * @brief Constructs the computer around a mean-estimation backend.
             *
             * @param backend Backend used to estimate normalized means.
             */
            explicit MeanIntegralComputer(IMeanEstimator& backend) : backend_(backend) {}

            /**
             * @brief Estimates the mean of three raw derivative sample vectors.
             *
             * Raw samples are normalized to [0, 1], passed to the backend, then
             * mapped back to the original derivative scale.
             *
             * @param f1_knots Raw samples of the first derivative component.
             * @param f2_knots Raw samples of the second derivative component.
             * @param f3_knots Raw samples of the third derivative component.
             * @param eps_target Requested additive error target.
             * @param delta_target Requested failure-probability target.
             * @return Mean estimates and normalization diagnostics.
             */
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
