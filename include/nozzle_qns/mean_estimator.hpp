#pragma once
/**
 * @file mean_estimator.hpp
 * @brief Defines the common interface for mean-estimation backends.
 */

#include <vector>

namespace nozzle_qns {

    /**
     * @brief Output of a backend mean-estimation call for normalized samples.
     */
    struct MeanEstimate {
        double G1{}, G2{}, G3{}; ///< Estimated means of normalized components in [0, 1].
        double eps{};            ///< Target or achieved additive error.
        double success_prob{};   ///< Backend-reported success probability.
    };

    /**
     * @brief Interface for classical or quantum mean-estimation backends.
     */
    class IMeanEstimator {
        public:
            virtual ~IMeanEstimator() = default;

            /**
             * @brief Estimates means of three normalized sample vectors.
             *
             * @param g1 First normalized component samples in [0, 1].
             * @param g2 Second normalized component samples in [0, 1].
             * @param g3 Third normalized component samples in [0, 1].
             * @param eps_target Requested additive error target.
             * @param delta_target Requested failure-probability target.
             * @return Component-wise mean estimate and diagnostics.
             */
            virtual MeanEstimate estimate_mean(
                const std::vector<double>& g1,
                const std::vector<double>& g2,
                const std::vector<double>& g3,
                double eps_target,
                double delta_target) = 0;
    };

} // namespace nozzle_qns
