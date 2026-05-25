#pragma once
/**
 * @file classical_mean_estimator.hpp
 * @brief Declares the exact classical baseline mean-estimation backend.
 */

#include "nozzle_qns/mean_estimator.hpp"
#include <cstdint>
#include <random>
#include <vector>

namespace nozzle_qns {

    /**
     * @brief Classical backend for debugging and verification.
     *
     * By default this estimator returns the exact discrete mean of the knot
     * samples. It can optionally add Gaussian noise for stress testing.
     */
    class ClassicalMeanEstimator final : public IMeanEstimator {
        public:
            /**
             * @brief Configuration for the classical estimator.
             */
            struct Options {
                bool add_noise = false;       ///< Whether to add Gaussian noise.
                double gaussian_sigma = 0.0;  ///< Standard deviation of additive noise.
                std::uint64_t seed = 0;       ///< Random seed for optional noise.
            };

            /**
             * @brief Constructs a classical mean estimator.
             */
            explicit ClassicalMeanEstimator(Options opt = {false, 0.0, 0}) : opt_(opt) {
                rng_.seed(opt_.seed);
                dist_ = std::normal_distribution<double>(0.0, opt_.gaussian_sigma);
            }

            /**
             * @brief Returns exact component-wise discrete means of normalized samples.
             */
            MeanEstimate estimate_mean(const std::vector<double>& g1,
                                    const std::vector<double>& g2,
                                    const std::vector<double>& g3,
                                    double eps_target,
                                    double delta_target) override;

        private:
            Options opt_;
            std::mt19937_64 rng_{};
            std::normal_distribution<double> dist_{0.0, 0.0};

            /**
             * @brief Computes the discrete mean of one sample vector.
             */
            double mean_of(const std::vector<double>& g);
    };

} // namespace nozzle_qns
