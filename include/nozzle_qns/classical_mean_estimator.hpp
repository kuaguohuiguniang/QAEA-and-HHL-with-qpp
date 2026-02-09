#pragma once

#include "nozzle_qns/mean_estimator.hpp"
#include <cstdint>
#include <random>
#include <vector>

namespace nozzle_qns {

    // Classical backend for debugging/verification.
    // Default: returns exact discrete mean of knot samples.
    // Optionally: add small Gaussian noise for stress-testing.
    class ClassicalMeanEstimator final : public IMeanEstimator {
        public:
            struct Options {
                bool add_noise = false;
                double gaussian_sigma = 0.0; // stddev of additive noise
                std::uint64_t seed = 0;
            };

            explicit ClassicalMeanEstimator(Options opt = {false, 0.0, 0}) : opt_(opt) {
                rng_.seed(opt_.seed);
                dist_ = std::normal_distribution<double>(0.0, opt_.gaussian_sigma);
            }

            MeanEstimate estimate_mean(const std::vector<double>& g1,
                                    const std::vector<double>& g2,
                                    const std::vector<double>& g3,
                                    double eps_target,
                                    double delta_target) override;

        private:
            Options opt_;
            std::mt19937_64 rng_{};
            std::normal_distribution<double> dist_{0.0, 0.0};

            double mean_of(const std::vector<double>& g);
    };

} // namespace nozzle_qns
