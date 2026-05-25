#pragma once
/**
 * @file qpp_mean_estimator.hpp
 * @brief Declares a Quantum++-backed QAEA-style mean-estimation backend.
 *
 * Inputs are normalized component samples in [0, 1]. The implementation in
 * source/qpp_mean_estimator.cpp estimates each component mean with an
 * amplitude-estimation-style state-vector simulation.
 */

#include "nozzle_qns/mean_estimator.hpp"
#include "nozzle_qns/nozzle_types.hpp"
#include <qpp/qpp.hpp>

#include <cstdint>
#include <random>
#include <string>

namespace nozzle_qns {
    using namespace qpp;

    /**
     * @brief Quantum++ simulation backend for Quantum Amplitude Estimation.
     */
    class QppMeanEstimator final : public IMeanEstimator {
        public:
            /**
             * @brief Configuration for the QAEA simulation.
             */
            struct Options {
                idx index_qubits = 0; ///< Number of index qubits; zero selects automatically.

                double default_eps = 1e-3;   ///< Default additive error target.
                double default_delta = 1e-3; ///< Default failure-probability target.

                idx max_iters = 0; ///< Reserved iteration cap for future variants.

                std::uint64_t seed = 0; ///< Random seed for measurement sampling.

                bool verbose = false; ///< Enables backend logging when used.
                std::string tag;      ///< Optional label for log messages.
            };

            /**
             * @brief Constructs a QAEA mean estimator.
             */
            explicit QppMeanEstimator(Options opt) : opt_(std::move(opt)) {
                rng_.seed(opt_.seed);
            }

            /**
             * @brief Estimates means of three normalized component sample vectors.
             */
            MeanEstimate estimate_mean(const std::vector<double>& g1,
                                    const std::vector<double>& g2,
                                    const std::vector<double>& g3,
                                    double eps_target,
                                    double delta_target) override;

            /**
             * @brief Returns the backend configuration.
             */
            const Options& options() const { return opt_; }

        private:
            Options opt_;
            std::mt19937_64 rng_{};
    };

} // namespace nozzle_qns
