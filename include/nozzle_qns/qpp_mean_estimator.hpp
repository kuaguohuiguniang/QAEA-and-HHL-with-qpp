// QPP-backed estimator placeholder for Quantum Amplitude Estimation (QAEA/QAE).
// This header defines the interface and configuration knobs; you implement the
// actual quantum circuit/simulation logic in source/qpp_mean_estimator.cpp.
//
// Design intent:
// - Inputs g1,g2,g3 are normalized samples in [0,1] at knots.
// - This backend should estimate the mean of each list (as an expectation value)
//   using an amplitude-estimation style routine (or a simulator equivalent).

#pragma once

#include "nozzle_qns/mean_estimator.hpp"
#include "nozzle_qns/nozzle_types.hpp"
#include <qpp/qpp.hpp>

#include <cstdint>
#include <random>
#include <string>

namespace nozzle_qns {
    using namespace qpp;

    class QppMeanEstimator final : public IMeanEstimator {
        public:
            struct Options {
                idx index_qubits = 0;

                // “Precision knobs” for your QAE implementation.
                double default_eps = 1e-3;
                double default_delta = 1e-3;

                idx max_iters = 0;

                // Randomness control (for measurement sampling / stochastic variants).
                std::uint64_t seed = 0;

                // Logging / debugging
                bool verbose = false;
                std::string tag; // optional label for logs
            };

            explicit QppMeanEstimator(Options opt) : opt_(std::move(opt)) {
                rng_.seed(opt_.seed);
            }

            MeanEstimate estimate_mean(const std::vector<double>& g1,
                                    const std::vector<double>& g2,
                                    const std::vector<double>& g3,
                                    double eps_target,
                                    double delta_target) override;

            const Options& options() const { return opt_; }

        private:
            Options opt_;
            std::mt19937_64 rng_{};
    };

} // namespace nozzle_qns
