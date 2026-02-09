#pragma once
#include <vector>

namespace nozzle_qns {

    // Output of the *backend* mean estimation for normalized g in [0,1]
    struct MeanEstimate {
        double G1{}, G2{}, G3{};       // mean of normalized g in [0,1]
        double eps{};                  // achieved/target error
        double success_prob{};         // estimator’s own success probability
    };

    // Quantum (qpp/QAEA) or classical (baseline) backend.
    class IMeanEstimator {
        public:
            virtual ~IMeanEstimator() = default;

            // Inputs: normalized samples g_a(u_k) in [0,1] on knots u_k.
            // Output: estimate of mean(g_a) over the knot set.
            virtual MeanEstimate estimate_mean(
                const std::vector<double>& g1,
                const std::vector<double>& g2,
                const std::vector<double>& g3,
                double eps_target,
                double delta_target) = 0;
    };

} // namespace nozzle_qns