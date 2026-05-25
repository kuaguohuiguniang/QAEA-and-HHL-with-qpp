#include "nozzle_qns/classical_mean_estimator.hpp"
/**
 * @file classical_mean_estimator.cpp
 * @brief Implements the exact classical mean-estimation backend.
 */

#include <numeric>      
#include <stdexcept>   

namespace nozzle_qns {

    /**
     * @brief Computes exact discrete means for three normalized sample vectors.
     */
    MeanEstimate ClassicalMeanEstimator::estimate_mean(const std::vector<double>& g1,
                                                        const std::vector<double>& g2,
                                                        const std::vector<double>& g3,
                                                        double eps_target,
                                                        double delta_target) {
        if (g1.empty() || g2.empty() || g3.empty())
            throw std::runtime_error("ClassicalMeanEstimator: empty input vector");

        if (g1.size() != g2.size() || g1.size() != g3.size())
            throw std::runtime_error("ClassicalMeanEstimator: size mismatch");

        /// Compute exact discrete means of the knot samples.
        double mean_g1 = mean_of(g1);
        double mean_g2 = mean_of(g2);
        double mean_g3 = mean_of(g3);

        return MeanEstimate{mean_g1, mean_g2, mean_g3, eps_target, std::max(1.0, delta_target)};
    }

    /**
     * @brief Computes the arithmetic mean of one vector.
     */
    double ClassicalMeanEstimator::mean_of(const std::vector<double>& g) {
        if (g.empty())
            throw std::runtime_error("ClassicalMeanEstimator: mean_of on empty vector");

        return std::accumulate(g.begin(), g.end(), 0.0) / g.size();
    }

} // namespace nozzle_qns
