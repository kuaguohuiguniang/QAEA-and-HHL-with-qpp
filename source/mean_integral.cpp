#include "nozzle_qns/mean_integral.hpp"

#include <algorithm> 
#include <stdexcept>
#include <vector>
#include <cmath>

namespace nozzle_qns {

    static inline double clamp01(double x) {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }

    IntegralEstimate MeanIntegralComputer::estimate_mean_f_from_knots(
        const std::vector<double>& f1_knots,
        const std::vector<double>& f2_knots,
        const std::vector<double>& f3_knots,
        double eps_target,
        double delta_target) const
    {
        if (f1_knots.empty() || f2_knots.empty() || f3_knots.empty())
            throw std::runtime_error("MeanIntegralComputer: knot vectors must be non-empty");

        if (f1_knots.size() != f2_knots.size() || f1_knots.size() != f3_knots.size())
            throw std::runtime_error("MeanIntegralComputer: knot vectors must have the same size");

        IntegralEstimate result;

        // 1. min/max per component
        result.f1_min = *std::min_element(f1_knots.begin(), f1_knots.end());
        result.f1_max = *std::max_element(f1_knots.begin(), f1_knots.end());
        result.f2_min = *std::min_element(f2_knots.begin(), f2_knots.end());
        result.f2_max = *std::max_element(f2_knots.begin(), f2_knots.end());
        result.f3_min = *std::min_element(f3_knots.begin(), f3_knots.end());
        result.f3_max = *std::max_element(f3_knots.begin(), f3_knots.end());

        const double df1 = result.f1_max - result.f1_min;
        const double df2 = result.f2_max - result.f2_min;
        const double df3 = result.f3_max - result.f3_min;

        const bool c1 = (df1 == 0.0);
        const bool c2 = (df2 == 0.0);
        const bool c3 = (df3 == 0.0);

        std::vector<double> g1, g2, g3;
        g1.resize(f1_knots.size());
        g2.resize(f2_knots.size());
        g3.resize(f3_knots.size());

        for (std::size_t i = 0; i < f1_knots.size(); ++i) {
            g1[i] = c1 ? 0.0 : clamp01((f1_knots[i] - result.f1_min) / df1);
            g2[i] = c2 ? 0.0 : clamp01((f2_knots[i] - result.f2_min) / df2);
            g3[i] = c3 ? 0.0 : clamp01((f3_knots[i] - result.f3_min) / df3);
        }

        // 2. backend mean estimate (only needed if at least one component is non-constant)
        if (!c1 || !c2 || !c3) {
            result.mean_diag = backend_.estimate_mean(g1, g2, g3, eps_target, delta_target);
        } else {
            result.mean_diag = MeanEstimate{0.0, 0.0, 0.0, eps_target, 1.0};
        }

        // 3. rescale back
        result.mean_f1 = c1 ? result.f1_min : (result.f1_min + df1 * result.mean_diag.G1);
        result.mean_f2 = c2 ? result.f2_min : (result.f2_min + df2 * result.mean_diag.G2);
        result.mean_f3 = c3 ? result.f3_min : (result.f3_min + df3 * result.mean_diag.G3);

        // Ensure finite
        if (!std::isfinite(result.mean_f1) || !std::isfinite(result.mean_f2) || !std::isfinite(result.mean_f3))
            throw std::runtime_error("MeanIntegralComputer: produced non-finite mean_f");

        return result;
    }

} // namespace nozzle_qns
