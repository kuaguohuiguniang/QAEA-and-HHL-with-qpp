#include "nozzle_qns/knot_policy.hpp"
/**
 * @file knot_policy.cpp
 * @brief Implements knot-grid generation policies.
 */

namespace nozzle_qns {

    /**
     * @brief Creates uniformly spaced knots including endpoints.
     */
    KnotGrid UniformKnotPolicy::make_knots(idx count) const {
        ensure(count >= 2, "UniformKnotPolicy: need at least 2 knots");
        KnotGrid kg;
        kg.u.resize(count);
        for (idx i = 0; i < count; ++i) {
            kg.u[i] = double(i) / double(count - 1); ///< Uniform in [0, 1].
        }
        return kg;
    }

} // namespace nozzle_qns
