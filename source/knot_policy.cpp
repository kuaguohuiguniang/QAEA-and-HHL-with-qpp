#include "nozzle_qns/knot_policy.hpp"

namespace nozzle_qns {

    KnotGrid UniformKnotPolicy::make_knots(idx count) const {
        ensure(count >= 2, "UniformKnotPolicy: need at least 2 knots");
        KnotGrid kg;
        kg.u.resize(count);
        for (idx i = 0; i < count; ++i) {
            kg.u[i] = double(i) / double(count - 1); // uniform in [0,1]
        }
        return kg;
    }

} // namespace nozzle_qns