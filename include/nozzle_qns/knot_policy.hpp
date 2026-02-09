#pragma once
#include "nozzle_types.hpp"
#include <vector>

namespace nozzle_qns {

    struct KnotGrid {
        std::vector<double> u; // u in [0,1]
    };

    class IKnotPolicy {
        public:
            virtual ~IKnotPolicy() = default;
            virtual KnotGrid make_knots(idx count) const = 0;
    };

    // Simple default: uniform knots including endpoints.
    class UniformKnotPolicy final : public IKnotPolicy {
        public:
            KnotGrid make_knots(idx count) const override;
    };

} // namespace nozzle_qns
