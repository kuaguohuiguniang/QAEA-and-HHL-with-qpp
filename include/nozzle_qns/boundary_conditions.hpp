#pragma once
#include "quasi1d_driver.hpp"

namespace nozzle_qns {

    enum class ExitType {
        SupersonicOutflow, // extrapolate U1,U2,U3
        SubsonicOutflow    // extrapolate U1,U2; set U3 from fixed pe
    };

    struct BoundaryConfig {
        ExitType exit_type = ExitType::SupersonicOutflow;
        double pe = 0.0; // used if SubsonicOutflow
    };

    class BoundaryApplier {
        public:
            BoundaryApplier(Grid1D grid, NozzleArea area, GasModel gas, BoundaryConfig cfg);

            void apply(std::vector<UVec>& U) const;

        private:
            void apply_inlet(std::vector<UVec>& U) const;   // subsonic inflow
            void apply_outlet(std::vector<UVec>& U) const;  // super/subsonic outflow

            Grid1D grid_;
            NozzleArea area_;
            GasModel gas_;
            BoundaryConfig cfg_;
    };

} // namespace nozzle_qns
