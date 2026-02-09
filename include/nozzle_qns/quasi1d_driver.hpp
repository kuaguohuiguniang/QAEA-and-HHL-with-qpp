#pragma once
#include "nozzle_geometry.hpp"
#include <cmath>

namespace nozzle_qns {

    // U -> primitive
    Primitive to_primitive(const UVec& U, double A, const GasModel& gas);

    // Inviscid flux F(U)
    FluxVec flux_inviscid_quasi1d(const UVec& U, const GasModel& gas);

    // Source J(U) with d/dx ln A
    SourceVec source_quasi1d(const UVec& U, double dlnA_dx, const GasModel& gas);

    // Central difference driver f(U)
    class Quasi1DInviscidDriver {
        public:
            Quasi1DInviscidDriver(Grid1D grid, NozzleArea area, GasModel gas);

            const Grid1D& grid() const;
            const GasModel& gas() const;

            // interior derivative only (j in 1..m-2)
            UVec dUdt(const std::vector<UVec>& U, idx j) const;

            // interior-only; boundaries left as zeros
            std::vector<UVec> dUdt_all(const std::vector<UVec>& U) const;

        private:
            Grid1D grid_;
            NozzleArea area_;
            GasModel gas_;
    };

} // namespace nozzle_qns
