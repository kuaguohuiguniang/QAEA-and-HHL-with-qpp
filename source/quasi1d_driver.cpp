#include "nozzle_qns/quasi1d_driver.hpp"
#include <utility>
#include <cmath>

namespace nozzle_qns {

    Primitive to_primitive(const UVec& U, double A, const GasModel& gas) {
        Primitive p;
        p.rho = U.U1 / A;
        p.v   = U.U2 / U.U1;
        double T = (gas.gamma - 1.0) *(U.U3 / U.U1 - 0.5 * gas.gamma * p.v * p.v);
        p.T = T;              // assuming R = 1
        p.p = p.rho * p.T;
        return p;
    }

    FluxVec flux_inviscid_quasi1d(const UVec& U, const GasModel& gas) {
        ensure(U.U1 > 0.0, "flux: U1 must be > 0");
        const double g = gas.gamma;
        FluxVec F;
        F.F1 = U.U2;
        F.F2 = (U.U2 * U.U2) / U.U1
            + (g - 1.0) / g * (U.U3 - (g / 2.0) * (U.U2 * U.U2) / U.U1);
        F.F3 = (g * U.U2 * U.U3) / U.U1
            - (g * (g - 1.0) / 2.0) * (U.U2 * U.U2 * U.U2) / (U.U1 * U.U1);
        return F;
    }


    SourceVec source_quasi1d(const UVec& U, double dlnA_dx, const GasModel& gas) {
        ensure(U.U1 > 0.0, "flux: U1 must be > 0");
        const double g = gas.gamma;
        SourceVec J;
        J.J1 = 0.0;
        J.J3 = 0.0;

        const double term = (g - 1.0) / g * (U.U3 - (g / 2.0) * (U.U2 * U.U2) / U.U1);
        J.J2 = term * dlnA_dx;

        return J;
    }


    Quasi1DInviscidDriver::Quasi1DInviscidDriver(Grid1D grid, NozzleArea area, GasModel gas)
        : grid_(std::move(grid)), area_(std::move(area)), gas_(std::move(gas)) {}

    const Grid1D& Quasi1DInviscidDriver::grid() const { return grid_; }
    const GasModel& Quasi1DInviscidDriver::gas() const { return gas_; }

    UVec Quasi1DInviscidDriver::dUdt(const std::vector<UVec>& U, idx j) const {
        ensure(U.size() == grid_.m, "Quasi1DInviscidDriver::dUdt: U size mismatch");
        ensure(j > 0 && j + 1 < grid_.m, "Quasi1DInviscidDriver::dUdt: j must be interior");

        const double Ajp = area_(grid_.x(j + 1));
        const double Ajm = area_(grid_.x(j - 1));

        const double dlnA_dx = (std::log(Ajp) - std::log(Ajm)) / (2.0 * grid_.dx);

        const FluxVec Fp = flux_inviscid_quasi1d(U[j + 1], gas_);
        const FluxVec Fm = flux_inviscid_quasi1d(U[j - 1], gas_);
        const SourceVec J  = source_quasi1d(U[j], dlnA_dx, gas_);

        UVec dU;
        dU.U1 = -(Fp.F1 - Fm.F1) / (2.0 * grid_.dx) + J.J1;
        dU.U2 = -(Fp.F2 - Fm.F2) / (2.0 * grid_.dx) + J.J2;
        dU.U3 = -(Fp.F3 - Fm.F3) / (2.0 * grid_.dx) + J.J3;
        return dU;
    }


    std::vector<UVec> Quasi1DInviscidDriver::dUdt_all(const std::vector<UVec>& U) const {
        ensure(U.size() == grid_.m, "Quasi1DInviscidDriver::dUdt_all: U size mismatch");
        std::vector<UVec> dUdt_vec(U.size());
        for (idx j = 1; j + 1 < U.size(); ++j) {
            dUdt_vec[j] = dUdt(U, j);
        }
        return dUdt_vec;
    }   
} // namespace nozzle_qns