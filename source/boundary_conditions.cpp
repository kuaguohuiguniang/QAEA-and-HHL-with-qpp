#include "nozzle_qns/boundary_conditions.hpp"
#include <utility>

namespace nozzle_qns {

    BoundaryApplier::BoundaryApplier(Grid1D grid, NozzleArea area, GasModel gas, BoundaryConfig cfg)
        : grid_(grid), area_(std::move(area)), gas_(gas), cfg_(cfg) {
        ensure(grid_.m >= 3, "BoundaryApplier: grid must have >= 3 points");
        if (cfg_.exit_type == ExitType::SubsonicOutflow) {
            ensure(cfg_.pe > 0.0, "BoundaryApplier: pe must be > 0 for subsonic outflow");
        }
    }

    void BoundaryApplier::apply_inlet(std::vector<UVec>& U) const {
        ensure(U.size() == grid_.m, "BoundaryApplier: U size mismatch");

        // Subsonic inflow at inlet (j=0): rho=1, T=1
        U[0].U1 = area_(grid_.x(0));           // rho*A = A
        U[0].U2 = 2.0 * U[1].U2 - U[2].U2;     // extrapolate mass flow

        const double gamma = gas_.gamma;
        const double v0 = U[0].U2 / U[0].U1;

        U[0].U3 = U[0].U1 * (1.0 / (gamma - 1.0) + (gamma / 2.0) * v0 * v0);
    }

    void BoundaryApplier::apply_outlet(std::vector<UVec>& U) const {
        ensure(U.size() == grid_.m, "BoundaryApplier: U size mismatch");

        const idx j = grid_.m - 1;

        if (cfg_.exit_type == ExitType::SupersonicOutflow) {
            U[j].U1 = 2.0 * U[j - 1].U1 - U[j - 2].U1;
            U[j].U2 = 2.0 * U[j - 1].U2 - U[j - 2].U2;
            U[j].U3 = 2.0 * U[j - 1].U3 - U[j - 2].U3;
            return;
        }

        // Subsonic outflow: extrapolate U1,U2; set U3 from fixed pe
        U[j].U1 = 2.0 * U[j - 1].U1 - U[j - 2].U1;
        U[j].U2 = 2.0 * U[j - 1].U2 - U[j - 2].U2;

        const double Ae = area_(grid_.x(j));
        const double gamma = gas_.gamma;

        U[j].U3 = (cfg_.pe * Ae) / (gamma - 1.0)
                + (gamma / 2.0) * (U[j].U2 * U[j].U2) / U[j].U1;
    }

    void BoundaryApplier::apply(std::vector<UVec>& U) const {
        apply_inlet(U);
        apply_outlet(U);
    }

} // namespace nozzle_qns
