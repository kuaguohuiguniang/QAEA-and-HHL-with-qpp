#include "nozzle_qns/boundary_conditions.hpp"
/**
 * @file boundary_conditions.cpp
 * @brief Implements nozzle inlet and outlet boundary conditions.
 */

#include <utility>

namespace nozzle_qns {

    /**
     * @brief Constructs a boundary-condition applier.
     */
    BoundaryApplier::BoundaryApplier(Grid1D grid, NozzleArea area, GasModel gas, BoundaryConfig cfg)
        : grid_(grid), area_(std::move(area)), gas_(gas), cfg_(cfg) {
        ensure(grid_.m >= 3, "BoundaryApplier: grid must have >= 3 points");
        if (cfg_.exit_type == ExitType::SubsonicOutflow) {
            ensure(cfg_.pe > 0.0, "BoundaryApplier: pe must be > 0 for subsonic outflow");
        }
    }

    /**
     * @brief Applies the subsonic inlet condition.
     */
    void BoundaryApplier::apply_inlet(std::vector<UVec>& U) const {
        ensure(U.size() == grid_.m, "BoundaryApplier: U size mismatch");

        /// Subsonic inflow at inlet j = 0: rho = 1 and T = 1.
        U[0].U1 = area_(grid_.x(0));       ///< rho * A = A at the inlet.
        U[0].U2 = 2.0 * U[1].U2 - U[2].U2; ///< Extrapolate mass flow.

        const double gamma = gas_.gamma;
        const double v0 = U[0].U2 / U[0].U1;

        U[0].U3 = U[0].U1 * (1.0 / (gamma - 1.0) + (gamma / 2.0) * v0 * v0);
    }

    /**
     * @brief Applies either supersonic or subsonic outlet conditions.
     */
    void BoundaryApplier::apply_outlet(std::vector<UVec>& U) const {
        ensure(U.size() == grid_.m, "BoundaryApplier: U size mismatch");

        const idx j = grid_.m - 1;

        if (cfg_.exit_type == ExitType::SupersonicOutflow) {
            U[j].U1 = 2.0 * U[j - 1].U1 - U[j - 2].U1;
            U[j].U2 = 2.0 * U[j - 1].U2 - U[j - 2].U2;
            U[j].U3 = 2.0 * U[j - 1].U3 - U[j - 2].U3;
            return;
        }

        /// Subsonic outflow: extrapolate U1 and U2; set U3 from fixed pe.
        U[j].U1 = 2.0 * U[j - 1].U1 - U[j - 2].U1;
        U[j].U2 = 2.0 * U[j - 1].U2 - U[j - 2].U2;

        const double Ae = area_(grid_.x(j));
        const double gamma = gas_.gamma;

        U[j].U3 = (cfg_.pe * Ae) / (gamma - 1.0)
                + (gamma / 2.0) * (U[j].U2 * U[j].U2) / U[j].U1;
    }

    /**
     * @brief Applies inlet and outlet boundary conditions.
     */
    void BoundaryApplier::apply(std::vector<UVec>& U) const {
        apply_inlet(U);
        apply_outlet(U);
    }

} // namespace nozzle_qns
