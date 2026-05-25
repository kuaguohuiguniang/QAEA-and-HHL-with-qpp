#pragma once
/**
 * @file boundary_conditions.hpp
 * @brief Declares inlet and outlet boundary-condition handling for nozzle flow.
 */

#include "quasi1d_driver.hpp"

namespace nozzle_qns {

    /**
     * @brief Outlet boundary-condition mode.
     */
    enum class ExitType {
        SupersonicOutflow, ///< Extrapolate U1, U2, and U3.
        SubsonicOutflow    ///< Extrapolate U1 and U2; set U3 from fixed exit pressure.
    };

    /**
     * @brief Boundary-condition configuration for the nozzle solver.
     */
    struct BoundaryConfig {
        ExitType exit_type = ExitType::SupersonicOutflow; ///< Outlet boundary type.
        double pe = 0.0; ///< Fixed exit pressure used for subsonic outflow.
    };

    /**
     * @brief Applies inlet and outlet boundary conditions to a nozzle state.
     */
    class BoundaryApplier {
        public:
            /**
             * @brief Constructs a boundary-condition applier for a fixed grid and geometry.
             *
             * @param grid Computational grid.
             * @param area Nozzle area model.
             * @param gas Gas model parameters.
             * @param cfg Boundary-condition configuration.
             */
            BoundaryApplier(Grid1D grid, NozzleArea area, GasModel gas, BoundaryConfig cfg);

            /**
             * @brief Applies inlet and outlet boundary conditions in place.
             *
             * @param U Conservative state vector over the grid.
             */
            void apply(std::vector<UVec>& U) const;

        private:
            /**
             * @brief Applies the subsonic inlet condition.
             */
            void apply_inlet(std::vector<UVec>& U) const;

            /**
             * @brief Applies the configured outlet condition.
             */
            void apply_outlet(std::vector<UVec>& U) const;

            Grid1D grid_;
            NozzleArea area_;
            GasModel gas_;
            BoundaryConfig cfg_;
    };

} // namespace nozzle_qns
