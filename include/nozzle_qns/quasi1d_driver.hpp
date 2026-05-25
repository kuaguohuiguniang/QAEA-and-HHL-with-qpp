#pragma once
/**
 * @file quasi1d_driver.hpp
 * @brief Declares quasi-one-dimensional inviscid nozzle-flow equations.
 */

#include "nozzle_geometry.hpp"
#include <cmath>

namespace nozzle_qns {

    /**
     * @brief Converts conservative variables to primitive variables.
     *
     * @param U Conservative state.
     * @param A Local nozzle area.
     * @param gas Gas model parameters.
     * @return Primitive flow variables.
     */
    Primitive to_primitive(const UVec& U, double A, const GasModel& gas);

    /**
     * @brief Computes the inviscid quasi-1D flux F(U).
     */
    FluxVec flux_inviscid_quasi1d(const UVec& U, const GasModel& gas);

    /**
     * @brief Computes the geometric source term J(U, d log A / dx).
     */
    SourceVec source_quasi1d(const UVec& U, double dlnA_dx, const GasModel& gas);

    /**
     * @brief Central-difference spatial driver for the quasi-1D ODE system.
     */
    class Quasi1DInviscidDriver {
        public:
            /**
             * @brief Constructs a driver over a fixed grid, area model, and gas model.
             */
            Quasi1DInviscidDriver(Grid1D grid, NozzleArea area, GasModel gas);

            /**
             * @brief Returns the computational grid.
             */
            const Grid1D& grid() const;

            /**
             * @brief Returns the gas model.
             */
            const GasModel& gas() const;

            /**
             * @brief Computes dU/dt at one interior grid point.
             *
             * @param U Current conservative state over all grid points.
             * @param j Interior grid-point index.
             * @return Time derivative at grid point j.
             */
            UVec dUdt(const std::vector<UVec>& U, idx j) const;

            /**
             * @brief Computes dU/dt over the full grid.
             *
             * Boundary entries are left as zero because boundary conditions are
             * applied separately.
             */
            std::vector<UVec> dUdt_all(const std::vector<UVec>& U) const;

        private:
            Grid1D grid_;
            NozzleArea area_;
            GasModel gas_;
    };

} // namespace nozzle_qns
