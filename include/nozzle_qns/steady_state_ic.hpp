#pragma once
/**
 * @file steady_state_ic.hpp
 * @brief Declares steady-state initial-condition generation for nozzle flow.
 */

#include "nozzle_geometry.hpp"
#include <vector>

namespace nozzle_qns {

/**
 * @brief Configuration for steady-state perturbed initial conditions.
 */
struct SteadyStateICConfig {
    bool shock_present = false; ///< Whether to use shock-case perturbation magnitudes.

    unsigned seed = 1u; ///< Random seed for reproducible perturbations.

    double rho0 = 1.0; ///< Reference dimensionless density.
    double T0   = 1.0; ///< Reference dimensionless temperature.

    /**
     * @brief Whether to enforce constant mass flow after perturbation.
     *
     * The SI-5D setup perturbs mass flow globally; this option keeps the
     * conservative state consistent with that choice.
     */
    bool enforce_constant_mdot = true;
};

/**
 * @brief Builds a perturbed steady-state initial condition.
 *
 * The construction starts from an isentropic steady solution and adds random
 * shifts following the SI-5D style setup.
 *
 * @param grid Computational grid.
 * @param area Nozzle area model.
 * @param gas Gas model parameters.
 * @param cfg Initial-condition configuration.
 * @return Conservative state values U[j] for all grid points.
 */
std::vector<UVec> make_initial_U_steady_state_SI5D(
    const Grid1D& grid,
    const NozzleArea& area,
    const GasModel& gas,
    const SteadyStateICConfig& cfg);

} // namespace nozzle_qns
