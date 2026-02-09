#pragma once
#include "nozzle_geometry.hpp" // Grid1D, NozzleArea, GasModel, UVec, idx
#include <vector>

namespace nozzle_qns {

struct SteadyStateICConfig {
    // Perturbation magnitudes differ if a shock is present.
    bool shock_present = false;

    // RNG seed for reproducibility
    unsigned seed = 1u;

    // Dimensionless variables.
    double rho0 = 1.0;
    double T0   = 1.0;

    // If true, enforce U2 (mass flow) to be constant across grid
    // after perturbation. SI-5D perturbs mdot globally; this keeps consistency.
    bool enforce_constant_mdot = true;
};

// Builds initial condition
// Start from exact isentropic steady solution, then add random shifts (SI-5D).
// Returns U[j] for j=0..m-1 (size m).
std::vector<UVec> make_initial_U_steady_state_SI5D(
    const Grid1D& grid,
    const NozzleArea& area,
    const GasModel& gas,
    const SteadyStateICConfig& cfg);

} // namespace nozzle_qns
