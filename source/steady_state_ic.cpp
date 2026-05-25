#include "nozzle_qns/steady_state_ic.hpp"
/**
 * @file steady_state_ic.cpp
 * @brief Implements perturbed isentropic steady-state initial conditions.
 */

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace nozzle_qns {

/**
 * @brief Computes A/A* as a function of Mach number for isentropic flow.
 */
static double area_mach_ratio(double M, double gamma) {
    const double g = gamma;
    const double term = (2.0 / (g + 1.0)) * (1.0 + (g - 1.0) * 0.5 * M * M);
    const double expo = (g + 1.0) / (g - 1.0);
    const double val2 = (1.0 / (M * M)) * std::pow(term, expo);
    return std::sqrt(val2);
}

enum class MachBranch { Subsonic, Supersonic };

/**
 * @brief Solves for Mach number from an area ratio on a selected branch.
 */
static double solve_mach_from_area_ratio(double A_over_Astar, double gamma, MachBranch branch) {
    /// Robust bisection on monotone subsonic and supersonic branches.
    if (!(A_over_Astar >= 1.0)) {
        A_over_Astar = 1.0;
    }

    auto f = [&](double M) {
        return area_mach_ratio(M, gamma) - A_over_Astar;
    };

    const double eps = 1e-12;

    double lo, hi;
    if (branch == MachBranch::Subsonic) {
        lo = 1e-10;
        hi = 1.0 - 1e-10;
    } else {
        lo = 1.0 + 1e-10;
        hi = 50.0; ///< Large enough for typical nozzle ratios.
    }

    /// Ensure the chosen interval brackets a root.
    double flo = f(lo);
    double fhi = f(hi);

    /// Expand the supersonic upper bound if needed.
    if (branch == MachBranch::Supersonic) {
        int expand = 0;
        while (fhi < 0.0 && expand < 30) {
            hi *= 1.5;
            fhi = f(hi);
            ++expand;
        }
    }

    /// If still not bracketed, return a near-sonic fallback.
    if (flo * fhi > 0.0) {
        return (branch == MachBranch::Subsonic) ? 0.999999 : 1.000001;
    }

    /// Bisection solve.
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double fmid = f(mid);
        if (std::abs(fmid) < 1e-13 || (hi - lo) < eps) {
            return mid;
        }
        if (flo * fmid > 0.0) {
            lo = mid;
            flo = fmid;
        } else {
            hi = mid;
            fhi = fmid;
        }
    }
    return 0.5 * (lo + hi);
}

/**
 * @brief Draws a uniform random number from [a, b].
 */
static inline double uni(std::mt19937& rng, double a, double b) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

/**
 * @brief Builds an SI-5D-style perturbed steady-state initial condition.
 */
std::vector<UVec> make_initial_U_steady_state_SI5D(
    const Grid1D& grid,
    const NozzleArea& area,
    const GasModel& gas,
    const SteadyStateICConfig& cfg)
{
    const idx m = grid.m;
    if (m < 3) throw std::runtime_error("make_initial_U_steady_state_SI5D: grid.m must be >= 3");
    if (!(gas.gamma > 1.0)) throw std::runtime_error("make_initial_U_steady_state_SI5D: gamma must be > 1");

    /// Evaluate areas and find the throat A* = min A.
    std::vector<double> A(m);
    idx j_throat = 0;
    double Astar = std::numeric_limits<double>::infinity();
    for (idx j = 0; j < m; ++j) {
        const double Aj = area(grid.x(j));
        if (!(Aj > 0.0)) throw std::runtime_error("make_initial_U_steady_state_SI5D: area must be > 0");
        A[j] = Aj;
        if (Aj < Astar) {
            Astar = Aj;
            j_throat = j;
        }
    }

    /// Steady isentropic solution in dimensionless variables.
    std::vector<double> M(m), Tss(m), rhoss(m);

    const double g = gas.gamma;
    for (idx j = 0; j < m; ++j) {
        if (j == j_throat) {
            M[j] = 1.0;
        } else if (j < j_throat) {
            M[j] = solve_mach_from_area_ratio(A[j] / Astar, g, MachBranch::Subsonic);
        } else {
            M[j] = solve_mach_from_area_ratio(A[j] / Astar, g, MachBranch::Supersonic);
        }

        const double denom = 1.0 + (g - 1.0) * 0.5 * M[j] * M[j];

        Tss[j] = cfg.T0 * (1.0 / denom);
        
        rhoss[j] = cfg.rho0 * std::pow(denom, -1.0 / (g - 1.0));
    }

    const double Tssmin  = *std::min_element(Tss.begin(), Tss.end());
    const double rhossmin = *std::min_element(rhoss.begin(), rhoss.end());

    /// Compute steady mass flow rate from the throat.
    const double Tt = Tss[j_throat];
    const double rhot = rhoss[j_throat];
    const double vt = std::sqrt(Tt); ///< Matches the SI-5D nondimensional sound speed.
    const double mdot_ss = rhot * A[j_throat] * vt; ///< M = 1 at the throat.

    /// Random shifts around the steady solution.
    const double dTmax   = (cfg.shock_present ? 0.01 : 0.02) * Tssmin;
    const double drhomax = (cfg.shock_present ? 0.02 : 0.10) * rhossmin;
    const double dmdotmax = 0.01 * mdot_ss;

    std::mt19937 rng(cfg.seed);

    std::vector<double> Tinit(m), rhoinit(m);
    for (idx j = 0; j < m; ++j) {
        const double dT = uni(rng, -dTmax, dTmax);
        const double dr = uni(rng, -drhomax, drhomax);
        Tinit[j] = Tss[j] + dT;
        rhoinit[j] = rhoss[j] + dr;

        if (Tinit[j] <= 0.0) Tinit[j] = 1e-8;
        if (rhoinit[j] <= 0.0) rhoinit[j] = 1e-8;
    }

    const double dmdot = uni(rng, -dmdotmax, dmdotmax);
    const double mdot_init = mdot_ss + dmdot;

    /// Build conservative U from rho_init, T_init, and mdot_init.
    std::vector<UVec> U(m);
    for (idx j = 0; j < m; ++j) {
        const double v = mdot_init / (rhoinit[j] * A[j]); ///< From mdot = rho * A * v.

        U[j].U1 = rhoinit[j] * A[j];
        U[j].U2 = rhoinit[j] * A[j] * v; ///< Equals mdot_init up to rounding.

        /// U3 = rho * A * (T / (g - 1) + g / 2 * v^2).
        U[j].U3 = rhoinit[j] * A[j] * (Tinit[j] / (g - 1.0) + (g * 0.5) * v * v);
    }

    /// Enforce U2 as the constant mass flow mdot_init.
    if (cfg.enforce_constant_mdot) {
        for (idx j = 0; j < m; ++j) {
            U[j].U2 = mdot_init;
        }
        /// Update U3 consistently with adjusted U2 by recomputing v.
        for (idx j = 0; j < m; ++j) {
            const double rho = U[j].U1 / A[j];
            const double v = U[j].U2 / U[j].U1;
            const double T = Tinit[j];
            U[j].U3 = rho * A[j] * (T / (g - 1.0) + (g * 0.5) * v * v);
        }
    }

    return U;
}

} // namespace nozzle_qns
