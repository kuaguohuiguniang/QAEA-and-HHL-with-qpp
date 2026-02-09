#include "nozzle_qns/quantum_ode_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nozzle_qns {

QuantumODESolver::QuantumODESolver(OneIntervalStepper stepper)
    : stepper_(std::move(stepper)) {
    const auto& cfg = stepper_.config();
    ensure(cfg.tp.n >= 1, "QuantumODESolver: tp.n must be >= 1");
    ensure(cfg.tp.k >= 1, "QuantumODESolver: tp.k must be >= 1");
    ensure(cfg.tp.Tfinal > 0.0, "QuantumODESolver: tp.Tfinal must be > 0");
}

Solution QuantumODESolver::solve(std::vector<UVec> U0) {
    // Basic sanity
    if (U0.empty())
        throw std::runtime_error("QuantumODESolver::solve: empty initial state");

    Solution sol;
    sol.U_final = std::move(U0);

    const auto& cfg = stepper_.config();
    const double h = cfg.tp.h();

    sol.diag.outer_steps_done = 0;
    sol.diag.t = 0.0;
    sol.diag.max_abs_U = 0.0;

    auto update_max_abs = [&](const std::vector<UVec>& U) {
        double mx = 0.0;
        for (const auto& u : U) {
            mx = std::max(mx, std::abs(u.U1));
            mx = std::max(mx, std::abs(u.U2));
            mx = std::max(mx, std::abs(u.U3));
        }
        sol.diag.max_abs_U = mx;
    };

    update_max_abs(sol.U_final);

    // Outer time stepping: n intervals
    for (idx i = 0; i < cfg.tp.n; ++i) {
        stepper_.advance_one_interval(sol.U_final);

        sol.diag.outer_steps_done = i + 1;
        sol.diag.t += h;
        update_max_abs(sol.U_final);
    }

    sol.diag.t = cfg.tp.Tfinal;

    return sol;
}

} // namespace nozzle_qns
