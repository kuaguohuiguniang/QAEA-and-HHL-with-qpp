#include "nozzle_qns/nozzle_app.hpp"

namespace nozzle_qns {

    std::vector<UVec> NozzleRunner::initial_U_simple() const {
        const Grid1D grid = make_grid();
        const NozzleArea area = make_area();
        const double gamma = cfg_.gas.gamma;

        ensure(grid.m >= 3, "NozzleRunner::initial_U: grid must have >= 3 points");
        ensure(gamma > 1.0, "NozzleRunner::initial_U: gamma must be > 1");

        // Simple IC (dimensionless):
        // rho = 1, T = 1 everywhere; v = 0 everywhere.
        // Then:
        // U1 = rho*A = A
        // U2 = rho*A*v = 0
        // U3 = rho*A*(T/(gamma-1) + gamma/2*v^2) = A/(gamma-1)
        std::vector<UVec> U0(grid.m);
        for (idx j = 0; j < grid.m; ++j) {
            const double A = area(grid.x(j));
            ensure(A > 0.0, "NozzleRunner::initial_U: area must be > 0");

            U0[j].U1 = A;
            U0[j].U2 = 0.0;
            U0[j].U3 = A / (gamma - 1.0);
        }

        // Enforce BCs once so boundary cells are consistent with your chosen BC policy.
        BoundaryApplier bc(grid, area, cfg_.gas, cfg_.bc);
        bc.apply(U0);
        return U0;
    }

    std::vector<UVec> NozzleRunner::initial_U() const {
        const Grid1D grid = make_grid();
        const NozzleArea area = make_area();

        SteadyStateICConfig icfg;
        icfg.shock_present = false;
        icfg.seed = 1u;
        icfg.rho0 = 1.0;
        icfg.T0   = 1.0;

        auto U0 = make_initial_U_steady_state_SI5D(grid, area, cfg_.gas, icfg);

        make_bc().apply(U0);

        return U0;
    }

    Solution NozzleRunner::run(IMeanEstimator& estimator_backend,
                        const IKnotPolicy& knot_policy) const{
        // Build core problem objects from config
        Quasi1DInviscidDriver driver = make_driver();
        BoundaryApplier bc = make_bc();

        // Bridge from classical knot samples -> (optional quantum) mean estimation
        MeanIntegralComputer mean_integrals(estimator_backend);

        // Stepper and solver
        OneIntervalStepper stepper(cfg_.stepper_cfg,
                                std::move(driver),
                                std::move(bc),
                                std::move(mean_integrals),
                                knot_policy);

        QuantumODESolver solver(std::move(stepper));

        // Initial state and run
        std::vector<UVec> U0 = initial_U();
        return solver.solve(std::move(U0));
    }

}