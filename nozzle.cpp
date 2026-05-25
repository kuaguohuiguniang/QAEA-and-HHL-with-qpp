/**
 * @file nozzle.cpp
 * @brief Example executable for comparing classical and QAEA nozzle backends.
 */

#include "nozzle_qns/nozzle_app.hpp"
#include "nozzle_qns/knot_policy.hpp"
#include "nozzle_qns/classical_mean_estimator.hpp"
#include "nozzle_qns/qpp_mean_estimator.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using nozzle_qns::idx;

namespace {

/**
 * @brief Prints representative primitive variables from a nozzle solution.
 */
void print_summary(const nozzle_qns::Solution& sol,
                   const nozzle_qns::Grid1D& grid,
                   const nozzle_qns::NozzleArea& area,
                   const nozzle_qns::GasModel& gas,
                   const std::string& tag) {
    std::cout << "\n==== " << tag << " ====\n";

    if (sol.U_final.empty()) {
        std::cout << "Solution U_final is empty.\n";
        return;
    }

    const idx m = static_cast<idx>(sol.U_final.size());
    const idx j_mid = m / 2;

    auto prim_mid = nozzle_qns::to_primitive(sol.U_final[j_mid], area(grid.x(j_mid)), gas);

    /// Mass-flow estimate: mdot = rho * v * A.
    const double mdot_mid = prim_mid.rho * prim_mid.v * area(grid.x(j_mid));

    std::cout << "Grid points m = " << m << "\n";
    std::cout << "Sample location x_mid = " << grid.x(j_mid) << "\n";
    std::cout << std::setprecision(10);
    std::cout << "Primitive@mid: rho=" << prim_mid.rho
              << ", v=" << prim_mid.v
              << ", T=" << prim_mid.T
              << ", p=" << prim_mid.p << "\n";
    std::cout << "mdot@mid = " << mdot_mid << "\n";

    /// Print inlet and outlet primitive variables as a rough sanity check.
    auto prim_in  = nozzle_qns::to_primitive(sol.U_final[0],       area(grid.x(0)),       gas);
    auto prim_out = nozzle_qns::to_primitive(sol.U_final[m - 1],   area(grid.x(m - 1)),   gas);

    std::cout << "Primitive@inlet : rho=" << prim_in.rho
              << ", v=" << prim_in.v
              << ", T=" << prim_in.T
              << ", p=" << prim_in.p << "\n";
    std::cout << "Primitive@outlet: rho=" << prim_out.rho
              << ", v=" << prim_out.v
              << ", T=" << prim_out.T
              << ", p=" << prim_out.p << "\n";

}

/**
 * @brief Builds a small demonstration configuration for the nozzle solver.
 */
nozzle_qns::NozzleConfig make_small_test_config() {
    nozzle_qns::NozzleConfig cfg;

    /// Domain and grid.
    cfg.x0 = 0.0;
    cfg.x1 = 3.0;
    cfg.Ngrid = 31; ///< Keep moderate for the demonstration run.

    /// Geometry: use the paper default nozzle.
    cfg.area_preset = nozzle_qns::NozzleAreaPreset::PaperDefault;

    /// Gas model.
    cfg.gas = nozzle_qns::GasModel{1.4};

    /// Boundary conditions for the shock-free case.
    cfg.bc.exit_type = nozzle_qns::ExitType::SupersonicOutflow;
    cfg.bc.pe = 0.6784; ///< Used only if SubsonicOutflow is selected.

    /// Time and algorithm parameters for a small demonstration run.
    cfg.stepper_cfg.tp.n = 4;        ///< Number of outer intervals.
    cfg.stepper_cfg.tp.k = 1;        ///< Hierarchy depth.
    cfg.stepper_cfg.tp.Tfinal = 0.2; ///< Final dimensionless time.

    cfg.stepper_cfg.knots_per_subsub = 8; ///< Number of mean-estimation knots.

    cfg.stepper_cfg.eps1 = 1e-2;
    cfg.stepper_cfg.delta = 1e-3;

    /// Taylor order in current code can only be 1.
    cfg.stepper_cfg.taylor.r = 1;

    return cfg;
}

} // namespace

/**
 * @brief Runs the classical and QAEA nozzle examples.
 */
int main() {
    try {
        /// Common configuration shared by both backends.
        nozzle_qns::NozzleConfig cfg = make_small_test_config();
        nozzle_qns::NozzleRunner runner(cfg);

        const auto grid = runner.make_grid();
        const auto area = runner.make_area();
        const auto gas  = cfg.gas;

        /// Knot policy: uniform knots in u in [0, 1].
        nozzle_qns::UniformKnotPolicy knot_policy;

        /// Classical backend baseline.

        {
            nozzle_qns::ClassicalMeanEstimator classical_backend;
            std::cout << "Running nozzle_qns with ClassicalMeanEstimator...\n";
            const nozzle_qns::Solution sol = runner.run(classical_backend, knot_policy);
            print_summary(sol, grid, area, gas, "Classical backend");
        }


        /// QppMeanEstimator backend using QAEA-style simulation.
        /// This can be slow depending on knots_per_subsub and eps.
        {
            nozzle_qns::QppMeanEstimator::Options opt;
            opt.index_qubits = 0;        ///< Auto-select from knot count.
            opt.default_eps = 1e-3;
            opt.default_delta = 1e-3;
            opt.seed = 12345;
            opt.verbose = false;
            opt.tag = "qaea";

            nozzle_qns::QppMeanEstimator qpp_backend(opt);

            std::cout << "\nRunning nozzle_qns with QppMeanEstimator (QAEA)...\n";
            const nozzle_qns::Solution sol = runner.run(qpp_backend, knot_policy);
            print_summary(sol, grid, area, gas, "Qpp/QAEA backend");
        }

        std::cout << "\nDone.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL] Exception: " << e.what() << "\n";
        return 1;
    }
}
