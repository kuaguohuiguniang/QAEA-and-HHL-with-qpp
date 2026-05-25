#pragma once
/**
 * @file nozzle_app.hpp
 * @brief Declares high-level configuration and orchestration for nozzle runs.
 */

#include "quantum_ode_solver.hpp"
#include "interval_stepper.hpp"
#include "mean_integral.hpp"
#include "knot_policy.hpp"
#include "boundary_conditions.hpp"
#include "quasi1d_driver.hpp"
#include "nozzle_geometry.hpp"
#include "nozzle_types.hpp"
#include "steady_state_ic.hpp"

#include <utility>

namespace nozzle_qns {

    /**
     * @brief Selects the nozzle area model.
     */
    enum class NozzleAreaPreset {
        PaperDefault, ///< A(x) = 1 + 2.2 * (x - 1.5)^2 on [0, 3].
        Custom        ///< Use user-provided NozzleArea.
    };

    /**
     * @brief Complete configuration for a nozzle simulation run.
     */
    struct NozzleConfig {
        double x0 = 0.0; ///< Left domain boundary.
        double x1 = 3.0; ///< Right domain boundary.
        idx Ngrid = 101; ///< Number of grid points.

        NozzleAreaPreset area_preset = NozzleAreaPreset::PaperDefault; ///< Built-in or custom area selector.
        NozzleArea custom_area = NozzleArea([](double) { return 1.0; }); ///< Area model used when area_preset is Custom.

        GasModel gas{1.4}; ///< Gas model parameters.

        /**
         * @brief Boundary-condition configuration.
         *
         * The inlet is a subsonic inflow with rho = 1 and T = 1 in BoundaryApplier.
         */
        BoundaryConfig bc{
            ExitType::SupersonicOutflow,
            0.6784
        };

        StepperConfig stepper_cfg{}; ///< ODE and mean-estimation control parameters.
    };

    /**
     * @brief High-level runner that wires together the nozzle solver components.
     */
    class NozzleRunner {
        public:
            /**
             * @brief Constructs a runner from a complete nozzle configuration.
             */
            explicit NozzleRunner(NozzleConfig cfg) : cfg_(std::move(cfg)) {
                ensure(cfg_.Ngrid >= 3, "NozzleRunner: Ngrid must be >= 3");
                ensure(cfg_.x1 > cfg_.x0, "NozzleRunner: x1 must be > x0");
            }

            /**
             * @brief Builds and runs the full nozzle solve pipeline.
             *
             * @param estimator_backend Classical or QAEA-style mean estimator.
             * @param knot_policy Policy used to choose normalized sampling knots.
             * @return Final solution and diagnostics.
             */
            Solution run(IMeanEstimator& estimator_backend,
                        const IKnotPolicy& knot_policy) const;

            /**
             * @brief Builds the configured computational grid.
             */
            Grid1D make_grid() const {
                return Grid1D::uniform(cfg_.Ngrid, cfg_.x0, cfg_.x1);
            }

            /**
             * @brief Builds the configured nozzle area model.
             */
            NozzleArea make_area() const {
                if (cfg_.area_preset == NozzleAreaPreset::PaperDefault) {
                    return NozzleArea([](double x) { return 1.0 + 2.2 * (x - 1.5) * (x - 1.5); });
                }
                return cfg_.custom_area;
            }

            /**
             * @brief Builds the quasi-1D inviscid flow driver.
             */
            Quasi1DInviscidDriver make_driver() const {
                return Quasi1DInviscidDriver(make_grid(), make_area(), cfg_.gas);
            }

            /**
             * @brief Builds the boundary-condition applier.
             */
            BoundaryApplier make_bc() const {
                return BoundaryApplier(make_grid(), make_area(), cfg_.gas, cfg_.bc);
            }

            /**
             * @brief Builds a simple uniform initial condition.
             */
            std::vector<UVec> initial_U_simple() const;

            /**
             * @brief Builds the default perturbed steady-state initial condition.
             */
            std::vector<UVec> initial_U() const;

            /**
             * @brief Returns the runner configuration.
             */
            const NozzleConfig& config() const { return cfg_; }

        private:
            NozzleConfig cfg_;
    };

} // namespace nozzle_qns
