#pragma once

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

    // Choose built-in nozzle area or user-supplied.
    enum class NozzleAreaPreset {
        PaperDefault,   // A(x) = 1 + 2.2 (x-1.5)^2 on [0,3]
        Custom          // use user-provided NozzleArea
    };

    struct NozzleConfig {
        // ---- Domain / grid ----
        double x0 = 0.0;
        double x1 = 3.0;
        idx Ngrid = 101; // number of grid points

        // ---- Geometry ----
        NozzleAreaPreset area_preset = NozzleAreaPreset::PaperDefault;
        // used if area_preset==Custom
        NozzleArea custom_area = NozzleArea([](double) { return 1.0; });

        // ---- Gas model ----
        GasModel gas{1.4};

        // ---- Boundary conditions ----
        // Inlet is subsonic inflow with rho=1,T=1 (hard-coded in BoundaryApplier)
        BoundaryConfig bc{
            ExitType::SupersonicOutflow,
            0.6784 // pe used only if SubsonicOutflow
        };

        // ---- Time/algorithm parameters ----
        // This is the “ODE + QAEA” layer control knobs.
        StepperConfig stepper_cfg{};
    };

    class NozzleRunner {
        public:
            explicit NozzleRunner(NozzleConfig cfg) : cfg_(std::move(cfg)) {
                ensure(cfg_.Ngrid >= 3, "NozzleRunner: Ngrid must be >= 3");
                ensure(cfg_.x1 > cfg_.x0, "NozzleRunner: x1 must be > x0");
            }

            // Build + run the full pipeline:
            // config -> grid/area/driver/bc -> MeanIntegralComputer -> stepper -> solver -> solution
            Solution run(IMeanEstimator& estimator_backend,
                        const IKnotPolicy& knot_policy) const;

            Grid1D make_grid() const {
                return Grid1D::uniform(cfg_.Ngrid, cfg_.x0, cfg_.x1);
            }

            NozzleArea make_area() const {
                if (cfg_.area_preset == NozzleAreaPreset::PaperDefault) {
                    // paper nozzle: A(x)=1+2.2(x-1.5)^2
                    return NozzleArea([](double x) { return 1.0 + 2.2 * (x - 1.5) * (x - 1.5); });
                }
                return cfg_.custom_area;
            }

            Quasi1DInviscidDriver make_driver() const {
                return Quasi1DInviscidDriver(make_grid(), make_area(), cfg_.gas);
            }

            BoundaryApplier make_bc() const {
                return BoundaryApplier(make_grid(), make_area(), cfg_.gas, cfg_.bc);
            }

            std::vector<UVec> initial_U_simple() const; // A naive initial condition

            std::vector<UVec> initial_U() const; // Steady-state initial condition with randonmized perturbations, non-shock version only

            const NozzleConfig& config() const { return cfg_; }

        private:
            NozzleConfig cfg_;
    };

} // namespace nozzle_qns
