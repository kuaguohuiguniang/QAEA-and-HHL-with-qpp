# include "nozzle_qns/interval_stepper.hpp"
#include <algorithm>
#include <cmath>
#include <utility>

namespace nozzle_qns {

    idx TimePartition::Nk() const {
        if (n == 0 || k == 0) return 0;
        idx r = 1;
        for (idx i = 0; i < k - 1; ++i) r *= n;
        return r;
    }

    double TimePartition::h() const {
        return n > 0 ? Tfinal / n : 0.0;
    }

    double TimePartition::hbar() const {
        idx Nk_ = Nk();
        return (n > 0 && Nk_ > 0) ? (Tfinal / (static_cast<double>(n) * static_cast<double>(Nk_))) : 0.0;
    }


    OneIntervalStepper::OneIntervalStepper(StepperConfig cfg,
                                        Quasi1DInviscidDriver driver,
                                        BoundaryApplier bc,
                                        MeanIntegralComputer mean_integrals,
                                        const IKnotPolicy& knot_policy)
        : cfg_(std::move(cfg)),
          driver_(std::move(driver)),
          bc_(std::move(bc)),
          mean_integrals_(std::move(mean_integrals)),
          knots_(knot_policy) {}
    
    void OneIntervalStepper::advance_one_interval(std::vector<UVec>& U) const {
        // --- Basic checks ---
        ensure(cfg_.tp.n >= 2, "OneIntervalStepper: tp.n must be >= 2");
        ensure(cfg_.tp.k >= 1, "OneIntervalStepper: tp.k must be >= 1");
        ensure(cfg_.tp.Tfinal > 0.0, "OneIntervalStepper: tp.Tfinal must be > 0");
        ensure(cfg_.knots_per_subsub >= 2, "OneIntervalStepper: need >=2 knots");
        ensure(U.size() == driver_.grid().m, "OneIntervalStepper: U size mismatch with grid");

        const idx m  = driver_.grid().m;
        const idx Nk = cfg_.tp.Nk();
        const double hbar = cfg_.tp.hbar();

        // Knot points u in [0,1]
        const KnotGrid kg = knots_.make_knots(cfg_.knots_per_subsub);
        ensure(!kg.u.empty(), "OneIntervalStepper: knot grid is empty");

        const double total_calls = 3.0 * static_cast<double>(cfg_.tp.n) *
                                    static_cast<double>(Nk) * static_cast<double>(m > 2 ? (m - 2) : 1);

        const double delta_per_call =
            cfg_.delta / std::max(1.0, total_calls);

        // Main loop over sub-subintervals
        for (idx l = 0; l < Nk; ++l) {
            // Enforce boundary conditions on current state
            bc_.apply(U);

            // Save starting state for this sub-subinterval
            const std::vector<UVec> U_start = U;

            // 1st-order Taylor base derivative at sub-subinterval start:
            // f0 = f(U_start)
            const std::vector<UVec> f0 = driver_.dUdt_all(U_start);

            // Storage: for each grid point j, collect knot samples of each component of f
            std::vector<std::vector<double>> f1_knots(m), f2_knots(m), f3_knots(m);
            for (idx j = 0; j < m; ++j) {
                f1_knots[j].reserve(kg.u.size());
                f2_knots[j].reserve(kg.u.size());
                f3_knots[j].reserve(kg.u.size());
            }

            // Evaluate f(U_approx(u)) at each knot
            for (double u : kg.u) {
                // --- Build Taylor-approximated state at this knot ---
                // Paper local model: l_i,l(j,u). Here: first-order Taylor in time.
                std::vector<UVec> U_knot = U_start;

                const double tau = u * hbar;

                // r=1 Taylor: U_knot = U_start + tau * f0
                for (idx j = 0; j < m; ++j) {
                    U_knot[j].U1 += tau * f0[j].U1;
                    U_knot[j].U2 += tau * f0[j].U2;
                    U_knot[j].U3 += tau * f0[j].U3;
                }

                bc_.apply(U_knot);

                // Evaluate f at this knot
                const std::vector<UVec> f_knot = driver_.dUdt_all(U_knot);

                // Store knot samples per grid point and component
                for (idx j = 0; j < m; ++j) {
                    f1_knots[j].push_back(f_knot[j].U1);
                    f2_knots[j].push_back(f_knot[j].U2);
                    f3_knots[j].push_back(f_knot[j].U3);
                }
            }

            // Compute mean(f) at each grid point using MeanIntegralComputer (backend inside)
            std::vector<UVec> mean_f(m); // boundaries remain 0; interior computed

            for (idx j = 1; j + 1 < m; ++j) {
                const IntegralEstimate est =
                    mean_integrals_.estimate_mean_f_from_knots(
                        f1_knots[j], f2_knots[j], f3_knots[j],
                        cfg_.eps1, delta_per_call);

                mean_f[j].U1 = est.mean_f1;
                mean_f[j].U2 = est.mean_f2;
                mean_f[j].U3 = est.mean_f3;
            }

            // Update: U_end = U_start + hbar * mean_f
            U = U_start;
            for (idx j = 1; j + 1 < m; ++j) {
                U[j].U1 += hbar * mean_f[j].U1;
                U[j].U2 += hbar * mean_f[j].U2;
                U[j].U3 += hbar * mean_f[j].U3;
            }

            // Enforce boundaries after update (optional but usually a good idea)
            bc_.apply(U);
        }
    }

} // namespace nozzle_qns