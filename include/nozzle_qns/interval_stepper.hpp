#pragma once
/**
 * @file interval_stepper.hpp
 * @brief Declares one-interval time stepping for the nozzle ODE solver.
 */

#include "boundary_conditions.hpp"
#include "knot_policy.hpp"
#include "mean_integral.hpp"

namespace nozzle_qns {

    /**
     * @brief Time-partition parameters for the nested ODE integration scheme.
     */
    struct TimePartition {
        idx n = 0;           ///< Number of outer intervals.
        idx k = 0;           ///< Hierarchy depth; Nk = n^(k - 1).
        double Tfinal = 0.0; ///< Final integration time.

        /**
         * @brief Returns n^(k - 1), the number of sub-subintervals per interval.
         */
        idx Nk() const;

        /**
         * @brief Returns the outer interval length T / n.
         */
        double h() const;

        /**
         * @brief Returns the sub-subinterval length T / n^k.
         */
        double hbar() const;
    };

    /**
     * @brief Taylor expansion configuration for local state approximation.
     */
    struct TaylorConfig {
        idx r = 1; ///< Taylor order; higher orders are reserved for future work.
    };

    /**
     * @brief Full configuration for one-interval stepping.
     */
    struct StepperConfig {
        TimePartition tp{};
        TaylorConfig taylor{};
        idx knots_per_subsub = 0; ///< Number of knots to sample in u in [0, 1].
        double eps1 = 1e-3;       ///< Target error for mean estimation.
        double delta = 1e-3;      ///< Failure-probability target.
    };

    /**
     * @brief Advances a nozzle state over one outer interval.
     */
    class OneIntervalStepper {
        public:
            /**
             * @brief Constructs a stepper from all numerical components.
             */
            OneIntervalStepper(StepperConfig cfg,
                            Quasi1DInviscidDriver driver,
                            BoundaryApplier bc,
                            MeanIntegralComputer mean_integrals,
                            const IKnotPolicy& knot_policy);

            /**
             * @brief Advances U in place over one outer interval.
             *
             * @param U Conservative state vector to update.
             */
            void advance_one_interval(std::vector<UVec>& U) const;

            /**
             * @brief Returns the stepper configuration.
             */
            const StepperConfig& config() const { return cfg_; }

        private:
            StepperConfig cfg_;
            Quasi1DInviscidDriver driver_;
            BoundaryApplier bc_;
            MeanIntegralComputer mean_integrals_;
            const IKnotPolicy& knots_;
    };

} // namespace nozzle_qns
