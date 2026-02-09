#pragma once
#include "boundary_conditions.hpp"
#include "knot_policy.hpp"
#include "mean_integral.hpp"

namespace nozzle_qns {

    // Time partition parameters
    struct TimePartition {
        idx n = 0;         // number of outer intervals
        idx k = 0;         // hierarchy depth (Nk = n^(k-1))
        double Tfinal = 0.0;

        idx Nk() const;      // n^(k-1)
        double h() const;    // T/n
        double hbar() const; // T/n^k (sub-subinterval length)
    };

    struct TaylorConfig {
        idx r = 1; // Taylor order; higher order to be implemented
    };

    struct StepperConfig {
        TimePartition tp{};
        TaylorConfig taylor{};
        idx knots_per_subsub = 0; // how many knots to sample in u∈[0,1]
        double eps1 = 1e-3;       // target error for mean estimation
        double delta = 1e-3;      // failure probability target
    };

    // Advances solution over one outer interval [t_i, t_{i+1}]
    class OneIntervalStepper {
        public:
            OneIntervalStepper(StepperConfig cfg,
                            Quasi1DInviscidDriver driver,
                            BoundaryApplier bc,
                            MeanIntegralComputer mean_integrals,
                            const IKnotPolicy& knot_policy);

            // Mutates U in place: U <- U at end of one outer interval
            void advance_one_interval(std::vector<UVec>& U) const;

            const StepperConfig& config() const { return cfg_; }

        private:
            StepperConfig cfg_;
            Quasi1DInviscidDriver driver_;
            BoundaryApplier bc_;
            MeanIntegralComputer mean_integrals_;
            const IKnotPolicy& knots_;
    };

} // namespace nozzle_qns
