#pragma once
/**
 * @file knot_policy.hpp
 * @brief Declares policies for choosing quadrature or sampling knots.
 */

#include "nozzle_types.hpp"
#include <vector>

namespace nozzle_qns {

    /**
     * @brief Collection of normalized knot locations.
     */
    struct KnotGrid {
        std::vector<double> u; ///< Knot positions u in [0, 1].
    };

    /**
     * @brief Interface for constructing knot grids.
     */
    class IKnotPolicy {
        public:
            virtual ~IKnotPolicy() = default;

            /**
             * @brief Creates a knot grid of the requested size.
             *
             * @param count Number of knots to generate.
             * @return Knot grid in the normalized interval [0, 1].
             */
            virtual KnotGrid make_knots(idx count) const = 0;
    };

    /**
     * @brief Uniform knot policy including both endpoints.
     */
    class UniformKnotPolicy final : public IKnotPolicy {
        public:
            /**
             * @brief Creates uniformly spaced knots in [0, 1].
             */
            KnotGrid make_knots(idx count) const override;
    };

} // namespace nozzle_qns
