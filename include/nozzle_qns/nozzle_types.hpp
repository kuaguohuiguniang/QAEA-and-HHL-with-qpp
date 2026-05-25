#pragma once
/**
 * @file nozzle_types.hpp
 * @brief Defines core scalar types and state vectors for the nozzle solver.
 */

#include <cstddef>
#include <stdexcept>
#include <vector>
#include <functional>

namespace nozzle_qns {

using idx = std::size_t;

    /**
     * @brief Throws a runtime error when a required condition is false.
     *
     * @param cond Condition to check.
     * @param msg Error message used if the condition fails.
     */
    inline void ensure(bool cond, const char* msg) {
        if (!cond) throw std::runtime_error(msg);
    }

    /**
     * @brief Conservative quasi-1D nozzle state vector.
     */
    struct UVec {
        double U1 = 0.0; ///< rho * A.
        double U2 = 0.0; ///< rho * A * v.
        double U3 = 0.0; ///< rho * A * (T / (gamma - 1) + gamma / 2 * v^2).
    };

    /**
     * @brief Inviscid flux vector for the conservative variables.
     */
    struct FluxVec { double F1{}, F2{}, F3{}; };

    /**
     * @brief Geometric source vector for the quasi-1D equations.
     */
    struct SourceVec { double J1{}, J2{}, J3{}; };

    /**
     * @brief Primitive flow variables derived from a conservative state.
     */
    struct Primitive {
        double rho{}; ///< Density.
        double v{};   ///< Flow velocity.
        double T{};   ///< Temperature.
        double p{};   ///< Pressure.
    };

    /**
     * @brief Gas-law parameters for the dimensionless nozzle model.
     */
    struct GasModel {
        double gamma = 1.4; ///< Ratio of specific heats.
    };

} // namespace nozzle_qns
