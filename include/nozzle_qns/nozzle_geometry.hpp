#pragma once
/**
 * @file nozzle_geometry.hpp
 * @brief Defines the one-dimensional nozzle grid and area model.
 */

#include "nozzle_types.hpp"
#include <cmath>

namespace nozzle_qns {

    /**
     * @brief Callable nozzle cross-sectional area model A(x).
     */
    class NozzleArea {
    public:
        using Func = std::function<double(double)>;

        /**
         * @brief Constructs an area model from a user-provided function.
         *
         * @param A_of_x Function returning the area at position x.
         */
        explicit NozzleArea(Func A_of_x) : A_(std::move(A_of_x)) {
            ensure(static_cast<bool>(A_), "NozzleArea: empty function");
        }

        /**
         * @brief Evaluates the nozzle area at a grid coordinate.
         *
         * @param x Axial coordinate.
         * @return Positive area A(x).
         */
        double operator()(double x) const {
            double A = A_(x);
            ensure(A > 0.0, "NozzleArea: A(x) must be > 0");
            return A;
        }

    private:
        Func A_;
    };

    /**
     * @brief Uniform one-dimensional grid over the nozzle domain.
     */
    struct Grid1D {
        idx m{};     ///< Number of grid points.
        double x0{}; ///< Left boundary coordinate.
        double x1{}; ///< Right boundary coordinate.
        double dx{}; ///< Uniform grid spacing.

        /**
         * @brief Creates a uniform grid on [x_begin, x_end].
         *
         * @param m_points Number of grid points.
         * @param x_begin Left boundary coordinate.
         * @param x_end Right boundary coordinate.
         * @return Constructed uniform grid.
         */
        static Grid1D uniform(idx m_points, double x_begin, double x_end) {
            ensure(m_points >= 3, "Grid1D: need >=3 points");
            ensure(x_end > x_begin, "Grid1D: x_end must be > x_begin");
            Grid1D g;
            g.m = m_points;
            g.x0 = x_begin;
            g.x1 = x_end;
            g.dx = (x_end - x_begin) / double(m_points - 1);
            return g;
        }

        /**
         * @brief Returns the coordinate of grid point j.
         *
         * @param j Grid-point index.
         * @return Coordinate x_j.
         */
        double x(idx j) const {
            ensure(j < m, "Grid1D: index out of range");
            return x0 + double(j) * dx;
        }
    };

} // namespace nozzle_qns
