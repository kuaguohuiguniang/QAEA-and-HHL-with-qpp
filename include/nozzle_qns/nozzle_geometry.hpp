#pragma once
#include "nozzle_types.hpp"
#include <cmath>

namespace nozzle_qns {

    class NozzleArea {
    public:
        using Func = std::function<double(double)>;

        explicit NozzleArea(Func A_of_x) : A_(std::move(A_of_x)) {
            ensure(static_cast<bool>(A_), "NozzleArea: empty function");
        }

        double operator()(double x) const {
            double A = A_(x);
            ensure(A > 0.0, "NozzleArea: A(x) must be > 0");
            return A;
        }

    private:
        Func A_;
    };

    struct Grid1D {
        idx m{};
        double x0{};
        double x1{};
        double dx{};

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

        double x(idx j) const {
            ensure(j < m, "Grid1D: index out of range");
            return x0 + double(j) * dx;
        }
    };

} // namespace nozzle_qns
