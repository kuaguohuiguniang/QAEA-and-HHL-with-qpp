#pragma once
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <functional>

namespace nozzle_qns {

using idx = std::size_t;

    inline void ensure(bool cond, const char* msg) {
        if (!cond) throw std::runtime_error(msg);
    }

    struct UVec {
        double U1 = 0.0; // rho*A
        double U2 = 0.0; // rho*A*v
        double U3 = 0.0; // rho*A*( e/(γ-1) + γ/2 v^2 )
    };

    struct FluxVec { double F1{}, F2{}, F3{}; };
    struct SourceVec { double J1{}, J2{}, J3{}; };

    struct Primitive {
        double rho{};
        double v{};
        double T{};
        double p{};
    };

    struct GasModel {
        double gamma = 1.4;
    };

} // namespace nozzle_qns
