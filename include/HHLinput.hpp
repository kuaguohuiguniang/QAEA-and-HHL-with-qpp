#pragma once
#include "qpp/qpp.hpp"

struct LinearSystem {
    qpp::cmat A;
    qpp::ket b;
    // Dimension N
    int N;
    // Sparsity d
    int d;
    void set_sparsity();
    // Tolerance for considering an element as non-zero
    static constexpr double tol = 1e-12;
    double norm_b;
    double A_scale;
};

class InputHandler {
public:
    // Reads dimension, matrix A, and vector b from file
    static LinearSystem load_from_file(std::string filename);

    // Verifies A is Hermitian (using C instead if not)
    // normalized A and b to have ||A||<= 1
    // Normalizes b to be a valid quantum state
    static void validate_system(LinearSystem& sys);
};
