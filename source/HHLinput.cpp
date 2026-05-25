#include "HHLinput.hpp"
#include <fstream>

/**
 * @brief Computes the maximum row sparsity of the current matrix A.
 */
void LinearSystem::set_sparsity() {
    d = 0;
    for (int i = 0; i < N; ++i) {
        int row_nz = 0;
        for (int j = 0; j < N; ++j) {
            if (std::abs(A(i, j)) > tol) {
                row_nz++;
            }
        }
        if (row_nz > d) {
            d = row_nz;
        }
    }
}

/**
 * @brief Loads, validates, and returns a linear system from a text file.
 *
 * The expected file format is the dimension N, followed by N x N complex
 * matrix entries and then N complex vector entries.
 */
LinearSystem InputHandler::load_from_file(std::string filename) {
    LinearSystem sys;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    else {
        std::cout << "Loading linear system from file: " << filename << std::endl;
        infile >> sys.N;
        sys.A = qpp::cmat::Zero(sys.N, sys.N);
        sys.b = qpp::ket::Zero(sys.N);
        double re, im;
        char skip = 0;
        for (int i = 0; i < sys.N; ++i) {
            for (int j = 0; j < sys.N; ++j) {
                infile >> re >> im >> skip; ///< Skip the matrix-entry separator.
                sys.A(i, j) = std::complex<double>(re, im);
            }
        }
        for (int i = 0; i < sys.N; ++i) {
            infile >> re >> im;
            sys.b(i) = std::complex<double>(re, im);
        }
    }
    infile.close();
    validate_system(sys);
    sys.set_sparsity();
    return sys;
}

/**
 * @brief Validates Hermiticity, scales A and b, and normalizes b.
 */
void InputHandler::validate_system(LinearSystem& sys) {
    /// Check Hermiticity with ||A - A^dagger|| and embed non-Hermitian A.
    double diff = qpp::norm(sys.A - qpp::adjoint(sys.A));
    if (diff > LinearSystem::tol) {
        std::cout << "A is not Hermitian. Using C=(0, A; A^dagger, 0) instead." << std::endl;
        qpp::cmat C = qpp::cmat::Zero(2 * sys.N, 2 * sys.N);
        C.topRightCorner(sys.N, sys.N) = sys.A;
        C.bottomLeftCorner(sys.N, sys.N) = qpp::adjoint(sys.A);
        sys.A = C;
        qpp::ket b_new = qpp::ket::Zero(2 * sys.N);
        b_new.head(sys.N) = sys.b;
        sys.b = b_new;
        sys.N = 2 * sys.N;
    }
    /// Scale A with its Frobenius norm so ||A||_2 <= ||A||_F <= 1 after scaling.
    double alpha = sys.A.norm(); ///< Frobenius norm computed by Eigen.
    if (alpha < LinearSystem::tol) {
        throw std::runtime_error("Error: Input matrix A is zero.");
    }

    /// Scale both A and b to keep Ax=b equivalent.
    sys.A /= alpha;
    sys.b /= alpha;
    sys.A_scale = alpha;

    std::cout << "Scaled (A,b) by 1/" << alpha
              << " so that ||A||_F <= 1 (and hence ||A||_2 <= 1)." << std::endl;
    sys.norm_b = qpp::norm(sys.b);
    if (sys.norm_b < LinearSystem::tol) {
         throw std::runtime_error("Error: Input vector b is zero.");
    }
    if (std::abs(sys.norm_b - 1.0) > LinearSystem::tol) {
        sys.b /= sys.norm_b;
        std::cout << "Normalized vector b to have unit norm." << std::endl;
    }
}
