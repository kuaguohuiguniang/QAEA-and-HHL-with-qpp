#include "HHLinput.hpp"
#include <fstream>

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
                infile >> re >> im >> skip; // Skip the comma and ;
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

void InputHandler::validate_system(LinearSystem& sys) {
    // Calculate the norm of (A - A^dagger) to check if A is Hermitian
    // If not, construct C = (0 A; A^dagger 0) and update b accordingly
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
    // Normalize/scale A so that ||A|| <= 1
    // Simple safe bound: Frobenius norm
    double alpha = sys.A.norm(); // Frobenius norm in Eigen
    if (alpha < LinearSystem::tol) {
        throw std::runtime_error("Error: Input matrix A is zero.");
    }

    // Scale both A and b to keep Ax=b equivalent
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