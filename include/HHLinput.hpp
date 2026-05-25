#pragma once
/**
 * @file HHLinput.hpp
 * @brief Declares input parsing and validation utilities for the HHL demo.
 */

#include "qpp/qpp.hpp"

/**
 * @brief Container for a linear system prepared for the HHL solver.
 *
 * The input handler normalizes and, when necessary, Hermitian-embeds the matrix
 * before this structure is passed to HHLSolver.
 */
struct LinearSystem {
    qpp::cmat A; ///< System matrix after validation and scaling.
    qpp::ket b;  ///< Right-hand side vector encoded as a normalized quantum state.
    int N;       ///< Dimension of the validated system matrix.
    int d;       ///< Maximum number of nonzero entries in any row of A.

    /**
     * @brief Computes the row sparsity parameter d for the current matrix A.
     */
    void set_sparsity();

    static constexpr double tol = 1e-12; ///< Numerical tolerance for zero checks.
    double norm_b;                       ///< Norm of b before final state normalization.
    double A_scale;                      ///< Frobenius scaling factor applied to A and b.
};

/**
 * @brief Reads and validates linear-system input for the HHL solver.
 */
class InputHandler {
public:
    /**
     * @brief Loads a linear system from a text file.
     *
     * @param filename Path to the input file.
     * @return Validated LinearSystem ready for HHLSolver.
     */
    static LinearSystem load_from_file(std::string filename);

    /**
     * @brief Validates and normalizes a linear system in place.
     *
     * If A is not Hermitian, the method embeds it into the block matrix
     * C = [[0, A], [A^dagger, 0]]. It then scales A and b by the Frobenius
     * norm of A and normalizes b so it can be used as a quantum state.
     *
     * @param sys Linear system to validate and normalize.
     */
    static void validate_system(LinearSystem& sys);
};
