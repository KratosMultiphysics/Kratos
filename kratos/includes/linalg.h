//   |   /
//   ' /   __| _` | __|  _ \   __|
//   . \  |   (   | |   (   | \__ `
// _|\_\_|  \__,_|\__|\___/ ____/
//         Multi-Physics
//
// License:         BSD License
//                  Kratos default license: kratos/license.txt
//
// Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <iostream>

// External includes
#include <experimental/linalg>

// Project includes

namespace Kratos::Future {

/// @name Linear Algebra Property Tags
/// @{

// Tag Types

/**
 * @brief Alias for linalg transpose operation tag type.
 */
using transpose_t = std::experimental::linalg::transpose_t;

/**
 * @brief Alias for linalg no-transpose operation tag type.
 */
using no_transpose_t = std::experimental::linalg::no_transpose_t;

/**
 * @brief Alias for linalg upper-triangle matrix property tag type.
 */
using upper_triangle_t = std::experimental::linalg::upper_triangle_t;

/**
 * @brief Alias for linalg lower-triangle matrix property tag type.
 */
using lower_triangle_t = std::experimental::linalg::lower_triangle_t;

/**
 * @brief Alias for linalg unit-diagonal matrix property tag type.
 */
using unit_diagonal_t = std::experimental::linalg::unit_diagonal_t;

/**
 * @brief Alias for linalg implicit-unit-diagonal matrix property tag type.
 */
using implicit_unit_diagonal_t = std::experimental::linalg::implicit_unit_diagonal_t;


// Tag Instances

/**
 * @brief Static instance of the transpose operation tag.
 */
inline constexpr transpose_t transpose = std::experimental::linalg::transpose;

/**
 * @brief Static instance of the no-transpose operation tag.
 */
inline constexpr no_transpose_t no_transpose = std::experimental::linalg::no_transpose;

/**
 * @brief Static instance of the upper-triangle matrix property tag.
 */
inline constexpr upper_triangle_t upper_triangle = std::experimental::linalg::upper_triangle;

/**
 * @brief Static instance of the lower-triangle matrix property tag.
 */
inline constexpr lower_triangle_t lower_triangle = std::experimental::linalg::lower_triangle;

/**
 * @brief Static instance of the unit-diagonal matrix property tag.
 */
inline constexpr unit_diagonal_t unit_diagonal = std::experimental::linalg::unit_diagonal;

/**
 * @brief Static instance of the implicit-unit-diagonal matrix property tag.
 */
inline constexpr implicit_unit_diagonal_t implicit_unit_diagonal = std::experimental::linalg::implicit_unit_diagonal;

/// @}
/// @name BLAS Level 1 Operations (Vector-Vector)
/// @{

/**
 * @brief Computes `y = alpha * x + y` (scaled vector addition).
 */
using std::experimental::linalg::axpy;

/**
 * @brief Computes the sum of absolute values (L1 norm).
 */
using std::experimental::linalg::asum;

/**
 * @brief Copies vector x into vector y.
 */
using std::experimental::linalg::copy;

/**
 * @brief Computes the dot product of two vectors.
 */
using std::experimental::linalg::dot;

/**
 * @brief Computes the dot product, conjugating the first vector.
 */
using std::experimental::linalg::dotc;

/**
 * @brief Finds the index of the element with the maximum absolute value.
 */
using std::experimental::linalg::iamax;

/**
 * @brief Computes the Euclidean (L2) norm of a vector.
 */
using std::experimental::linalg::nrm2;

/**
 * @brief Multiplies all elements of a vector by a scalar.
 */
using std::experimental::linalg::scale;

/**
 * @brief Swaps the contents of two vectors.
 */
using std::experimental::linalg::swap_elements;

/// @}
/// @name BLAS Level 2 Operations (Matrix-Vector)
/// @{

/**
 * @brief Computes `A = alpha * x * y^T + A` (general rank-1 update).
 */
using std::experimental::linalg::ger;

/**
 * @brief Computes `A = alpha * x * y^H + A` (general rank-1 update, conjugating y).
 */
using std::experimental::linalg::gerc;

/**
 * @brief Computes `y = A * x` for a Hermitian matrix A.
 */
using std::experimental::linalg::hemv;

/**
 * @brief Computes `y = A * x` (general matrix-vector product).
 */
using std::experimental::linalg::matvec;

/**
 * @brief Computes `y = A * x` for a symmetric matrix A.
 */
using std::experimental::linalg::symv;

/**
 * @brief Computes `x = A * x` for a triangular matrix A.
 */
using std::experimental::linalg::trmv;

/**
 * @brief Solves `A * x_new = x_old` for a triangular matrix A.
 */
using std::experimental::linalg::trsv;

/// @}
/// @name BLAS Level 3 Operations (Matrix-Matrix)
/// @{

/**
 * @brief Computes `C = A * B` (general matrix-matrix product).
 */
using std::experimental::linalg::matmul;

/**
 * @brief Performs a rank-k update for a Hermitian matrix: `C = A * A^H + C`.
 */
using std::experimental::linalg::herk;

/**
 * @brief Performs a rank-2k update for a Hermitian matrix: `C = A * B^H + B * A^H + C`.
 */
using std::experimental::linalg::her2k;

/**
 * @brief Performs a rank-k update for a symmetric matrix: `C = A * A^T + C`.
 */
using std::experimental::linalg::syrk;

/**
 * @brief Performs a rank-2k update for a symmetric matrix: `C = A * B^T + B * A^T + C`.
 */
using std::experimental::linalg::syr2k;

/**
 * @brief Computes `B = A * B` where A is a triangular matrix.
 */
using std::experimental::linalg::trmm;

/**
 * @brief Solves `A * X = B` where A is a triangular matrix and X, B are dense matrices.
 */
using std::experimental::linalg::trsm;

/// @}

} // namespace Kratos::Future