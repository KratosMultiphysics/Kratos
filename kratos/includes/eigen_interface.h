//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Riccardo Rossi
//

#pragma once

// The Eigen linear-algebra backend (KRATOS_LINEAR_ALGEBRA_BACKEND=eigen, i.e.
// KRATOS_USE_EIGEN_BACKEND): every Kratos dense and sparse linear-algebra
// alias is Eigen-backed, and the uBLAS free-function idiom (prod, noalias,
// trans, row, project, ...) is provided in namespace Kratos with a pure-Eigen
// implementation. This header is the Eigen arm of includes/ublas_interface.h
// and is not meant to be included directly.
//
// The Kratos CMake defines the two Eigen plugin headers globally when the
// Eigen backend is selected (the uBLAS-style size1()/size2()/clear() on every
// dense expression and resize(..., preserve) on every plain object are part
// of the backend contract), so refuse to build without them.
#if !defined(EIGEN_MATRIXBASE_PLUGIN) || !defined(EIGEN_PLAINOBJECTBASE_PLUGIN)
#error "KRATOS_USE_EIGEN_BACKEND requires EIGEN_MATRIXBASE_PLUGIN and EIGEN_PLAINOBJECTBASE_PLUGIN to be defined globally (kratos/includes/eigen_matrixbase_plugin.h and eigen_plainobjectbase_plugin.h); configure with KRATOS_LINEAR_ALGEBRA_BACKEND=eigen through the Kratos CMake."
#endif

// System includes
#include <complex>
#include <cstddef>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>
// The boost uBLAS headers are still made available so that the few
// uBLAS-only helpers spelling boost::numeric::ublas:: explicitly keep
// compiling; nothing in namespace Kratos refers to them under this backend.
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

// Project includes
#include "includes/eigen_dense_types.h"
#include "includes/eigen_sparse_types.h"
#include "includes/eigen_operations.h"

namespace Kratos
{

///@name Type Definitions
///@{

// Dense containers
/// Dense matrix of the Eigen backend. Alias of EigenMatrix<TDataType>, a dynamic, storage-owning matrix that is
/// row-major by default.
template <typename TDataType> using DenseMatrix = EigenMatrix<TDataType>;
/// Dense column vector of the Eigen backend. Alias of EigenVector<TDataType>, a dynamic, storage-owning column vector.
template <typename TDataType> using DenseVector = EigenVector<TDataType>;
/// Fixed-size dense matrix of the Eigen backend. Alias of EigenBoundedMatrix<TDataType, TSize1, TSize2>, a matrix whose
/// dimensions are known at compile time.
template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using BoundedMatrix = EigenBoundedMatrix<TDataType, TSize1, TSize2>;
/// Fixed-size dense column vector of the Eigen backend. Alias of EigenBoundedVector<TDataType, TSize>, a vector whose
/// length is known at compile time.
template <typename TDataType, std::size_t TSize> using BoundedVector = EigenBoundedVector<TDataType, TSize>;

/// Dense vector of doubles. Default instance of DenseVector, the uBLAS-style vector type of the Kratos Eigen backend.
using Vector = DenseVector<double>;
/// Dense matrix of doubles. Default instance of DenseMatrix, the uBLAS-style matrix type of the Kratos Eigen backend.
using Matrix = DenseMatrix<double>;

// Lazy factories (no storage, as the uBLAS ones)
/// Lazy unit vector of doubles. A vector with one entry equal to one and all others zero. It is evaluated on access and
/// has no storage, as boost::numeric::ublas::unit_vector.
using UnitVector = unit_vector<double>;
/// Lazy zero vector of doubles. A vector with all entries zero. It is evaluated on access and has no storage, as
/// boost::numeric::ublas::zero_vector.
using ZeroVector = zero_vector<double>;
/// Lazy vector of doubles with one constant value. A vector with all entries equal to a given value. It is evaluated on
/// access and has no storage, as boost::numeric::ublas::scalar_vector.
using ScalarVector = scalar_vector<double>;
/// Lazy identity matrix of doubles. A square matrix with ones on the diagonal and zeros elsewhere. It is evaluated on
/// access and has no storage, as boost::numeric::ublas::identity_matrix.
using IdentityMatrix = identity_matrix<double>;
/// Lazy zero matrix of doubles. A matrix with all entries zero. It is evaluated on access and has no storage, as
/// boost::numeric::ublas::zero_matrix.
using ZeroMatrix = zero_matrix<double>;
/// Lazy matrix of doubles with one constant value. A matrix with all entries equal to a given value. It is evaluated on
/// access and has no storage, as boost::numeric::ublas::scalar_matrix.
using ScalarMatrix = scalar_matrix<double>;

// Proxies (views over the dense containers)
/// Writable view of a vector subrange. A contiguous range of a Vector, as boost::numeric::ublas::vector_range.
using VectorRange = vector_range<Vector>;
/// Writable strided view of a vector. A slice over the contiguous storage of a Vector, as
/// boost::numeric::ublas::vector_slice.
using VectorSlice = vector_slice<Vector>;
/// Writable view of one matrix column. A column of a Matrix seen as a vector, as boost::numeric::ublas::matrix_column.
using MatrixColumn = matrix_column<Matrix>;
/// Writable view of a matrix block. A rectangular block of a Matrix, as boost::numeric::ublas::matrix_range.
using MatrixRange = matrix_range<Matrix>;
/// Writable view of one matrix row. Alias template of matrix_row, a column-shaped view of one row of a matrix, as
/// boost::numeric::ublas::matrix_row.
template <typename TExpressionType> using MatrixRow = matrix_row<TExpressionType>;

// Sparse container (CSR)
/// Compressed sparse row matrix of doubles. Alias of EigenCompressedMatrix<double>, the sparse matrix type of the Eigen
/// backend in compressed sparse row (CSR) storage.
using CompressedMatrix = EigenCompressedMatrix<double>;

// Storage-order tags, as boost::numeric::ublas::row_major/column_major: a
// handful of call sites request a column-major dense matrix explicitly (e.g.
// a local buffer handed to a Fortran/LAPACK-style routine such as FEAST).
/**
 * @brief Row-major storage-order tag
 * @details Selects row-major storage for a dense matrix type, as boost::numeric::ublas::row_major. The matching Eigen
 * option is looked up through Internals::EigenLayoutOf.
 */
struct row_major {};
/**
 * @brief Column-major storage-order tag
 * @details Selects column-major storage for a dense matrix type, as boost::numeric::ublas::column_major. The matching
 * Eigen option is looked up through Internals::EigenLayoutOf.
 */
struct column_major {};

namespace Internals
{
/**
 * @brief Maps a storage-order tag onto the matching Eigen storage option
 * @details Only the row_major and column_major specializations below are defined. Their static member value holds
 * Eigen::RowMajor or Eigen::ColMajor, respectively.
 * @tparam TLayout The storage-order tag, row_major or column_major.
 */
template <class TLayout> struct EigenLayoutOf;
/**
 * @brief Eigen storage option of the row_major tag
 * @details Specialization of EigenLayoutOf for row_major.
 */
template <> struct EigenLayoutOf<row_major>
{
    /**
     * @brief Eigen row-major storage option
     * @details Equal to Eigen::RowMajor.
     */
    static constexpr int value = Eigen::RowMajor;
};
/**
 * @brief Eigen storage option of the column_major tag
 * @details Specialization of EigenLayoutOf for column_major.
 */
template <> struct EigenLayoutOf<column_major>
{
    /**
     * @brief Eigen column-major storage option
     * @details Equal to Eigen::ColMajor.
     */
    static constexpr int value = Eigen::ColMajor;
};
} // namespace Internals

// Lowercase boost::numeric::ublas names that Kratos code uses unqualified
// (under the uBLAS backend they resolve through the using-directive of
// ublas_interface.h); mapped onto the Eigen-backed types. matrix<T> defaults
// to row_major, as ublas::matrix<T> does; matrix<T, column_major> selects the
// column-major EigenMatrix instantiation.
/// uBLAS-style dense matrix. Alias of EigenMatrix<TDataType> with the storage order given by TLayout. It is row_major
/// by default, as boost::numeric::ublas::matrix.
template <typename TDataType, typename TLayout = row_major> using matrix = EigenMatrix<TDataType, Internals::EigenLayoutOf<TLayout>::value>;
/// uBLAS-style dense vector. Alias of EigenVector<TDataType>.
template <typename TDataType> using vector = EigenVector<TDataType>;
/// uBLAS-style fixed-size dense matrix. Alias of EigenBoundedMatrix<TDataType, TSize1, TSize2>.
template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using bounded_matrix = EigenBoundedMatrix<TDataType, TSize1, TSize2>;
/// uBLAS-style fixed-size dense vector. Alias of EigenBoundedVector<TDataType, TSize>.
template <typename TDataType, std::size_t TSize> using bounded_vector = EigenBoundedVector<TDataType, TSize>;
/// uBLAS-style compressed sparse matrix. Alias of EigenCompressedMatrix<TDataType>, a compressed sparse row matrix.
template <typename TDataType> using compressed_matrix = EigenCompressedMatrix<TDataType>;

// The uBLAS expression-parameter idiom (f(const vector_expression<E>& e) and
// e() to reach the concrete expression) maps onto Eigen's dense expression
// base; the zero-argument call operator is provided by the MatrixBase plugin.
/// Base of vector expressions. Alias of Eigen::MatrixBase<TExpressionType>, so that a function taking a
/// vector_expression accepts any dense Eigen expression.
template <typename TExpressionType> using vector_expression = Eigen::MatrixBase<TExpressionType>;
/// Base of matrix expressions. Alias of Eigen::MatrixBase<TExpressionType>, the same base as vector_expression.
template <typename TExpressionType> using matrix_expression = Eigen::MatrixBase<TExpressionType>;

///@}

}  // namespace Kratos.
