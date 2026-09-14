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
template <typename TDataType> using DenseMatrix = EigenMatrix<TDataType>;
template <typename TDataType> using DenseVector = EigenVector<TDataType>;
template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using BoundedMatrix = EigenBoundedMatrix<TDataType, TSize1, TSize2>;
template <typename TDataType, std::size_t TSize> using BoundedVector = EigenBoundedVector<TDataType, TSize>;

typedef DenseVector<double> Vector;
typedef DenseMatrix<double> Matrix;

// Lazy factories (no storage, as the uBLAS ones)
typedef unit_vector<double> UnitVector;
typedef zero_vector<double> ZeroVector;
typedef scalar_vector<double> ScalarVector;
typedef identity_matrix<double> IdentityMatrix;
typedef zero_matrix<double> ZeroMatrix;
typedef scalar_matrix<double> ScalarMatrix;

// Proxies (views over the dense containers)
typedef vector_range<Vector> VectorRange;
typedef vector_slice<Vector> VectorSlice;
typedef matrix_column<Matrix> MatrixColumn;
typedef matrix_range<Matrix> MatrixRange;
template <typename TExpressionType> using MatrixRow = matrix_row<TExpressionType>;

// Sparse container (CSR)
typedef EigenCompressedMatrix<double> CompressedMatrix;

// Lowercase boost::numeric::ublas names that Kratos code uses unqualified
// (under the uBLAS backend they resolve through the using-directive of
// ublas_interface.h); mapped onto the Eigen-backed types.
template <typename TDataType> using matrix = EigenMatrix<TDataType>;
template <typename TDataType> using vector = EigenVector<TDataType>;
template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using bounded_matrix = EigenBoundedMatrix<TDataType, TSize1, TSize2>;
template <typename TDataType, std::size_t TSize> using bounded_vector = EigenBoundedVector<TDataType, TSize>;
template <typename TDataType> using compressed_matrix = EigenCompressedMatrix<TDataType>;

// The uBLAS expression-parameter idiom (f(const vector_expression<E>& e) and
// e() to reach the concrete expression) maps onto Eigen's dense expression
// base; the zero-argument call operator is provided by the MatrixBase plugin.
template <typename TExpressionType> using vector_expression = Eigen::MatrixBase<TExpressionType>;
template <typename TExpressionType> using matrix_expression = Eigen::MatrixBase<TExpressionType>;

///@}

}  // namespace Kratos.
