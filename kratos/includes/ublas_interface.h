//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes
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
#include "includes/define.h"
#ifdef KRATOS_USE_EIGEN_BACKEND
#include "includes/eigen_bounded_types.h"
#endif

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    using namespace boost::numeric::ublas;

    template <typename TDataType> using DenseMatrix = boost::numeric::ublas::matrix<TDataType>;
    template <typename TDataType> using DenseVector = boost::numeric::ublas::vector<TDataType>;

#ifdef KRATOS_USE_EIGEN_BACKEND
    // Under the Eigen backend the fixed-size dense types are Eigen-backed
    // (see eigen_bounded_types.h); the dynamic dense types above remain uBLAS
    // and the mixed idioms are provided by eigen_ublas_compat_operations.h,
    // included at the end of this header.
    template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using BoundedMatrix = EigenBoundedMatrix<TDataType, TSize1, TSize2>;
    template <typename TDataType, std::size_t TSize> using BoundedVector = EigenBoundedVector<TDataType, TSize>;
#else
    template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using BoundedMatrix = boost::numeric::ublas::bounded_matrix<TDataType, TSize1, TSize2>;
    template <typename TDataType, std::size_t TSize> using BoundedVector = boost::numeric::ublas::bounded_vector<TDataType, TSize>;
#endif


    typedef boost::numeric::ublas::vector<double> Vector;
    typedef unit_vector<double> UnitVector;
    typedef zero_vector<double> ZeroVector;
    typedef scalar_vector<double> ScalarVector;
    //typedef sparse_vector<double> SparseVector;
    typedef mapped_vector<double> SparseVector;

    typedef compressed_vector<double> CompressedVector;
    typedef coordinate_vector<double> CoordinateVector;
    typedef vector_range<Vector> VectorRange;
    typedef vector_slice<Vector> VectorSlice;

    typedef matrix<double> Matrix;
    typedef identity_matrix<double> IdentityMatrix;
    typedef zero_matrix<double> ZeroMatrix;
    typedef scalar_matrix<double> ScalarMatrix;
    typedef triangular_matrix<double> TriangularMatrix;
    typedef symmetric_matrix<double> SymmetricMatrix;
    typedef hermitian_matrix<double> HermitianMatrix;
    typedef banded_matrix<double> BandedMatrix;
    //typedef sparse_matrix<double> SparseMatrix;
    typedef mapped_matrix<double> SparseMatrix;
    typedef coordinate_matrix<double> CoordinateMatrix;
    typedef matrix_column<Matrix> MatrixColumn;
    typedef matrix_vector_range<Matrix> MatrixVectorRange;
    typedef matrix_vector_slice<Matrix> MatrixVectorSlice;
    typedef matrix_range<Matrix> MatrixRange;
    typedef matrix_slice<Matrix> MatrixSlice;

	template <typename TExpressionType> using MatrixRow = matrix_row<TExpressionType>;

    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrix;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}
}  // namespace Kratos.

#ifdef KRATOS_USE_EIGEN_BACKEND
// Mixed uBLAS/Eigen idioms (noalias, prod, inner_prod, ... across the two
// worlds). Included after the namespace so the header sees the type aliases
// defined above.
#include "includes/eigen_ublas_compat_operations.h"
#endif