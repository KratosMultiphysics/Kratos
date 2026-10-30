//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre
//

#pragma once

#ifdef KRATOS_USE_EIGEN_BACKEND

// System includes
#include <complex>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/default_interface.h"

namespace Kratos
{

///@name Type Definitions
///@{
// The complex counterparts of the Eigen-backed real types (the same templates
// instantiated on std::complex<double>).

using ComplexVector = DenseVector<std::complex<double>>;
using ComplexUnitVector = unit_vector<std::complex<double>>;
using ComplexZeroVector = zero_vector<std::complex<double>>;
using ComplexScalarVector = scalar_vector<std::complex<double>>;
using ComplexVectorRange = vector_range<ComplexVector>;
using ComplexVectorSlice = vector_slice<ComplexVector>;
using ComplexMatrix = DenseMatrix<std::complex<double>>;
using ComplexIdentityMatrix = identity_matrix<std::complex<double>>;
using ComplexZeroMatrix = zero_matrix<std::complex<double>>;
using ComplexScalarMatrix = scalar_matrix<std::complex<double>>;
using ComplexCompressedMatrix = EigenCompressedMatrix<std::complex<double>>;
using ComplexMatrixRow = matrix_row<ComplexMatrix>;
using ComplexMatrixColumn = matrix_column<ComplexMatrix>;
using ComplexMatrixRange = matrix_range<ComplexMatrix>;

///@}

}  // namespace Kratos.

#else // uBLAS backend

// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

///@name Type Definitions
///@{

using namespace boost::numeric::ublas;

using ComplexVector = DenseVector<std::complex<double>>;
using ComplexUnitVector = unit_vector<std::complex<double>>;
using ComplexZeroVector = zero_vector<std::complex<double>>;
using ComplexScalarVector = scalar_vector<std::complex<double>>;
using ComplexSparseVector = mapped_vector<std::complex<double>>;
using ComplexCompressedVector = compressed_vector<std::complex<double>>;
using ComplexCoordinateVector = coordinate_vector<std::complex<double>>;
using ComplexVectorRange = vector_range<ComplexVector>;
using ComplexVectorSlice = vector_slice<ComplexVector>;
using ComplexMatrix = DenseMatrix<std::complex<double>>;
using ComplexIdentityMatrix = identity_matrix<std::complex<double>>;
using ComplexZeroMatrix = zero_matrix<std::complex<double>>;
using ComplexScalarMatrix = scalar_matrix<std::complex<double>>;
using ComplexTriangularMatrix = triangular_matrix<std::complex<double>>;
using ComplexSymmetricMatrix = symmetric_matrix<std::complex<double>>;
using ComplexHermitianMatrix = hermitian_matrix<std::complex<double>>;
using ComplexBandedMatrix = banded_matrix<std::complex<double>>;
using ComplexSparseMatrix = mapped_matrix<std::complex<double>>;
using ComplexCompressedMatrix = compressed_matrix<std::complex<double>>;
using ComplexCoordinateMatrix = coordinate_matrix<std::complex<double>>;
using ComplexMatrixRow = matrix_row<ComplexMatrix>;
using ComplexMatrixColumn = matrix_column<ComplexMatrix>;
using ComplexMatrixVectorRange = matrix_vector_range<ComplexMatrix>;
using ComplexMatrixVectorSlice = matrix_vector_slice<ComplexMatrix>;
using ComplexMatrixRange = matrix_range<ComplexMatrix>;
using ComplexMatrixSlice = matrix_slice<ComplexMatrix>;

///@}

}  // namespace Kratos.

#endif // KRATOS_USE_EIGEN_BACKEND
