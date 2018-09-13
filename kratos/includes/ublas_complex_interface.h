//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre
//

#if !defined(KRATOS_UBLAS_COMPLEX_INTERFACE_H_INCLUDED )
#define  KRATOS_UBLAS_COMPLEX_INTERFACE_H_INCLUDED

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

namespace Kratos
{

///@name Type Definitions
///@{

using namespace boost::numeric::ublas;

typedef DenseVector<std::complex<double>> ComplexVector;
typedef unit_vector<std::complex<double>> ComplexUnitVector;
typedef zero_vector<std::complex<double>> ComplexZeroVector;
typedef scalar_vector<std::complex<double>> ComplexScalarVector;
typedef mapped_vector<std::complex<double>> ComplexSparseVector;
typedef compressed_vector<std::complex<double>> ComplexCompressedVector;
typedef coordinate_vector<std::complex<double>> ComplexCoordinateVector;
typedef vector_range<ComplexVector> ComplexVectorRange;
typedef vector_slice<ComplexVector> ComplexVectorSlice;
typedef DenseMatrix<std::complex<double>> ComplexMatrix;
typedef identity_matrix<std::complex<double>> ComplexIdentityMatrix;
typedef zero_matrix<std::complex<double>> ComplexZeroMatrix;
typedef scalar_matrix<std::complex<double>> ComplexScalarMatrix;
typedef triangular_matrix<std::complex<double>> ComplexTriangularMatrix;
typedef symmetric_matrix<std::complex<double>> ComplexSymmetricMatrix;
typedef hermitian_matrix<std::complex<double>> ComplexHermitianMatrix;
typedef banded_matrix<std::complex<double>> ComplexBandedMatrix;
typedef mapped_matrix<std::complex<double>> ComplexSparseMatrix;
typedef compressed_matrix<std::complex<double>> ComplexCompressedMatrix;
typedef coordinate_matrix<std::complex<double>> ComplexCoordinateMatrix;
typedef matrix_row<ComplexMatrix> ComplexMatrixRow;
typedef matrix_column<ComplexMatrix> ComplexMatrixColumn;
typedef matrix_vector_range<ComplexMatrix> ComplexMatrixVectorRange;
typedef matrix_vector_slice<ComplexMatrix> ComplexMatrixVectorSlice;
typedef matrix_range<ComplexMatrix> ComplexMatrixRange;
typedef matrix_slice<ComplexMatrix> ComplexMatrixSlice;

///@}

}  // namespace Kratos.

#endif // KRATOS_UBLAS_COMPLEX_INTERFACE_H_INCLUDED  defined 
