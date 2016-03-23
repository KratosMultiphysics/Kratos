//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//











#if !defined(KRATOS_UBLAS_INTERFACE_H_INCLUDED )
#define  KRATOS_UBLAS_INTERFACE_H_INCLUDED



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

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

using namespace boost::numeric::ublas;

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
typedef compressed_matrix<double> CompressedMatrix;
typedef coordinate_matrix<double> CoordinateMatrix;
typedef matrix_row<Matrix> MatrixRow;
typedef matrix_column<Matrix> MatrixColumn;
typedef matrix_vector_range<Matrix> MatrixVectorRange;
typedef matrix_vector_slice<Matrix> MatrixVectorSlice;
typedef matrix_range<Matrix> MatrixRange;
typedef matrix_slice<Matrix> MatrixSlice;

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

}  // namespace Kratos.

#endif // KRATOS_UBLAS_INTERFACE_H_INCLUDED  defined 



