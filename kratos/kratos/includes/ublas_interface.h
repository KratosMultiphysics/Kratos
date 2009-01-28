/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.3 $
//
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


