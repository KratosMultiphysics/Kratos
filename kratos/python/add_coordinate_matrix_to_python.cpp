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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.3 $
//
//

#if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE

// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/matrix_python_interface.h"
#include "python/matrix_scalar_operator_python.h"
#include "python/matrix_scalar_assignment_operator_python.h"
#include "python/matrix_vector_operator_python.h"

namespace Kratos
{
namespace Python
{
		
  using namespace boost::python;
	
	
  void  AddCoordinateMatrixToPython()
  {
	  MatrixPythonInterface<coordinate_matrix<double> >::CreateInterface("CoordinateMatrix")
		  .def(init<coordinate_matrix<double>::size_type, coordinate_matrix<double>::size_type>())
		  .def(MatrixScalarOperatorPython<coordinate_matrix<double>, double>())
		  .def(MatrixScalarAssignmentOperatorPython<coordinate_matrix<double>, double>())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, zero_matrix<double>, coordinate_matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, identity_matrix<double>, coordinate_matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, scalar_matrix<double>, coordinate_matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, banded_matrix<double>, banded_matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, triangular_matrix<double, upper>, matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, triangular_matrix<double, lower>, matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, symmetric_matrix<double, upper>, matrix<double> >())
#if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, hermitian_matrix<double, upper>, matrix<double> >())
#endif
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, mapped_matrix<double>, mapped_matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, compressed_matrix<double>, compressed_matrix<double> >())
	          .def(MatrixMatrixOperatorPython<coordinate_matrix<double>, matrix<double>, matrix<double> >())
		  ;
  }
	
}  // namespace Python.

} // Namespace Kratos

#endif  // if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
