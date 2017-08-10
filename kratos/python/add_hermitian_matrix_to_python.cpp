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
//                   Riccardo Rossi
//


#if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/matrix_python_interface.h"
#include "python/matrix_scalar_operator_python.h"
#include "python/matrix_scalar_assignment_operator_python.h"
#include "python/matrix_matrix_operator_python.h"

namespace Kratos
{
namespace Python
{

using namespace boost::python;


void  AddHermitianMatrixToPython()
{
// 	  MatrixPythonInterface<hermitian_matrix<double, upper>, upper>::CreateInterface("HermitianMatrix")
// 		  .def(init<hermitian_matrix<double, upper>::size_type>())
// 		  .def(init<hermitian_matrix<double, upper>::size_type, hermitian_matrix<double, upper>::size_type>())
// 		  .def(MatrixScalarOperatorPython<hermitian_matrix<double, upper>, double>())
// 		  .def(MatrixScalarAssignmentOperatorPython<hermitian_matrix<double, upper>, double>())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, zero_matrix<double>, hermitian_matrix<double, upper> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, identity_matrix<double>, hermitian_matrix<double, upper> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, scalar_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, banded_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, triangular_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, triangular_matrix<double, lower>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, symmetric_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, mapped_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, compressed_matrix<double>, matrix<double> >())
// #if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
// 	          .def(MatrixMatrixOperatorPython<hermitian_matrix<double, upper>, coordinate_matrix<double>, matrix<double> >())
// #endif
// 		  ;

}

}  // namespace Python.

} // Namespace Kratos

#endif // if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE
