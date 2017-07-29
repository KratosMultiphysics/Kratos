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


void  AddTriangularMatrixToPython()
{

// 	  MatrixPythonInterface<triangular_matrix<double, upper>, upper>::CreateInterface("UpperTriangularMatrix")
// 		  .def(init<triangular_matrix<double, upper>::size_type, triangular_matrix<double, upper>::size_type>())
// 		  .def(MatrixScalarOperatorPython<triangular_matrix<double, upper>, double, matrix<double> >())
// 		  ;

// 	  MatrixPythonInterface<triangular_matrix<double, lower>, lower>::CreateInterface("LowerTriangularMatrix")
// 		  .def(init<triangular_matrix<double, lower>::size_type, triangular_matrix<double, lower>::size_type>())
// 		  .def(MatrixScalarOperatorPython<triangular_matrix<double, lower>, double, matrix<double> >())
// 		  ;

}

}  // namespace Python.

} // Namespace Kratos

