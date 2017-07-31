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


void  AddMatrixToPython()
{
    MatrixPythonInterface<matrix<double> >::CreateInterface("Matrix")
    .def(init<matrix<double>::size_type, matrix<double>::size_type>())
    /*		  .def(MatrixScalarOperatorPython<matrix<double>, double>())
      .def(MatrixScalarAssignmentOperatorPython<matrix<double>, double>())
          .def(MatrixMatrixOperatorPython<matrix<double>, zero_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, identity_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, scalar_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, banded_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, triangular_matrix<double, upper>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, triangular_matrix<double, lower>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, symmetric_matrix<double, upper>, matrix<double> >())
    #if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE
          .def(MatrixMatrixOperatorPython<matrix<double>, hermitian_matrix<double, upper>, matrix<double> >())
    #endif
          .def(MatrixMatrixOperatorPython<matrix<double>, mapped_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, compressed_matrix<double>, matrix<double> >())
    #if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
          .def(MatrixMatrixOperatorPython<matrix<double>, coordinate_matrix<double>, matrix<double> >())
    #endif*/
    ;
}

}  // namespace Python.

} // Namespace Kratos

