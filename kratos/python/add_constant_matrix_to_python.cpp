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
#include "python/matrix_vector_operator_python.h"

namespace Kratos
{
namespace Python
{

using namespace boost::python;


void  AddConstantMatrixToPython()
{
// 	  ReadonlyMatrixPythonInterface<identity_matrix<double> >::CreateInterface("IdentityMatrix")
// 		  .def(init<identity_matrix<double>::size_type>())
// 		  .def(init<identity_matrix<double>::size_type, identity_matrix<double>::size_type>())
// 		  .def(MatrixScalarOperatorPython<identity_matrix<double>, double, matrix<double>, matrix<double> >())
// 		  .def(MatrixVectorOperatorPython<identity_matrix<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<identity_matrix<double>, zero_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<identity_matrix<double>, unit_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<identity_matrix<double>, scalar_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<identity_matrix<double>, mapped_vector<double>  >())
// 		  .def(MatrixVectorOperatorPython<identity_matrix<double>, compressed_vector<double> >())
// 		  .def(MatrixVectorOperatorPython<identity_matrix<double>, coordinate_vector<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, zero_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, identity_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, scalar_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, banded_matrix<double>, banded_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, triangular_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, triangular_matrix<double, lower>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, symmetric_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, hermitian_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, mapped_matrix<double>, mapped_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, compressed_matrix<double>, compressed_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<identity_matrix<double>, coordinate_matrix<double>, coordinate_matrix<double> >())
// 		  ;
//
// 	  ReadonlyMatrixPythonInterface<zero_matrix<double> >::CreateInterface("ZeroMatrix")
// 		  .def(init<zero_matrix<double>::size_type>())
// 		  .def(init<zero_matrix<double>::size_type, zero_matrix<double>::size_type>())
// 		  .def(MatrixScalarOperatorPython<zero_matrix<double>, double, matrix<double>, matrix<double> >())
// 		  .def(MatrixVectorOperatorPython<zero_matrix<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<zero_matrix<double>, zero_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<zero_matrix<double>, unit_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<zero_matrix<double>, scalar_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<zero_matrix<double>, mapped_vector<double>  >())
// 		  .def(MatrixVectorOperatorPython<zero_matrix<double>, compressed_vector<double> >())
// 		  .def(MatrixVectorOperatorPython<zero_matrix<double>, coordinate_vector<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, zero_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, identity_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, scalar_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, banded_matrix<double>, banded_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, triangular_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, triangular_matrix<double, lower>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, symmetric_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, hermitian_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, mapped_matrix<double>, mapped_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, compressed_matrix<double>, compressed_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<zero_matrix<double>, coordinate_matrix<double>, coordinate_matrix<double> >())
// 		  ;
//
// 	  ReadonlyMatrixPythonInterface<scalar_matrix<double> >::CreateInterface("ScalarMatrix")
// 		  .def(init<scalar_matrix<double>::size_type, scalar_matrix<double>::size_type>())
// 		  .def(init<scalar_matrix<double>::size_type, scalar_matrix<double>::size_type, scalar_matrix<double>::value_type>())
// 		  .def(MatrixScalarOperatorPython<scalar_matrix<double>, double, matrix<double>, matrix<double> >())
// 		  .def(MatrixVectorOperatorPython<scalar_matrix<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<scalar_matrix<double>, zero_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<scalar_matrix<double>, unit_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<scalar_matrix<double>, scalar_vector<double>, vector<double> >())
// 		  .def(MatrixVectorOperatorPython<scalar_matrix<double>, mapped_vector<double>  >())
// 		  .def(MatrixVectorOperatorPython<scalar_matrix<double>, compressed_vector<double> >())
// 		  .def(MatrixVectorOperatorPython<scalar_matrix<double>, coordinate_vector<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, zero_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, identity_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, scalar_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, banded_matrix<double>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, triangular_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, triangular_matrix<double, lower>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, symmetric_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, hermitian_matrix<double, upper>, matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, mapped_matrix<double>, mapped_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, compressed_matrix<double>, compressed_matrix<double> >())
// 	          .def(MatrixMatrixOperatorPython<scalar_matrix<double>, coordinate_matrix<double>, coordinate_matrix<double> >())
// 		  ;

}

}  // namespace Python.

} // Namespace Kratos

