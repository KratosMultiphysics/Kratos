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


#if !defined(KRATOS_ADD_MATRIX_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_MATRIX_TO_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes


namespace Kratos
{

namespace Python
{

void  AddMatrixToPython(pybind11::module& m);
// void  AddBandedMatrixToPython();
// void  AddSymmetricMatrixToPython();
// #if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE
// void  AddHermitianMatrixToPython();
// #endif
// void  AddZeroMatrixToPython();
// void  AddIdentityMatrixToPython();
// void  AddScalarMatrixToPython();
// void  AddTriangularMatrixToPython();
// void  AddSparseMatrixToPython();
// // // // // // // // // // // void  AddCompressedMatrixToPython(pybind11::module& m);
// #if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
// void  AddCoordinateMatrixToPython();
// #endif

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_MATRIX_TO_PYTHON_H_INCLUDED  defined 
