//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_TRILINOS_ADD_SPACE_TO_PYTHON_H_INCLUDED )
#define  KRATOS_TRILINOS_ADD_SPACE_TO_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"



namespace Kratos
{
namespace Python
{
void  AddBasicOperations(pybind11::module& m);
}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_ADD_SPACE_TO_PYTHON_H_INCLUDED  defined 
