//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Carlos A. Roig
//
//

#if !defined(KRATOS_ADD_TESTING_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_TESTING_TO_PYTHON_H_INCLUDED


// Project includes
#include <pybind11/pybind11.h>

namespace Kratos
{

namespace Python
{

void  AddTestingToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_TESTING_TO_PYTHON_H_INCLUDED  defined
