//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_CUSTOM_RESPONSE_FUNCTIONS_PYTHON_H_INCLUDED)
#define KRATOS_CUSTOM_RESPONSE_FUNCTIONS_PYTHON_H_INCLUDED

#include "pybind11/pybind11.h"

namespace Kratos
{

namespace Python
{

  void AddCustomResponseFunctionsToPython(pybind11::module& m);

} // namespace Python

} // namespace Kratos

#endif // KRATOS_CUSTOM_RESPONSE_FUNCTIONS_PYTHON_H_INCLUDED
