//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre
//

#if !defined(KRATOS_RESPONSE_FUNCTIONS_PYTHON_H_INCLUDED)
#define KRATOS_RESPONSE_FUNCTIONS_PYTHON_H_INCLUDED

#include "pybind11/pybind11.h"

namespace Kratos
{

namespace Python
{

  void AddResponseFunctionsToPython(pybind11::module& m);

} // namespace Python

} // namespace Kratos

#endif // KRATOS_RESPONSE_FUNCTIONS_PYTHON_H_INCLUDED
