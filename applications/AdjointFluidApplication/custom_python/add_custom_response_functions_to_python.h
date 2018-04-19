//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
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
