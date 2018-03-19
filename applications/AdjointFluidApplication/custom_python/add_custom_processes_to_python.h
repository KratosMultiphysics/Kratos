//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_CUSTOM_PROCESSES_PYTHON_H_INCLUDED )
#define  KRATOS_CUSTOM_PROCESSES_PYTHON_H_INCLUDED

// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m);

}  /* namespace Python */

}  /* namespace Kratos */

#endif /* KRATOS_CUSTOM_PROCESSES_PYTHON_H_INCLUDED defined */
