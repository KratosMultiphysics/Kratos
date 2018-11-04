//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//
//


#if !defined(KRATOS_STRATEGIES_PYTHON_H_INCLUDED )
#define  KRATOS_STRATEGIES_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

namespace Python
{

void  AddCustomStrategiesToPython(pybind11::module& m);
/* { */
/*   using namespace boost::python; */
/* } */

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_STRATEGIES_TO_PYTHON_H_INCLUDED  defined 
