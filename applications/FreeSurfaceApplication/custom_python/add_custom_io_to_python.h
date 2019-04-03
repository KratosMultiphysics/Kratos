//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//

#if !defined(KRATOS_ADD_CUSTOM_IO_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_IO_TO_PYTHON_H_INCLUDED



// System includes
#include<pybind11/pybind11.h>

// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

void  AddCustomIOToPython(pybind11::module& pymodule);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_IO_TO_PYTHON_H_INCLUDED  defined 
