//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//
//

#if !defined(KRATOS_ADD_CUSTOM_RANS_PROCESSES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_CUSTOM_RANS_PROCESSES_TO_PYTHON_H_INCLUDED

// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/define.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_RANS_PROCESSES_TO_PYTHON_H_INCLUDED defined
