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
