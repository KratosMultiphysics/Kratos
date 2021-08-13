//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#ifndef KRATOS_HDF5APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H
#define KRATOS_HDF5APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/define.h"


namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& rModule);

} // namespace Python
} // namespace Kratos


#endif // KRATOS_HDF5APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H