//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

namespace Kratos::Python {

/**
 * @brief Register KaHIPApplication modeler classes with the Python module.
 * @details Exposes:
 *          - KaHIPPartitioningModeler
 * @param m  The pybind11 module object
 */
void AddCustomModelerToPython(pybind11::module& m);

} // namespace Kratos::Python
