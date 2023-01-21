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

#pragma once

// --- External Includes ---
#include "pybind11/pybind11.h"


namespace Kratos::Python
{


void AddPredicatesToPython(pybind11::module& rModule);


} // namespace Kratos::Python
