//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// External includes
#include <pybind11/pybind11.h>


namespace Kratos::Python {


void AddCustomAdjointsToPython(pybind11::module& rModule);


} // namespace Kratos::Python
