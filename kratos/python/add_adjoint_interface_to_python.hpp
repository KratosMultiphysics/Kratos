//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// --- External Includes ---
#include <pybind11/pybind11.h>


namespace Kratos::Python {


void AddAdjointInterfaceToPython(pybind11::module_& rModule);


} // namespace Kratos::Python
