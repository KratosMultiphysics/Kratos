//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes

// Application includes
#include "custom_python/add_custom_spatial_methods_to_python.h"
#include "custom_python/add_custom_temporal_methods_to_python.h"

// Include base h
#include "custom_python/add_custom_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomMethodsToPython(pybind11::module& m)
{
    AddCustomSpatialMethodsToPython(m);
    AddCustomTemporalMethodsToPython(m);
}

} // namespace Python.
} // Namespace Kratos
