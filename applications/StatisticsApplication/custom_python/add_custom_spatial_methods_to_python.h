//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

#if !defined(KRATOS_ADD_SPATIAL_METHODS_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_SPATIAL_METHODS_TO_PYTHON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

// Application includes
#include "custom_python/add_custom_spatial_historical_methods_to_python.h"
#include "custom_python/add_custom_spatial_non_historical_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomSpatialMethodsToPython(pybind11::module& m)
{
    auto spatial_method_module = m.def_submodule("SpatialMethods");
    AddCustomSpatialHistoricalMethodsToPython(spatial_method_module);
    AddCustomSpatialNonHistoricalMethodsToPython(spatial_method_module);
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_ADD_SPATIAL_METHODS_TO_PYTHON_H_INCLUDED  defined
