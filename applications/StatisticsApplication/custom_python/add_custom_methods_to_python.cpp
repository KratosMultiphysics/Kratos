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

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes

// Application includes

// spatial methods
#include "custom_python/add_custom_spatial_historical_nodal_methods_to_python.h"
#include "custom_python/add_custom_spatial_non_historical_nodal_methods_to_python.h"
#include "custom_python/add_custom_spatial_non_historical_condition_methods_to_python.h"
#include "custom_python/add_custom_spatial_non_historical_element_methods_to_python.h"

// temporal methods
#include "custom_python/add_custom_temporal_methods_to_python.h"

// Include base h
#include "custom_python/add_custom_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto spatial_method_module = m.def_submodule("SpatialMethods");

    auto spatial_historical_method_module = spatial_method_module.def_submodule("Historical");
    AddCustomSpatialHistoricalNodalMethodsToPython(spatial_historical_method_module);

    auto spatial_non_historical_method_module = spatial_method_module.def_submodule("NonHistorical");
    AddCustomSpatialNonHistoricalNodalMethodsToPython(spatial_non_historical_method_module);
    AddCustomSpatialNonHistoricalConditionMethodsToPython(spatial_non_historical_method_module);
    AddCustomSpatialNonHistoricalElementMethodsToPython(spatial_non_historical_method_module);
}

} // namespace Python.
} // Namespace Kratos
