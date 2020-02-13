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
#include "custom_methods/spatial_methods.h"

// Include base h
#include "custom_python/add_custom_spatial_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomSpatialMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto spatial_method_module = m.def_submodule("SpatialMethods");

    auto spatial_historical_method_module = spatial_method_module.def_submodule("Historical");
    spatial_historical_method_module.def_submodule("ValueMethods");
    spatial_historical_method_module.def_submodule("NormMethods");

    auto spatial_non_historical_method_module = spatial_method_module.def_submodule("NonHistorical");
    auto spatial_non_historical_nodal_method_module = spatial_non_historical_method_module.def_submodule("Nodes");
    spatial_non_historical_nodal_method_module.def_submodule("ValueMethods");
    spatial_non_historical_nodal_method_module.def_submodule("NormMethods");
    auto spatial_non_historical_condition_method_module = spatial_non_historical_method_module.def_submodule("Conditions");
    spatial_non_historical_condition_method_module.def_submodule("ValueMethods");
    spatial_non_historical_condition_method_module.def_submodule("NormMethods");
    auto spatial_non_historical_element_method_module = spatial_non_historical_method_module.def_submodule("Elements");
    spatial_non_historical_element_method_module.def_submodule("ValueMethods");
    spatial_non_historical_element_method_module.def_submodule("NormMethods");

    ADD_KRATOS_STATISTICS_SPATIAL_VALUE_METHOD_PYTHON_INTERFACE(CalculateSum, "Sum", m)
    ADD_KRATOS_STATISTICS_SPATIAL_VALUE_METHOD_PYTHON_INTERFACE(CalculateMean, "Mean", m)
    ADD_KRATOS_STATISTICS_SPATIAL_VALUE_METHOD_PYTHON_INTERFACE(CalculateVariance, "Variance", m)
    ADD_KRATOS_STATISTICS_SPATIAL_VALUE_METHOD_PYTHON_INTERFACE(CalculateRootMeanSquare, "RootMeanSquare", m)

    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(CalculateNormSum, "Sum", m)
    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(CalculateNormMean, "Mean", m)
    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(CalculateNormVariance, "Variance", m)
    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(CalculateNormRootMeanSquare, "RootMeanSquare", m)
    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(GetNormMin, "Min", m)
    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(GetNormMax, "Max", m)
    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(GetNormMedian, "Median", m)
    ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(GetNormDistribution, "Distribution", m)
}

} // namespace Python.
} // Namespace Kratos
