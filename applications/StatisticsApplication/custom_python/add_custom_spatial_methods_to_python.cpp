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

// MPI methods doesn't support vector and matrix SumAll methods,
// therefore only double and fixed array templates are added
#ifndef ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(         \
    method, method_name, value_method_module, container)                               \
    {                                                                                  \
        using type_double = double;                                                    \
        using type_array = array_1d<double, 3>;                                        \
        using current_container = SpatialMethods::container;                           \
        value_method_module.def(method_name, &current_container::method<type_double>); \
        value_method_module.def(method_name, &current_container::method<type_array>);  \
    }
#endif

#ifndef ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(        \
    method, method_name, norm_method_module, container)                              \
    {                                                                                \
        using type_double = double;                                                  \
        using type_array = array_1d<double, 3>;                                      \
        using type_vector = Vector;                                                  \
        using type_matrix = Matrix;                                                  \
        using current_container = SpatialMethods::container;                         \
        norm_method_module.def(method_name, &current_container::method<type_double>, \
                               py::arg("model_part"), py::arg("variable"),           \
                               py::arg("norm_type"),                                 \
                               py::arg("parameters") = Parameters(R"({})"));         \
        norm_method_module.def(method_name, &current_container::method<type_array>,  \
                               py::arg("model_part"), py::arg("variable"),           \
                               py::arg("norm_type"),                                 \
                               py::arg("parameters") = Parameters(R"({})"));         \
        norm_method_module.def(method_name, &current_container::method<type_vector>, \
                               py::arg("model_part"), py::arg("variable"),           \
                               py::arg("norm_type"),                                 \
                               py::arg("parameters") = Parameters(R"({})"));         \
        norm_method_module.def(method_name, &current_container::method<type_matrix>, \
                               py::arg("model_part"), py::arg("variable"),           \
                               py::arg("norm_type"),                                 \
                               py::arg("parameters") = Parameters(R"({})"));         \
    }
#endif

#ifndef ADD_KRATOS_STATISTICS_SPATIAL_VALUE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_SPATIAL_VALUE_METHOD_PYTHON_INTERFACE(                                      \
    method, method_name, python_application_module)                                                       \
    {                                                                                                     \
        py::module spatial_module =                                                                       \
            (py::module)python_application_module.attr("SpatialMethods");                                 \
        py::module historical_module =                                                                    \
            (py::module)spatial_module.attr("Historical");                                                \
        py::module historical_value_module =                                                              \
            (py::module)historical_module.attr("ValueMethods");                                           \
        py::module non_historical_module =                                                                \
            (py::module)spatial_module.attr("NonHistorical");                                             \
        py::module non_historical_nodal_module =                                                          \
            (py::module)non_historical_module.attr("Nodes");                                              \
        py::module non_historical_nodal_value_module =                                                    \
            (py::module)non_historical_nodal_module.attr("ValueMethods");                                 \
        py::module non_historical_condition_module =                                                      \
            (py::module)non_historical_module.attr("Conditions");                                         \
        py::module non_historical_condition_value_module =                                                \
            (py::module)non_historical_condition_module.attr("ValueMethods");                             \
        py::module non_historical_element_module =                                                        \
            (py::module)non_historical_module.attr("Elements");                                           \
        py::module non_historical_element_value_module =                                                  \
            (py::module)non_historical_element_module.attr("ValueMethods");                               \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, historical_value_module, HistoricalSpatialMethods)                       \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, non_historical_nodal_value_module, NodalNonHistoricalSpatialMethods)     \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, non_historical_condition_value_module,                                   \
            ConditionNonHistoricalSpatialMethods)                                                         \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, non_historical_element_value_module, ElementNonHistoricalSpatialMethods) \
    }
#endif

#ifndef ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_SPATIAL_NORM_METHOD_PYTHON_INTERFACE(                                      \
    method, method_name, python_application_module)                                                      \
    {                                                                                                    \
        py::module spatial_module =                                                                      \
            (py::module)python_application_module.attr("SpatialMethods");                                \
        py::module historical_module =                                                                   \
            (py::module)spatial_module.attr("Historical");                                               \
        py::module historical_norm_module =                                                              \
            (py::module)historical_module.attr("NormMethods");                                           \
        py::module non_historical_module =                                                               \
            (py::module)spatial_module.attr("NonHistorical");                                            \
        py::module non_historical_nodal_module =                                                         \
            (py::module)non_historical_module.attr("Nodes");                                             \
        py::module non_historical_nodal_norm_module =                                                    \
            (py::module)non_historical_nodal_module.attr("NormMethods");                                 \
        py::module non_historical_condition_module =                                                     \
            (py::module)non_historical_module.attr("Conditions");                                        \
        py::module non_historical_condition_norm_module =                                                \
            (py::module)non_historical_condition_module.attr("NormMethods");                             \
        py::module non_historical_element_module =                                                       \
            (py::module)non_historical_module.attr("Elements");                                          \
        py::module non_historical_element_norm_module =                                                  \
            (py::module)non_historical_element_module.attr("NormMethods");                               \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, historical_norm_module, HistoricalSpatialMethods)                       \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, non_historical_nodal_norm_module, NodalNonHistoricalSpatialMethods)     \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, non_historical_condition_norm_module,                                   \
            ConditionNonHistoricalSpatialMethods)                                                        \
        ADD_KRATOS_STATISTICS_SPATIAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(                            \
            method, method_name, non_historical_element_norm_module, ElementNonHistoricalSpatialMethods) \
    }
#endif

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
