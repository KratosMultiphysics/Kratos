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
#include "custom_methods/temporal_methods.h"

// Include base h
#include "custom_python/add_custom_temporal_methods_to_python.h"

// macros
#ifndef KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE
#define KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE(type) const Variable<type>&
#endif

#ifndef KRATOS_STATISTICS_TWO_OUTPUT_VALUE_TYPE
#define KRATOS_STATISTICS_TWO_OUTPUT_VALUE_TYPE(type) \
    const Variable<type>&, const Variable<type>&
#endif

#ifndef KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE
#define KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE(type) const Variable<double>&
#endif

#ifndef KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE
#define KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE(type) \
    const Variable<double>&, const Variable<double>&
#endif

#ifndef KRATOS_STATISTICS_DEFAULT_INPUTS
#define KRATOS_STATISTICS_DEFAULT_INPUTS(type) \
    ModelPart&, const std::string&, const Variable<type>&, const int
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(                         \
    method, method_name, value_method_module, container, output_type)                                   \
    {                                                                                                   \
        using type_double = double;                                                                     \
        using type_array = array_1d<double, 3>;                                                         \
        using type_vector = Vector;                                                                     \
        using type_matrix = Matrix;                                                                     \
        auto method_module = value_method_module.def_submodule(method_name);                            \
        using current_method = TemporalMethods::container::method;                                      \
        py::class_<current_method::ValueMethod<type_double>,                                            \
                   current_method::ValueMethod<type_double>::Pointer, TemporalMethods::TemporalMethod>( \
            method_module, "Double")                                                                    \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_double), output_type(type_double)>());  \
        py::class_<current_method::ValueMethod<type_array>,                                             \
                   current_method::ValueMethod<type_array>::Pointer, TemporalMethods::TemporalMethod>(  \
            method_module, "Array")                                                                     \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_array), output_type(type_array)>());    \
        py::class_<current_method::ValueMethod<type_vector>,                                            \
                   current_method::ValueMethod<type_vector>::Pointer, TemporalMethods::TemporalMethod>( \
            method_module, "Vector")                                                                    \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_vector), output_type(type_vector)>());  \
        py::class_<current_method::ValueMethod<type_matrix>,                                            \
                   current_method::ValueMethod<type_matrix>::Pointer, TemporalMethods::TemporalMethod>( \
            method_module, "Matrix")                                                                    \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_matrix), output_type(type_matrix)>());  \
    }
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(                              \
    method, method_name, norm_method_module, container, output_type)                                        \
    {                                                                                                       \
        using type_double = double;                                                                         \
        using type_array = array_1d<double, 3>;                                                             \
        using type_vector = Vector;                                                                         \
        using type_matrix = Matrix;                                                                         \
        auto method_module = norm_method_module.def_submodule(method_name);                                 \
        using current_method = TemporalMethods::container::method;                                          \
        py::class_<current_method::NormMethod<type_double>,                                                 \
                   current_method::NormMethod<type_double>::Pointer, TemporalMethods::TemporalMethod>(      \
            method_module, "Double")                                                                        \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_double), output_type(type_double)>());      \
        py::class_<current_method::NormMethod<type_array>, current_method::NormMethod<type_array>::Pointer, \
                   TemporalMethods::TemporalMethod>(method_module, "Array")                                 \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_array), output_type(type_array)>());        \
        py::class_<current_method::NormMethod<type_vector>,                                                 \
                   current_method::NormMethod<type_vector>::Pointer, TemporalMethods::TemporalMethod>(      \
            method_module, "Vector")                                                                        \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_vector), output_type(type_vector)>());      \
        py::class_<current_method::NormMethod<type_matrix>,                                                 \
                   current_method::NormMethod<type_matrix>::Pointer, TemporalMethods::TemporalMethod>(      \
            method_module, "Matrix")                                                                        \
            .def(py::init<KRATOS_STATISTICS_DEFAULT_INPUTS(type_matrix), output_type(type_matrix)>());      \
    }
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(                                                \
    method, method_name, python_application_module)                                                                   \
    {                                                                                                                 \
        py::module temporal_module =                                                                                  \
            (py::module)python_application_module.attr("TemporalMethods");                                            \
        py::module historical_module =                                                                                \
            (py::module)temporal_module.attr("Historical");                                                           \
        py::module historical_historical_output_module =                                                              \
            (py::module)historical_module.attr("HistoricalOutput");                                                   \
        py::module historical_non_historical_output_module =                                                          \
            (py::module)historical_module.attr("NonHistoricalOutput");                                                \
        py::module non_historical_module =                                                                            \
            (py::module)temporal_module.attr("NonHistorical");                                                        \
        py::module non_historical_nodal_module =                                                                      \
            (py::module)non_historical_module.attr("Nodes");                                                          \
        py::module non_historical_condition_module =                                                                  \
            (py::module)non_historical_module.attr("Conditions");                                                     \
        py::module non_historical_element_module =                                                                    \
            (py::module)non_historical_module.attr("Elements");                                                       \
        historical_historical_output_module.def(                                                                      \
            method_name,                                                                                              \
            &TemporalMethods::HistoricalInputHistoricalOutputTemporalMethods::method::CreateTemporalMethodObject);    \
        historical_non_historical_output_module.def(                                                                  \
            method_name,                                                                                              \
            &TemporalMethods::HistoricalInputNonHistoricalOutputTemporalMethods::method::CreateTemporalMethodObject); \
        non_historical_nodal_module.def(                                                                              \
            method_name, &TemporalMethods::NodalNonHistoricalTemporalMethods::method::CreateTemporalMethodObject);    \
        non_historical_condition_module.def(                                                                          \
            method_name,                                                                                              \
            &TemporalMethods::ConditionNonHistoricalTemporalMethods::method::CreateTemporalMethodObject);             \
        non_historical_element_module.def(                                                                            \
            method_name,                                                                                              \
            &TemporalMethods::ElementNonHistoricalTemporalMethods::method::CreateTemporalMethodObject);               \
    }
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(           \
    method, method_name, python_application_module, output_type)                \
    {                                                                           \
        py::module temporal_module =                                            \
            (py::module)python_application_module.attr("TemporalMethods");      \
        py::module historical_module =                                          \
            (py::module)temporal_module.attr("Historical");                     \
        py::module historical_historical_output_module =                        \
            (py::module)historical_module.attr("HistoricalOutput");             \
        py::module historical_historical_output_value_module =                  \
            (py::module)historical_historical_output_module.attr(               \
                "ValueMethods");                                                \
        py::module historical_non_historical_output_module =                    \
            (py::module)historical_module.attr("NonHistoricalOutput");          \
        py::module historical_non_historical_output_value_module =              \
            (py::module)historical_non_historical_output_module.attr(           \
                "ValueMethods");                                                \
        py::module non_historical_module =                                      \
            (py::module)temporal_module.attr("NonHistorical");                  \
        py::module non_historical_nodal_module =                                \
            (py::module)non_historical_module.attr("Nodes");                    \
        py::module non_historical_nodal_value_module =                          \
            (py::module)non_historical_nodal_module.attr("ValueMethods");       \
        py::module non_historical_condition_module =                            \
            (py::module)non_historical_module.attr("Conditions");               \
        py::module non_historical_condition_value_module =                      \
            (py::module)non_historical_condition_module.attr("ValueMethods");   \
        py::module non_historical_element_module =                              \
            (py::module)non_historical_module.attr("Elements");                 \
        py::module non_historical_element_value_module =                        \
            (py::module)non_historical_element_module.attr("ValueMethods");     \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE( \
            method, method_name, historical_historical_output_value_module,     \
            HistoricalInputHistoricalOutputTemporalMethods, output_type)        \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE( \
            method, method_name, historical_non_historical_output_value_module, \
            HistoricalInputNonHistoricalOutputTemporalMethods, output_type)     \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE( \
            method, method_name, non_historical_nodal_value_module,             \
            NodalNonHistoricalTemporalMethods, output_type)                     \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE( \
            method, method_name, non_historical_condition_value_module,         \
            ConditionNonHistoricalTemporalMethods, output_type)                 \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE( \
            method, method_name, non_historical_element_value_module,           \
            ElementNonHistoricalTemporalMethods, output_type)                   \
    }
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(            \
    method, method_name, python_application_module, output_type)                \
    {                                                                           \
        py::module temporal_module =                                            \
            (py::module)python_application_module.attr("TemporalMethods");      \
        py::module historical_module =                                          \
            (py::module)temporal_module.attr("Historical");                     \
        py::module historical_historical_output_module =                        \
            (py::module)historical_module.attr("HistoricalOutput");             \
        py::module historical_historical_output_value_module =                  \
            (py::module)historical_historical_output_module.attr(               \
                "NormMethods");                                                 \
        py::module historical_non_historical_output_module =                    \
            (py::module)historical_module.attr("NonHistoricalOutput");          \
        py::module historical_non_historical_output_value_module =              \
            (py::module)historical_non_historical_output_module.attr(           \
                "NormMethods");                                                 \
        py::module non_historical_module =                                      \
            (py::module)temporal_module.attr("NonHistorical");                  \
        py::module non_historical_nodal_module =                                \
            (py::module)non_historical_module.attr("Nodes");                    \
        py::module non_historical_nodal_value_module =                          \
            (py::module)non_historical_nodal_module.attr("NormMethods");        \
        py::module non_historical_condition_module =                            \
            (py::module)non_historical_module.attr("Conditions");               \
        py::module non_historical_condition_value_module =                      \
            (py::module)non_historical_condition_module.attr("NormMethods");    \
        py::module non_historical_element_module =                              \
            (py::module)non_historical_module.attr("Elements");                 \
        py::module non_historical_element_value_module =                        \
            (py::module)non_historical_element_module.attr("NormMethods");      \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(  \
            method, method_name, historical_historical_output_value_module,     \
            HistoricalInputHistoricalOutputTemporalMethods, output_type)        \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(  \
            method, method_name, historical_non_historical_output_value_module, \
            HistoricalInputNonHistoricalOutputTemporalMethods, output_type)     \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(  \
            method, method_name, non_historical_nodal_value_module,             \
            NodalNonHistoricalTemporalMethods, output_type)                     \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(  \
            method, method_name, non_historical_condition_value_module,         \
            ConditionNonHistoricalTemporalMethods, output_type)                 \
        ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(  \
            method, method_name, non_historical_element_value_module,           \
            ElementNonHistoricalTemporalMethods, output_type)                   \
    }
#endif

namespace Kratos
{
namespace Python
{
void AddCustomTemporalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto temporal_method_module = m.def_submodule("TemporalMethods");

    py::class_<TemporalMethods::TemporalMethod, TemporalMethods::TemporalMethod::Pointer>(temporal_method_module,"TemporalMethod")
        .def(py::init<ModelPart&, const int>())
        .def("GetModelPart", &TemporalMethods::TemporalMethod::GetModelPart)
        .def("GetTotalTime", &TemporalMethods::TemporalMethod::GetTotalTime)
        .def("GetEchoLevel", &TemporalMethods::TemporalMethod::GetEchoLevel)
        .def("InitializeStatisticsMethod", &TemporalMethods::TemporalMethod::InitializeStatisticsMethod)
        .def("CalculateStatistics", &TemporalMethods::TemporalMethod::CalculateStatistics);

    auto temporal_historical_method = temporal_method_module.def_submodule("Historical");
    auto temporal_historical_historical_method = temporal_historical_method.def_submodule("HistoricalOutput");
    temporal_historical_historical_method.def_submodule("ValueMethods");
    temporal_historical_historical_method.def_submodule("NormMethods");
    auto temporal_historical_non_historical_method = temporal_historical_method.def_submodule("NonHistoricalOutput");
    temporal_historical_non_historical_method.def_submodule("ValueMethods");
    temporal_historical_non_historical_method.def_submodule("NormMethods");
    auto temporal_non_historical_method = temporal_method_module.def_submodule("NonHistorical");
    auto temporal_non_historical_nodal_method = temporal_non_historical_method.def_submodule("Nodes");
    temporal_non_historical_nodal_method.def_submodule("ValueMethods");
    temporal_non_historical_nodal_method.def_submodule("NormMethods");
    auto temporal_non_historical_condition_method = temporal_non_historical_method.def_submodule("Conditions");
    temporal_non_historical_condition_method.def_submodule("ValueMethods");
    temporal_non_historical_condition_method.def_submodule("NormMethods");
    auto temporal_non_historical_element_method = temporal_non_historical_method.def_submodule("Elements");
    temporal_non_historical_element_method.def_submodule("ValueMethods");
    temporal_non_historical_element_method.def_submodule("NormMethods");

    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(SumMethod, "Sum", m, KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(MeanMethod, "Mean", m, KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(VarianceMethod, "Variance", m, KRATOS_STATISTICS_TWO_OUTPUT_VALUE_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(RootMeanSquareMethod, "RootMeanSquare", m, KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE)

    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(SumMethod, "Sum", m, KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(MeanMethod, "Mean", m, KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(VarianceMethod, "Variance", m, KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(RootMeanSquareMethod, "RootMeanSquare", m, KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(MinMethod, "Min", m, KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(MaxMethod, "Max", m, KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE)

    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(SumMethod, "Sum", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(MeanMethod, "Mean", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(VarianceMethod, "Variance", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(RootMeanSquareMethod, "RootMeanSquare", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(MinMethod, "Min", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(MaxMethod, "Max", m)
}

} // namespace Python.
} // Namespace Kratos
