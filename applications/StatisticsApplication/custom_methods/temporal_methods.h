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

#if !defined(KRATOS_TEMPORAL_METHODS_H_INCLUDED)
#define KRATOS_TEMPORAL_METHODS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

// Application method includes
#include "custom_methods/temporal_max_method.h"
#include "custom_methods/temporal_mean_method.h"
#include "custom_methods/temporal_min_method.h"
#include "custom_methods/temporal_sum_method.h"
#include "custom_methods/temporal_variance_method.h"

// macros
#ifndef KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE
#define KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE(type) \
const Variable<type>&
#endif

#ifndef KRATOS_STATISTICS_TWO_OUTPUT_VALUE_TYPE
#define KRATOS_STATISTICS_TWO_OUTPUT_VALUE_TYPE(type) \
const Variable<type>&, const Variable<type>&
#endif

#ifndef KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE
#define KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE(type) \
const Variable<double>&
#endif

#ifndef KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE
#define KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE(type) \
const Variable<double>&, const Variable<double>&
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(method, method_name, value_method_module, container, output_type)  \
{                                                                                                                                       \
    using type_double = double;                                                                                                         \
    using type_array = array_1d<double, 3>;                                                                                             \
    using type_vector = Vector;                                                                                                         \
    using type_matrix = Matrix;                                                                                                         \
    auto method_module = value_method_module.def_submodule(method_name);                                                                \
    using current_method = TemporalMethods::container::method;                                                                          \
    py::class_<current_method::ValueMethod<type_double>, current_method::ValueMethod<type_double>::Pointer,                             \
                TemporalMethods::TemporalMethod>(method_module, "Double")                                                               \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,                                                     \
                        output_type(type_double)>());                                                                                   \
    py::class_<current_method::ValueMethod<type_array>, current_method::ValueMethod<type_array>::Pointer,                               \
                TemporalMethods::TemporalMethod>(method_module, "Array")                                                                \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,                                                      \
                        output_type(type_array)>());                                                                                    \
    py::class_<current_method::ValueMethod<type_vector>, current_method::ValueMethod<type_vector>::Pointer,                             \
                TemporalMethods::TemporalMethod>(method_module, "Vector")                                                               \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,                                                     \
                        output_type(type_vector)>());                                                                                   \
    py::class_<current_method::ValueMethod<type_matrix>, current_method::ValueMethod<type_matrix>::Pointer,                             \
                TemporalMethods::TemporalMethod>(method_module, "Matrix")                                                               \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,                                                     \
                        output_type(type_matrix)>());                                                                                   \
}
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(method, method_name, norm_method_module, container, output_type)    \
{                                                                                                                                       \
    using type_double = double;                                                                                                         \
    using type_array = array_1d<double, 3>;                                                                                             \
    using type_vector = Vector;                                                                                                         \
    using type_matrix = Matrix;                                                                                                         \
    auto method_module = norm_method_module.def_submodule(method_name);                                                                 \
    using current_method = TemporalMethods::container::method;                                                                          \
    py::class_<current_method::NormMethod<type_double>, current_method::NormMethod<type_double>::Pointer,                               \
                TemporalMethods::TemporalMethod>(method_module, "Double")                                                               \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,                                                     \
                        output_type(type_double)>());                                                                                                \
    py::class_<current_method::NormMethod<type_array>, current_method::NormMethod<type_array>::Pointer,                                 \
                TemporalMethods::TemporalMethod>(method_module, "Array")                                                                \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,                                                      \
                        output_type(type_array)>());                                                                                                \
    py::class_<current_method::NormMethod<type_vector>, current_method::NormMethod<type_vector>::Pointer,                               \
                TemporalMethods::TemporalMethod>(method_module, "Vector")                                                               \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,                                                     \
                        output_type(type_vector)>());                                                                                                \
    py::class_<current_method::NormMethod<type_matrix>, current_method::NormMethod<type_matrix>::Pointer,                               \
                TemporalMethods::TemporalMethod>(method_module, "Matrix")                                                               \
        .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,                                                     \
                        output_type(type_matrix)>());                                                                                                \
}
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(method, method_name, python_application_module) \
{                                                                                                               \
    py::module temporal_module = (py::module) python_application_module.attr("TemporalMethods"); \
    py::module historical_module = (py::module) temporal_module.attr("Historical"); \
    py::module historical_historical_output_module = (py::module) historical_module.attr("HistoricalOutput"); \
    py::module historical_non_historical_output_module = (py::module) historical_module.attr("NonHistoricalOutput"); \
    py::module non_historical_module = (py::module) temporal_module.attr("NonHistorical"); \
    py::module non_historical_nodal_module = (py::module) non_historical_module.attr("Nodes"); \
    py::module non_historical_condition_module = (py::module) non_historical_module.attr("Conditions"); \
    py::module non_historical_element_module = (py::module) non_historical_module.attr("Elements"); \
    historical_historical_output_module.def(method_name, &TemporalMethods::HistoricalInputHistoricalOutputTemporalMethods::method::CreateTemporalMethodObject); \
    historical_non_historical_output_module.def(method_name, &TemporalMethods::HistoricalInputNonHistoricalOutputTemporalMethods::method::CreateTemporalMethodObject); \
    non_historical_nodal_module.def(method_name, &TemporalMethods::NodalNonHistoricalTemporalMethods::method::CreateTemporalMethodObject); \
    non_historical_condition_module.def(method_name, &TemporalMethods::ConditionNonHistoricalTemporalMethods::method::CreateTemporalMethodObject); \
    non_historical_element_module.def(method_name, &TemporalMethods::ElementNonHistoricalTemporalMethods::method::CreateTemporalMethodObject); \
}
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(method, method_name, python_application_module, output_type)       \
{                                                                                                                                       \
    py::module temporal_module = (py::module) python_application_module.attr("TemporalMethods"); \
    py::module historical_module = (py::module) temporal_module.attr("Historical"); \
    py::module historical_historical_output_module = (py::module) historical_module.attr("HistoricalOutput"); \
    py::module historical_historical_output_value_module = (py::module) historical_historical_output_module.attr("ValueMethods"); \
    py::module historical_non_historical_output_module = (py::module) historical_module.attr("NonHistoricalOutput"); \
    py::module historical_non_historical_output_value_module = (py::module) historical_non_historical_output_module.attr("ValueMethods"); \
    py::module non_historical_module = (py::module) temporal_module.attr("NonHistorical"); \
    py::module non_historical_nodal_module = (py::module) non_historical_module.attr("Nodes"); \
    py::module non_historical_nodal_value_module = (py::module) non_historical_nodal_module.attr("ValueMethods"); \
    py::module non_historical_condition_module = (py::module) non_historical_module.attr("Conditions"); \
    py::module non_historical_condition_value_module = (py::module) non_historical_condition_module.attr("ValueMethods"); \
    py::module non_historical_element_module = (py::module) non_historical_module.attr("Elements"); \
    py::module non_historical_element_value_module = (py::module) non_historical_element_module.attr("ValueMethods"); \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(method, method_name, historical_historical_output_value_module, HistoricalInputHistoricalOutputTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(method, method_name, historical_non_historical_output_value_module, HistoricalInputNonHistoricalOutputTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(method, method_name, non_historical_nodal_value_module, NodalNonHistoricalTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(method, method_name, non_historical_condition_value_module, ConditionNonHistoricalTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_VALUE_METHOD_PYTHON_INTERFACE(method, method_name, non_historical_element_value_module, ElementNonHistoricalTemporalMethods, output_type) \
}
#endif

#ifndef ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE
#define ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(method, method_name, python_application_module, output_type)       \
{                                                                                                                                       \
    py::module temporal_module = (py::module) python_application_module.attr("TemporalMethods"); \
    py::module historical_module = (py::module) temporal_module.attr("Historical"); \
    py::module historical_historical_output_module = (py::module) historical_module.attr("HistoricalOutput"); \
    py::module historical_historical_output_value_module = (py::module) historical_historical_output_module.attr("NormMethods"); \
    py::module historical_non_historical_output_module = (py::module) historical_module.attr("NonHistoricalOutput"); \
    py::module historical_non_historical_output_value_module = (py::module) historical_non_historical_output_module.attr("NormMethods"); \
    py::module non_historical_module = (py::module) temporal_module.attr("NonHistorical"); \
    py::module non_historical_nodal_module = (py::module) non_historical_module.attr("Nodes"); \
    py::module non_historical_nodal_value_module = (py::module) non_historical_nodal_module.attr("NormMethods"); \
    py::module non_historical_condition_module = (py::module) non_historical_module.attr("Conditions"); \
    py::module non_historical_condition_value_module = (py::module) non_historical_condition_module.attr("NormMethods"); \
    py::module non_historical_element_module = (py::module) non_historical_module.attr("Elements"); \
    py::module non_historical_element_value_module = (py::module) non_historical_element_module.attr("NormMethods"); \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(method, method_name, historical_historical_output_value_module, HistoricalInputHistoricalOutputTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(method, method_name, historical_non_historical_output_value_module, HistoricalInputNonHistoricalOutputTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(method, method_name, non_historical_nodal_value_module, NodalNonHistoricalTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(method, method_name, non_historical_condition_value_module, ConditionNonHistoricalTemporalMethods, output_type) \
    ADD_KRATOS_STATISTICS_TEMPORAL_CONTAINER_NORM_METHOD_PYTHON_INTERFACE(method, method_name, non_historical_element_value_module, ElementNonHistoricalTemporalMethods, output_type) \
}
#endif

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

namespace TemporalMethods
{
template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor, template <typename T> typename TDataStorageFunctor>
class ContainerTemporalMethods
{
public:
    using SumMethod = TemporalSumMethod<TContainerType, TContainerItemType, TDataRetrievalFunctor, TDataStorageFunctor>;
    using MeanMethod = TemporalMeanMethod<TContainerType, TContainerItemType, TDataRetrievalFunctor, TDataStorageFunctor>;
    using VarianceMethod = TemporalVarianceMethod<TContainerType, TContainerItemType, TDataRetrievalFunctor, TDataStorageFunctor>;
    using MinMethod = TemporalMinMethod<TContainerType, TContainerItemType, TDataRetrievalFunctor, TDataStorageFunctor>;
    using MaxMethod = TemporalMaxMethod<TContainerType, TContainerItemType, TDataRetrievalFunctor, TDataStorageFunctor>;
};

using NodeType = ModelPart::NodeType;
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;

using NodesContainerType = ModelPart::NodesContainerType;
using ElementsContainerType = ModelPart::ElementsContainerType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;

template <template <typename T> typename TDataStorageFunctor>
class HistoricalTemporalMethods
    : public ContainerTemporalMethods<NodesContainerType, NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor, TDataStorageFunctor>
{
};

class HistoricalInputHistoricalOutputTemporalMethods
    : public HistoricalTemporalMethods<MethodsUtilities::HistoricalDataValueRetrievalFunctor>
{
};

class HistoricalInputNonHistoricalOutputTemporalMethods
    : public HistoricalTemporalMethods<MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class NodalNonHistoricalTemporalMethods
    : public ContainerTemporalMethods<NodesContainerType, NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class ConditionNonHistoricalTemporalMethods
    : public ContainerTemporalMethods<ConditionsContainerType, ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class ElementNonHistoricalTemporalMethods
    : public ContainerTemporalMethods<ElementsContainerType, ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

} // namespace TemporalMethods
} // namespace Kratos
#endif // KRATOS_TEMPORAL_METHODS_H_INCLUDED