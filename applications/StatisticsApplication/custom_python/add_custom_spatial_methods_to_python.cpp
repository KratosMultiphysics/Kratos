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
#include <pybind11/stl.h>

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

template <class T>
struct VariantPointer {};

template<class... TArgs>
struct VariantPointer<std::variant<TArgs...>> {
    using type = std::variant<typename TArgs::Pointer...>;
};

template<class TDataType>
void AddSpatialMethodsForVariableType(pybind11::module& m)
{
    using DataLocation = Globals::DataLocation;

    using NormType = typename Norms::NormType<TDataType>::type;

    using VariantPointer = typename VariantPointer<NormType>::type;

    namespace py = pybind11;

    std::string distribution_info_name;
    if constexpr(std::is_same_v<TDataType, double>) {
        distribution_info_name = "DoubleDistributionInfo";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        distribution_info_name = "Array3DistributionInfo";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 4>>) {
        distribution_info_name = "Array4DistributionInfo";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 6>>) {
        distribution_info_name = "Array6DistributionInfo";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 9>>) {
        distribution_info_name = "Array9DistributionInfo";
    } else if constexpr(std::is_same_v<TDataType, Vector>) {
        distribution_info_name = "VectorDistributionInfo";
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        distribution_info_name = "MatrixDistributionInfo";
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
    }

    using distribution_info_type = SpatialMethods::DistributionInfo<TDataType>;
    py::class_<distribution_info_type, typename distribution_info_type::Pointer>(m, distribution_info_name.c_str())
        .def("GetMin", &distribution_info_type::GetMin)
        .def("GetMax", &distribution_info_type::GetMax)
        .def("GetGroupUpperValues", &distribution_info_type::GetGroupUpperValues)
        .def("GetGroupNumberOfValues", &distribution_info_type::GetGroupNumberOfValues)
        .def("GetGroupValueDistributionPercentage", &distribution_info_type::GetGroupValueDistributionPercentage)
        .def("GetGroupMeans", &distribution_info_type::GetGroupMeans)
        .def("GetGroupVariances", &distribution_info_type::GetGroupVariances)
        ;

    m.def("Sum", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&>(&SpatialMethods::Sum<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"));
    m.def("Sum", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::Sum<TDataType>(rModelPart, rVariable, rDataLocation, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("norm"));

    m.def("Mean", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&>(&SpatialMethods::Mean<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"));
    m.def("Mean", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::Mean<TDataType>(rModelPart, rVariable, rDataLocation, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("norm"));

    m.def("RootMeanSquare", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&>(&SpatialMethods::RootMeanSquare<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"));
    m.def("RootMeanSquare", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::RootMeanSquare<TDataType>(rModelPart, rVariable, rDataLocation, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("norm"));

    m.def("Variance", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&>(&SpatialMethods::Variance<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"));
    m.def("Variance", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::Variance<TDataType>(rModelPart, rVariable, rDataLocation, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("norm"));

    m.def("Min", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&>(&SpatialMethods::Min<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"));
    m.def("Min", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::Min<TDataType>(rModelPart, rVariable, rDataLocation, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("norm"));

    m.def("Max", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&>(&SpatialMethods::Max<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"));
    m.def("Max", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::Max<TDataType>(rModelPart, rVariable, rDataLocation, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("norm"));

    m.def("Median", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&>(&SpatialMethods::Median<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"));
    m.def("Median", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::Median<TDataType>(rModelPart, rVariable, rDataLocation, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("norm"));

    m.def("Distribution", py::overload_cast<const ModelPart&,const Variable<TDataType>&,const DataLocation&, Parameters>(&SpatialMethods::Distribution<TDataType>), py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("parameters"));
    m.def("Distribution", [](const ModelPart& rModelPart, const Variable<TDataType>& rVariable, const DataLocation& rDataLocation, Parameters Params, const VariantPointer& rNormPointer) { return std::visit([&](auto& pNorm) { return SpatialMethods::Distribution<TDataType>(rModelPart, rVariable, rDataLocation, Params, *pNorm); }, rNormPointer);}, py::arg("model_part"), py::arg("variable"), py::arg("data_location"), py::arg("parameters"), py::arg("norm"));
}

void AddCustomSpatialMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using DataLocation = Globals::DataLocation;

    auto spatial_method_module = m.def_submodule("SpatialMethods");
    spatial_method_module.def("Sum", [](const ModelPart& rModelPart, const Flags& rFlag, const DataLocation& rLocation) { return SpatialMethods::Sum(rModelPart, rFlag, rLocation); }, py::arg("model_part"), py::arg("flag"), py::arg("data_location"));
    AddSpatialMethodsForVariableType<double>(spatial_method_module);
    AddSpatialMethodsForVariableType<array_1d<double, 3>>(spatial_method_module);
    AddSpatialMethodsForVariableType<array_1d<double, 4>>(spatial_method_module);
    AddSpatialMethodsForVariableType<array_1d<double, 6>>(spatial_method_module);
    AddSpatialMethodsForVariableType<array_1d<double, 9>>(spatial_method_module);
    AddSpatialMethodsForVariableType<Vector>(spatial_method_module);
    AddSpatialMethodsForVariableType<Matrix>(spatial_method_module);

    spatial_method_module.def("Sum", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Sum(*pExpression, rDataCommunicator); });
    spatial_method_module.def("Sum", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Sum(*pExpression, rDataCommunicator, *pNorm); }, rNorm); });
    spatial_method_module.def("Sum", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Sum(*pContainerExpression); }, rContainerExpression); });
    spatial_method_module.def("Sum", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Sum(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); });

    spatial_method_module.def("Mean", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Mean(*pExpression, rDataCommunicator); });
    spatial_method_module.def("Mean", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Mean(*pExpression, rDataCommunicator, *pNorm); }, rNorm); });
    spatial_method_module.def("Mean", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Mean(*pContainerExpression); }, rContainerExpression); });
    spatial_method_module.def("Mean", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Mean(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); });

    spatial_method_module.def("RootMeanSquare", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::RootMeanSquare(*pExpression, rDataCommunicator); });
    spatial_method_module.def("RootMeanSquare", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::RootMeanSquare(*pExpression, rDataCommunicator, *pNorm); }, rNorm); });
    spatial_method_module.def("RootMeanSquare", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::RootMeanSquare(*pContainerExpression); }, rContainerExpression); });
    spatial_method_module.def("RootMeanSquare", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::RootMeanSquare(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); });

    spatial_method_module.def("Variance", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Variance(*pExpression, rDataCommunicator); });
    spatial_method_module.def("Variance", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Variance(*pExpression, rDataCommunicator, *pNorm); }, rNorm); });
    spatial_method_module.def("Variance", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Variance(*pContainerExpression); }, rContainerExpression); });
    spatial_method_module.def("Variance", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Variance(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); });

    spatial_method_module.def("Min", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Min(*pExpression, rDataCommunicator); });
    spatial_method_module.def("Min", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Min(*pExpression, rDataCommunicator, *pNorm); }, rNorm); });
    spatial_method_module.def("Min", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Min(*pContainerExpression); }, rContainerExpression); });
    spatial_method_module.def("Min", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Min(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); });

    spatial_method_module.def("Max", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Max(*pExpression, rDataCommunicator); });
    spatial_method_module.def("Max", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Max(*pExpression, rDataCommunicator, *pNorm); }, rNorm); });
    spatial_method_module.def("Max", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Max(*pContainerExpression); }, rContainerExpression); });
    spatial_method_module.def("Max", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Max(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); });

    spatial_method_module.def("Median", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Median(*pExpression, rDataCommunicator); });
    spatial_method_module.def("Median", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Median(*pExpression, rDataCommunicator, *pNorm); }, rNorm); });
    spatial_method_module.def("Median", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Median(*pContainerExpression); }, rContainerExpression); });
    spatial_method_module.def("Median", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Median(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); });

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
