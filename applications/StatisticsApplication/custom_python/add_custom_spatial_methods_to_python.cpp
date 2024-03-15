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

    spatial_method_module.def("Sum", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Sum(*pExpression, rDataCommunicator); }, py::arg("expression"), py::arg("data_communicator"));
    spatial_method_module.def("Sum", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Sum(*pExpression, rDataCommunicator, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("norm"));
    spatial_method_module.def("Sum", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Sum(*pContainerExpression); }, rContainerExpression); }, py::arg("container_expression"));
    spatial_method_module.def("Sum", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Sum(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("norm"));

    spatial_method_module.def("Mean", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Mean(*pExpression, rDataCommunicator); }, py::arg("expression"), py::arg("data_communicator"));
    spatial_method_module.def("Mean", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Mean(*pExpression, rDataCommunicator, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("norm"));
    spatial_method_module.def("Mean", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Mean(*pContainerExpression); }, rContainerExpression); }, py::arg("container_expression"));
    spatial_method_module.def("Mean", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Mean(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("norm"));

    spatial_method_module.def("RootMeanSquare", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::RootMeanSquare(*pExpression, rDataCommunicator); }, py::arg("expression"), py::arg("data_communicator"));
    spatial_method_module.def("RootMeanSquare", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::RootMeanSquare(*pExpression, rDataCommunicator, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("norm"));
    spatial_method_module.def("RootMeanSquare", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::RootMeanSquare(*pContainerExpression); }, rContainerExpression); }, py::arg("container_expression"));
    spatial_method_module.def("RootMeanSquare", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::RootMeanSquare(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("norm"));

    spatial_method_module.def("Variance", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Variance(*pExpression, rDataCommunicator); }, py::arg("expression"), py::arg("data_communicator"));
    spatial_method_module.def("Variance", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Variance(*pExpression, rDataCommunicator, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("norm"));
    spatial_method_module.def("Variance", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Variance(*pContainerExpression); }, rContainerExpression); }, py::arg("container_expression"));
    spatial_method_module.def("Variance", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Variance(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("norm"));

    spatial_method_module.def("Min", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Min(*pExpression, rDataCommunicator); }, py::arg("expression"), py::arg("data_communicator"));
    spatial_method_module.def("Min", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Min(*pExpression, rDataCommunicator, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("norm"));
    spatial_method_module.def("Min", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Min(*pContainerExpression); }, rContainerExpression); }, py::arg("container_expression"));
    spatial_method_module.def("Min", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Min(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("norm"));

    spatial_method_module.def("Max", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Max(*pExpression, rDataCommunicator); }, py::arg("expression"), py::arg("data_communicator"));
    spatial_method_module.def("Max", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Max(*pExpression, rDataCommunicator, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("norm"));
    spatial_method_module.def("Max", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Max(*pContainerExpression); }, rContainerExpression); }, py::arg("container_expression"));
    spatial_method_module.def("Max", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Max(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("norm"));

    spatial_method_module.def("Median", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator) { return SpatialMethods::Median(*pExpression, rDataCommunicator); }, py::arg("expression"), py::arg("data_communicator"));
    spatial_method_module.def("Median", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Median(*pExpression, rDataCommunicator, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("norm"));
    spatial_method_module.def("Median", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression) { return std::visit([](auto& pContainerExpression) { return SpatialMethods::Median(*pContainerExpression); }, rContainerExpression); }, py::arg("container_expression"));
    spatial_method_module.def("Median", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Median(*pContainerExpression, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("norm"));

    spatial_method_module.def("Distribution", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, Parameters Params) { return SpatialMethods::Distribution(*pExpression, rDataCommunicator, Params); }, py::arg("expression"), py::arg("data_communicator"), py::arg("parameters"));
    spatial_method_module.def("Distribution", [](const Expression::ConstPointer pExpression, const DataCommunicator& rDataCommunicator, Parameters Params, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pNorm) {return SpatialMethods::Distribution(*pExpression, rDataCommunicator, Params, *pNorm); }, rNorm); }, py::arg("expression"), py::arg("data_communicator"), py::arg("parameters"), py::arg("norm"));
    spatial_method_module.def("Distribution", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, Parameters Params) { return std::visit([&](auto& pContainerExpression) { return SpatialMethods::Distribution(*pContainerExpression, Params); }, rContainerExpression); }, py::arg("container_expression"), py::arg("parameters"));
    spatial_method_module.def("Distribution", [](const typename VariantPointer<SpatialMethods::ContainerExpressionType>::type& rContainerExpression, Parameters Params, const typename VariantPointer<Norms::AllNormTypes>::type& rNorm) { return std::visit([&](auto& pContainerExpression, auto& pNorm) { return SpatialMethods::Distribution(*pContainerExpression, Params, *pNorm); }, rContainerExpression, rNorm); }, py::arg("container_expression"), py::arg("parameters"), py::arg("norm"));

}

} // namespace Python.
} // Namespace Kratos
