//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cstdint>

// External includes
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "python/numpy_utils.h"
#include "utilities/container_io_utils.h"
#include "utilities/parallel_utilities.h"

// Tensor adaptors
#include "tensor_adaptors/combined_tensor_adaptor.h"
#include "tensor_adaptors/connectivity_ids_tensor_adaptor.h"
#include "tensor_adaptors/equation_ids_tensor_adaptor.h"
#include "tensor_adaptors/flags_tensor_adaptor.h"
#include "tensor_adaptors/gauss_point_variable_tensor_adaptor.h"
#include "tensor_adaptors/historical_variable_tensor_adaptor.h"
#include "tensor_adaptors/node_position_tensor_adaptor.h"
#include "tensor_adaptors/geometries_tensor_adaptor.h"
#include "tensor_adaptors/tensor_adaptor.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"
#include "tensor_adaptors/nodal_neighbour_count_tensor_adaptor.h"
#include "tensor_adaptors/geometry_metrics_tensor_adaptor.h"

// Include base h
#include "add_tensor_adaptors_to_python.h"

namespace Kratos::Python {

namespace Detail {

template <class TDataType>
void AddBaseTensorAdaptor(
     pybind11::module& rModule,
     const std::string& rName)
{
    // add the base tensor adaptor
    using tensor_adaptor = TensorAdaptor<TDataType>;
    pybind11::class_<tensor_adaptor, typename tensor_adaptor::Pointer>(
        rModule, (rName + "Adaptor").c_str())
        .def(pybind11::init<typename tensor_adaptor::ContainerPointerType, typename NDData<TDataType>::Pointer, const bool>(),
             pybind11::arg("container"),
             pybind11::arg("nd_data"),
             pybind11::arg("copy") = true)
        .def(pybind11::init<const tensor_adaptor&, const bool>(),
             pybind11::arg("tensor_adaptor"),
             pybind11::arg("copy") = true)
        .def("Check", &tensor_adaptor::Check)
        .def("CollectData", &tensor_adaptor::CollectData)
        .def("StoreData", &tensor_adaptor::StoreData)
        .def("GetContainer", &tensor_adaptor::GetContainer)
        .def("HasContainer", &tensor_adaptor::HasContainer)
        .def("Shape", &tensor_adaptor::Shape)
        .def("DataShape", &tensor_adaptor::DataShape)
        .def("Size", &tensor_adaptor::Size)
        .def("__str__", PrintObject<tensor_adaptor>)
        .def("ViewData", [](tensor_adaptor& rSelf) { return GetPybindArray<TDataType>(*rSelf.pGetStorage()); })
        .def("SetData", [](tensor_adaptor& rSelf, const pybind11::array& rArray) { SetPybindArray<TDataType>(*rSelf.pGetStorage(), rArray); },
            pybind11::arg("array").noconvert())
        .def_property("data",
            [](tensor_adaptor& rSelf) {
                return GetPybindArray<TDataType>(*rSelf.pGetStorage());
            },
            pybind11::cpp_function(
                [](tensor_adaptor& rSelf, const pybind11::array& rArray) {
                    SetPybindArray<TDataType>(*rSelf.pGetStorage(), rArray);
                },
                pybind11::arg("self"),
                // no convert makes sure that the numpy arrays are
                // not converted, hence nothing will be copied. numpy
                // array will be passed as it is to the SetPybindArray
                // method.
                pybind11::arg("array").noconvert()));
}

template <class TDataType>
void AddCombinedTensorAdaptor(
     pybind11::module &rModule,
     const std::string &rName)
{
     using combined_ta_type = CombinedTensorAdaptor<TDataType>;
     pybind11::class_<combined_ta_type, typename combined_ta_type::Pointer, typename combined_ta_type::BaseType>(rModule, rName.c_str())
          .def(pybind11::init<const typename combined_ta_type::TensorAdaptorVectorType &, const bool, const bool>(),
               pybind11::arg("list_of_tensor_adaptors"),
               pybind11::arg("perform_collect_data_recursively") = true,
               pybind11::arg("perform_store_data_recursively") = true) // reveling ctor
          .def(pybind11::init<const typename combined_ta_type::TensorAdaptorVectorType &, const unsigned int, const bool, const bool>(),
               pybind11::arg("list_of_tensor_adaptors"), pybind11::arg("axis"),
               pybind11::arg("perform_collect_data_recursively") = true,
               pybind11::arg("perform_store_data_recursively") = true) // axis based ctor
          .def(pybind11::init<const combined_ta_type &, const bool, const bool, const bool>(),
               pybind11::arg("list_of_tensor_adaptors"),
               pybind11::arg("perform_collect_data_recursively") = true,
               pybind11::arg("perform_store_data_recursively") = true,
               pybind11::arg("copy") = true)
          .def("GetTensorAdaptors", &combined_ta_type::GetTensorAdaptors);
}

} // namespace Detail

void AddTensorAdaptorsToPython(pybind11::module& m)
{
     namespace py = pybind11;

     auto tensor_adaptor_sub_module = m.def_submodule("TensorAdaptors");
     Detail::AddBaseTensorAdaptor<bool>(tensor_adaptor_sub_module, "BoolTensor");
     Detail::AddBaseTensorAdaptor<int>(tensor_adaptor_sub_module, "IntTensor");
     Detail::AddBaseTensorAdaptor<double>(tensor_adaptor_sub_module, "DoubleTensor");

     Detail::AddCombinedTensorAdaptor<bool>(tensor_adaptor_sub_module, "BoolCombinedTensorAdaptor");
     Detail::AddCombinedTensorAdaptor<int>(tensor_adaptor_sub_module, "IntCombinedTensorAdaptor");
     Detail::AddCombinedTensorAdaptor<double>(tensor_adaptor_sub_module, "DoubleCombinedTensorAdaptor");

     py::class_<HistoricalVariableTensorAdaptor, HistoricalVariableTensorAdaptor::Pointer, HistoricalVariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "HistoricalVariableTensorAdaptor")
          .def(py::init<ModelPart::NodesContainerType::Pointer, HistoricalVariableTensorAdaptor::VariablePointerType, const int>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("step_index") = 0)
          .def(py::init<ModelPart::NodesContainerType::Pointer, TensorAdaptorUtils::VariablePointerType, const std::vector<unsigned int>&, const int>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("data_shape"),
               py::arg("step_index") = 0)
          .def(py::init<const HistoricalVariableTensorAdaptor::BaseType&, HistoricalVariableTensorAdaptor::VariablePointerType, const int, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("variable"),
               py::arg("step_index") = 0,
               py::arg("copy") = true);

     py::class_<VariableTensorAdaptor, VariableTensorAdaptor::Pointer, VariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "VariableTensorAdaptor")
          .def(py::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(),
               py::arg("container"),
               py::arg("variable"))
          .def(py::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("data_shape"))
          .def(py::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(),
               py::arg("container"),
               py::arg("variable"))
          .def(py::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("data_shape"))
          .def(py::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(),
               py::arg("container"),
               py::arg("variable"))
          .def(py::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("data_shape"))
          .def(py::init<ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(),
               py::arg("container"),
               py::arg("variable"))
          .def(py::init<ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("data_shape"))
          .def(py::init<ModelPart::GeometryContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(),
               py::arg("container"),
               py::arg("variable"))
          .def(py::init<ModelPart::GeometryContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("data_shape"))
          .def(py::init<ModelPart::MasterSlaveConstraintContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(),
               py::arg("container"),
               py::arg("variable"))
          .def(py::init<ModelPart::MasterSlaveConstraintContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("data_shape"))
          .def(py::init<const VariableTensorAdaptor::BaseType&, VariableTensorAdaptor::VariablePointerType, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("variable"),
               py::arg("copy") = true);

     py::class_<GaussPointVariableTensorAdaptor, GaussPointVariableTensorAdaptor::Pointer, GaussPointVariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "GaussPointVariableTensorAdaptor")
          .def(py::init<ModelPart::ConditionsContainerType::Pointer, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("process_info"))
          .def(py::init<ModelPart::ElementsContainerType::Pointer, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer>(),
               py::arg("container"),
               py::arg("variable"),
               py::arg("process_info"))
          .def(py::init<const GaussPointVariableTensorAdaptor::BaseType&, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("variable"),
               py::arg("process_info"),
               py::arg("copy") = true);

     py::class_<EquationIdsTensorAdaptor, EquationIdsTensorAdaptor::Pointer, EquationIdsTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "EquationIdsTensorAdaptor")
          .def(py::init<ModelPart::ConditionsContainerType::Pointer, ProcessInfo::Pointer>(),
               py::arg("container"),
               py::arg("process_info"))
          .def(py::init<ModelPart::ElementsContainerType::Pointer, ProcessInfo::Pointer>(),
               py::arg("container"),
               py::arg("process_info"))
          .def(py::init<const EquationIdsTensorAdaptor::BaseType&, ProcessInfo::Pointer, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("process_info"),
               py::arg("copy") = true);

     py::class_<FlagsTensorAdaptor, FlagsTensorAdaptor::Pointer, FlagsTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "FlagsTensorAdaptor")
          .def(py::init<ModelPart::NodesContainerType::Pointer, const Flags&>(),
               py::arg("container"),
               py::arg("flag"))
          .def(py::init<ModelPart::ConditionsContainerType::Pointer, const Flags&>(),
               py::arg("container"),
               py::arg("flag"))
          .def(py::init<ModelPart::ElementsContainerType::Pointer, const Flags&>(),
               py::arg("container"),
               py::arg("flag"))
          .def(py::init<const FlagsTensorAdaptor::BaseType&, const Flags&, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("flag"), py::arg("copy") = true);

     py::class_<NodePositionTensorAdaptor, NodePositionTensorAdaptor::Pointer, NodePositionTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "NodePositionTensorAdaptor")
          .def(py::init<ModelPart::NodesContainerType::Pointer, Globals::Configuration>(),
               py::arg("container"),
               py::arg("configuration"))
          .def(py::init<ModelPart::NodesContainerType::Pointer, Globals::Configuration, const std::vector<unsigned int>&>(),
               py::arg("container"),
               py::arg("configuration"),
               py::arg("data_shape"))
          .def(py::init<const NodePositionTensorAdaptor::BaseType&, Globals::Configuration, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("configuration"),
               py::arg("copy") = true);

     using GeometriesTensorAdaptorType = GeometriesTensorAdaptor;
     py::class_<GeometriesTensorAdaptorType, GeometriesTensorAdaptorType::Pointer, GeometriesTensorAdaptorType::BaseType> geometries_tensor_adaptor(tensor_adaptor_sub_module, "GeometriesTensorAdaptor");

     py::enum_<GeometriesTensorAdaptorType::DatumType>(geometries_tensor_adaptor,"DatumType")
          .value("ShapeFunctions", GeometriesTensorAdaptorType::DatumType::ShapeFunctions)
          .value("ShapeFunctionDerivatives", GeometriesTensorAdaptorType::DatumType::ShapeFunctionDerivatives)
          .value("Jacobians", GeometriesTensorAdaptorType::DatumType::Jacobians)
          .value("IntegrationWeights", GeometriesTensorAdaptorType::DatumType::IntegrationWeights)
          .export_values();

     geometries_tensor_adaptor
          .def(py::init<ModelPart::GeometryContainerType::Pointer, GeometriesTensorAdaptorType::DatumType,GeometryData::IntegrationMethod>(),
               py::arg("container"),
               py::arg("datum"),
               py::arg("integration_method") = GeometryData::IntegrationMethod::NumberOfIntegrationMethods)
          .def(py::init<ModelPart::ElementsContainerType::Pointer, GeometriesTensorAdaptorType::DatumType, GeometryData::IntegrationMethod>(),
               py::arg("container"),
               py::arg("datum"),
               py::arg("integration_method") = GeometryData::IntegrationMethod::NumberOfIntegrationMethods)
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, GeometriesTensorAdaptorType::DatumType, GeometryData::IntegrationMethod>(),
               py::arg("container"),
               py::arg("datum"),
               py::arg("integration_method") = GeometryData::IntegrationMethod::NumberOfIntegrationMethods)
          .def(py::init<const GeometriesTensorAdaptorType::BaseType&, ModelPart::GeometryContainerType::Pointer, GeometriesTensorAdaptorType::DatumType, GeometryData::IntegrationMethod, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("container"),
               py::arg("datum"),
               py::arg("integration_method"),
               py::arg("copy") = true)
          .def(py::init<const GeometriesTensorAdaptorType::BaseType&, ModelPart::ElementsContainerType::Pointer, GeometriesTensorAdaptorType::DatumType, GeometryData::IntegrationMethod, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("container"),
               py::arg("datum"),
               py::arg("integration_method"),
               py::arg("copy") = true)
          .def(py::init<const GeometriesTensorAdaptorType::BaseType&, ModelPart::ConditionsContainerType::Pointer, GeometriesTensorAdaptorType::DatumType, GeometryData::IntegrationMethod, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("container"),
               py::arg("datum"),
               py::arg("integration_method"),
               py::arg("copy") = true);

     py::class_<ConnectivityIdsTensorAdaptor, ConnectivityIdsTensorAdaptor::Pointer, ConnectivityIdsTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "ConnectivityIdsTensorAdaptor")
          .def(py::init<ModelPart::GeometryContainerType::Pointer>(),
               py::arg("container"))
          .def(py::init<ModelPart::ElementsContainerType::Pointer>(),
               py::arg("container"))
          .def(py::init<ModelPart::ConditionsContainerType::Pointer>(),
               py::arg("container"))
          .def(py::init<const ConnectivityIdsTensorAdaptor::BaseType&, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("copy") = true);

     py::class_<NodalNeighbourCountTensorAdaptor, NodalNeighbourCountTensorAdaptor::Pointer, NodalNeighbourCountTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "NodalNeighbourCountTensorAdaptor")
          .def(py::init<const NodalNeighbourCountTensorAdaptor::BaseType&, ModelPart::ConditionsContainerType::Pointer, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("entity_container"),
               py::arg("copy") = true)
          .def(py::init<const NodalNeighbourCountTensorAdaptor::BaseType&, ModelPart::ElementsContainerType::Pointer, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("entity_container"),
               py::arg("copy") = true)
          .def(py::init<ModelPart::NodesContainerType::Pointer, ModelPart::ConditionsContainerType::Pointer>(),
               py::arg("nodes_container"),
               py::arg("entity_container"))
          .def(py::init<ModelPart::NodesContainerType::Pointer, ModelPart::ElementsContainerType::Pointer>(),
               py::arg("nodes_container"),
               py::arg("entity_container"));

     py::class_<GeometryMetricsTensorAdaptor, GeometryMetricsTensorAdaptor::Pointer, GeometryMetricsTensorAdaptor::BaseType> geometric_tensor_adaptor(tensor_adaptor_sub_module, "GeometryMetricsTensorAdaptor");

     py::enum_<GeometryMetricsTensorAdaptor::DatumType>(geometric_tensor_adaptor,"DatumType")
          .value("DomainSize", GeometryMetricsTensorAdaptor::DatumType::DomainSize)
          .export_values();

     geometric_tensor_adaptor
          .def(py::init<const GeometryMetricsTensorAdaptor::BaseType&, GeometryMetricsTensorAdaptor::DatumType, const bool>(),
               py::arg("tensor_adaptor"),
               py::arg("datum"),
               py::arg("copy") = true)
          .def(py::init<ModelPart::GeometryContainerType::Pointer, GeometryMetricsTensorAdaptor::DatumType>(),
               py::arg("container"),
               py::arg("datum"))
          .def(py::init<ModelPart::ConditionsContainerType::Pointer, GeometryMetricsTensorAdaptor::DatumType>(),
               py::arg("container"),
               py::arg("datum"))
          .def(py::init<ModelPart::ElementsContainerType::Pointer, GeometryMetricsTensorAdaptor::DatumType>(),
               py::arg("container"),
               py::arg("datum"));
}

} // namespace Kratos::Python.