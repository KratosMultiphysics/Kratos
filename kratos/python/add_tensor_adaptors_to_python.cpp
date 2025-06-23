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

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define_python.h"
#include "utilities/container_io_utils.h"
#include "tensor_adaptors/tensor_adaptor.h"
#include "tensor_adaptors/variant_variable_tensor_adaptor.h"

// Include base h
#include "add_tensor_adaptors_to_python.h"

namespace Kratos::Python {

namespace Detail {

template <class TContainerType, class TDataType>
class TensorAdaptorTrampoline final : public TensorAdaptor<TContainerType, TDataType> {
public:
    using BaseType = TensorAdaptor<TContainerType, TDataType>;
    void CollectData() override
    {
        PYBIND11_OVERRIDE_PURE(void,        /*return type*/
                               BaseType,    /*base type*/
                               CollectData  /*function name*/
        );
    }

    void StoreData() override
    {
        PYBIND11_OVERRIDE_PURE(void,        /*return type*/
                               BaseType,    /*base type*/
                               StoreData    /*function name*/
        );
    }
}; // class ExpressionTrampoline

template<class TTensorAdapterType>
void AddBaseTensorAdaptor(
    pybind11::module& rModule,
    const std::string& AdaptorName)
{
    using container_type = typename TTensorAdapterType::ContainerType;
    using primitive_data_type = typename TTensorAdapterType::PrimitiveDataType;
    pybind11::class_<TTensorAdapterType, typename TTensorAdapterType::Pointer, Detail::TensorAdaptorTrampoline<container_type, primitive_data_type>>(rModule, AdaptorName.c_str())
        .def("CollectData", &TTensorAdapterType::CollectData)
        .def("StoreData", &TTensorAdapterType::StoreData)
        .def("GetContainer", &TTensorAdapterType::GetContainer)
        .def("Shape", &TTensorAdapterType::Shape)
        .def("ViewData", pybind11::overload_cast<>(&TTensorAdapterType::ViewData))
        .def("MoveData", &TTensorAdapterType::MoveData)
        .def("__str__", PrintObject<TTensorAdapterType>);
    ;
}

template<class TContainerType>
void AddContainerBaseTensorAdaptors(
    pybind11::module& rModule,
    const std::string& AdaptorNamePrefix)
{
    AddBaseTensorAdaptor<TensorAdaptor<TContainerType, int>>(rModule, AdaptorNamePrefix + "IntTensorAdaptor");
    AddBaseTensorAdaptor<TensorAdaptor<TContainerType, double>>(rModule, AdaptorNamePrefix + "DoubleTensorAdaptor");
    AddBaseTensorAdaptor<TensorAdaptor<TContainerType, bool>>(rModule, AdaptorNamePrefix + "BoolTensorAdaptor");
}

template<class TContainerType>
void AddNonHistoricalTensorAdaptors(
    pybind11::module& rModule,
    const std::string& rAdaptorName)
{
    using non_historical_variant_variable_ta_type = VariantVariableTensorAdaptor<TContainerType, NonHistoricalIO>;
    pybind11::class_<non_historical_variant_variable_ta_type, typename non_historical_variant_variable_ta_type::Pointer, typename non_historical_variant_variable_ta_type::BaseType>(rModule, rAdaptorName.c_str())
        .def(pybind11::init<typename non_historical_variant_variable_ta_type::ContainerType::Pointer, typename non_historical_variant_variable_ta_type::VariableType>(), pybind11::arg("container"), pybind11::arg("variable"));

}

} // namespace Detail


void AddTensorAdaptorsToPython(pybind11::module& m)
{
    auto tensor_adaptor_sub_module = m.def_submodule("TensorAdaptors");

    auto node_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("NodeTensorAdaptors");
    Detail::AddContainerBaseTensorAdaptors<ModelPart::NodesContainerType>(node_tensor_adaptors_sub_module, "Node");

    auto condition_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("ConditionTensorAdaptors");
    Detail::AddContainerBaseTensorAdaptors<ModelPart::ConditionsContainerType>(condition_tensor_adaptors_sub_module, "Condition");

    auto element_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("ElementTensorAdaptors");
    Detail::AddContainerBaseTensorAdaptors<ModelPart::ElementsContainerType>(element_tensor_adaptors_sub_module, "Element");

    auto properties_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("PropertyTensorAdaptors");
    Detail::AddContainerBaseTensorAdaptors<ModelPart::PropertiesContainerType>(properties_tensor_adaptors_sub_module, "Property");

    auto geometries_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("GeometryTensorAdaptors");
    Detail::AddContainerBaseTensorAdaptors<ModelPart::GeometryContainerType::GeometriesMapType>(geometries_tensor_adaptors_sub_module, "Geometry");

    // adding non historical variable tensor adaptors
    Detail::AddNonHistoricalTensorAdaptors<ModelPart::NodesContainerType>(node_tensor_adaptors_sub_module, "NodeNonHistoricalVariableTensorAdaptor");
    Detail::AddNonHistoricalTensorAdaptors<ModelPart::ConditionsContainerType>(condition_tensor_adaptors_sub_module, "ConditionVariableTensorAdaptor");
    Detail::AddNonHistoricalTensorAdaptors<ModelPart::ElementsContainerType>(element_tensor_adaptors_sub_module, "ElementVariableTensorAdaptor");
    Detail::AddNonHistoricalTensorAdaptors<ModelPart::PropertiesContainerType>(properties_tensor_adaptors_sub_module, "PropertyVariableTensorAdaptor");
    Detail::AddNonHistoricalTensorAdaptors<ModelPart::GeometryContainerType::GeometriesMapType>(geometries_tensor_adaptors_sub_module, "GeometryVariableTensorAdaptor");

    // adding non-historical variable tensor adaptors
    using node_historical_variant_variable_ta_type = VariantVariableTensorAdaptor<ModelPart::NodesContainerType, HistoricalIO, const int>;
    pybind11::class_<node_historical_variant_variable_ta_type, node_historical_variant_variable_ta_type::Pointer, node_historical_variant_variable_ta_type::BaseType>(node_tensor_adaptors_sub_module, "NodeHistoricalVariableTensorAdaptor")
        .def(pybind11::init<typename node_historical_variant_variable_ta_type::ContainerType::Pointer, typename node_historical_variant_variable_ta_type::VariableType, const int>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("step_index"));

    // adding gauss point calculation tensor adaptors
    using condition_gp_variant_variable_ta_type = VariantVariableTensorAdaptor<ModelPart::ConditionsContainerType, GaussPointIO, const ProcessInfo&>;
    pybind11::class_<condition_gp_variant_variable_ta_type, condition_gp_variant_variable_ta_type::Pointer, condition_gp_variant_variable_ta_type::BaseType>(condition_tensor_adaptors_sub_module, "ConditionGaussPointVariableTensorAdaptor")
        .def(pybind11::init<typename condition_gp_variant_variable_ta_type::ContainerType::Pointer, typename condition_gp_variant_variable_ta_type::VariableType, const ProcessInfo&>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("process_info"));

    using element_gp_variant_variable_ta_type = VariantVariableTensorAdaptor<ModelPart::ElementsContainerType, GaussPointIO, const ProcessInfo&>;
    pybind11::class_<element_gp_variant_variable_ta_type, element_gp_variant_variable_ta_type::Pointer, element_gp_variant_variable_ta_type::BaseType>(element_tensor_adaptors_sub_module, "ElementGaussPointVariableTensorAdaptor")
        .def(pybind11::init<typename element_gp_variant_variable_ta_type::ContainerType::Pointer, typename element_gp_variant_variable_ta_type::VariableType, const ProcessInfo&>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("process_info"));

}

} // namespace Kratos::Python.