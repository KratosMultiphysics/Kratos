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
#include "tensor_adaptors/tensor_adaptor.h"
#include "includes/model_part.h"
#include "includes/define_python.h"

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

} // namespace Detail


void AddTensorAdaptorsToPython(pybind11::module& m)
{
    auto tensor_adaptor_sub_module = m.def_submodule("TensorAdaptors");

    auto node_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("NodeTensorAdaptors");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::NodesContainerType, int>>(node_tensor_adaptors_sub_module, "NodeIntTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::NodesContainerType, double>>(node_tensor_adaptors_sub_module, "NodeDoubleTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::NodesContainerType, bool>>(node_tensor_adaptors_sub_module, "NodeBoolTensorAdaptor");

    auto condition_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("ConditionTensorAdaptors");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::ConditionsContainerType, int>>(condition_tensor_adaptors_sub_module, "ConditionIntTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::ConditionsContainerType, double>>(condition_tensor_adaptors_sub_module, "ConditionDoubleTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::ConditionsContainerType, bool>>(condition_tensor_adaptors_sub_module, "ConditionBoolTensorAdaptor");

    auto element_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("ElementTensorAdaptors");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::ElementsContainerType, int>>(element_tensor_adaptors_sub_module, "ElementIntTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::ElementsContainerType, double>>(element_tensor_adaptors_sub_module, "ElementDoubleTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::ElementsContainerType, bool>>(element_tensor_adaptors_sub_module, "ElementBoolTensorAdaptor");

    auto properties_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("PropertyTensorAdaptors");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::PropertiesContainerType, int>>(properties_tensor_adaptors_sub_module, "PropertyIntTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::PropertiesContainerType, double>>(properties_tensor_adaptors_sub_module, "PropertyDoubleTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::PropertiesContainerType, bool>>(properties_tensor_adaptors_sub_module, "PropertyBoolTensorAdaptor");

    auto geometries_tensor_adaptors_sub_module = tensor_adaptor_sub_module.def_submodule("GeometryTensorAdaptors");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::GeometryContainerType, int>>(geometries_tensor_adaptors_sub_module, "GeometryIntTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::GeometryContainerType, double>>(geometries_tensor_adaptors_sub_module, "GeometryDoubleTensorAdaptor");
    Detail::AddBaseTensorAdaptor<TensorAdaptor<ModelPart::GeometryContainerType, bool>>(geometries_tensor_adaptors_sub_module, "GeometryBoolTensorAdaptor");

}

} // namespace Kratos::Python.