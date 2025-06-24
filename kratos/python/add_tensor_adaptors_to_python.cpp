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
#include <pybind11/stl.h>

// Project includes
#include "includes/model_part.h"
#include "includes/define_python.h"
#include "utilities/container_io_utils.h"
#include "tensor_adaptors/tensor_adaptor.h"
#include "tensor_adaptors/variant_variable_tensor_adaptor.h"
#include "numpy_utils.h"
#include "utilities/parallel_utilities.h"

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

    std::string Info() const override
    {
        PYBIND11_OVERRIDE_PURE(std::string,  /*return type*/
                               BaseType,     /*base type*/
                               Info          /*function name*/
        );
    }
}; // class ExpressionTrampoline

template<class TTensorAdapterType>
pybind11::array_t<typename TTensorAdapterType::PrimitiveDataType, pybind11::array::c_style> GetNumpyArray(TTensorAdapterType& rTensorAdaptor)
{
    using primitive_data_type = typename TTensorAdapterType::PrimitiveDataType;
    using numpy_array_type = pybind11::array_t<primitive_data_type, pybind11::array::c_style>;

    const auto& r_shape = rTensorAdaptor.Shape();

    std::vector<std::size_t> c_shape(r_shape.size());
    std::copy(r_shape.begin(), r_shape.end(), c_shape.begin());
    std::vector<std::size_t> strides(c_shape.size());

    std::size_t stride_items = 1;
    for (int i = c_shape.size() - 1; i >= 0; --i) {
        strides[i] = sizeof(primitive_data_type) * stride_items;
        stride_items *= c_shape[i];
    }

    // do nothing in the release of the numpy array since the ownership is not passed
    // the ownership is kept with the TensorAdaptor.
    pybind11::capsule release(rTensorAdaptor.ViewData().data(), [](void* a) {});

    return numpy_array_type(
        c_shape,
        strides,
        rTensorAdaptor.ViewData().data(),
        release
    );
}

template<class TTensorAdapterType>
void SetNumpyArray(
    TTensorAdapterType& rTensorAdaptor,
    const pybind11::array_t<typename TTensorAdapterType::PrimitiveDataType, pybind11::array::c_style>& rArray)
{
    KRATOS_ERROR_IF(rArray.ndim() == 0)
        << "Passed data is not compatible [ array = "
        << rArray << ", Tensor adaptor = " << rTensorAdaptor << " ].\n";

    std::vector<int> shape(rArray.ndim());
    std::copy(rArray.shape(), rArray.shape() + rArray.ndim(), shape.begin());

    const auto& r_shape = rTensorAdaptor.Shape();

    KRATOS_ERROR_IF_NOT(shape.size() == r_shape.size())
        << "Dimensions mismatch. [ Tensor dimensions = " << r_shape.size()
        << ", numpy array dimensions = " << shape.size()
        << ", Tensor adaptor = " << rTensorAdaptor << " ].\n";

    for (unsigned int i = 0; i < shape.size(); ++i) {
        KRATOS_ERROR_IF_NOT(r_shape[i] == shape[i])
            << "Shape mismatch. [ Tensor shape = " << rTensorAdaptor.Shape()
            << ", numpy array shape = " << shape
            << ", Tensor adaptor = " << rTensorAdaptor << " ].\n";
    }

    // copy data from the input to the Adaptor
    IndexPartition<IndexType>(rArray.size()).for_each([&rArray, &rTensorAdaptor](const auto Index) {
        rTensorAdaptor.ViewData()[Index] = rArray.data()[Index];
    });
}

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
        .def_property("data", &GetNumpyArray<TTensorAdapterType>, &SetNumpyArray<TTensorAdapterType>)
        .def("ViewData", &GetNumpyArray<TTensorAdapterType>)
        .def("MoveData", [](TTensorAdapterType& rSelf){
                using numpy_array_type = pybind11::array_t<primitive_data_type, pybind11::array::c_style>;

                const auto& r_shape = rSelf.Shape();

                std::vector<std::size_t> c_shape(r_shape.size());
                std::copy(r_shape.begin(), r_shape.end(), c_shape.begin());
                std::vector<std::size_t> strides(c_shape.size());

                std::size_t stride_items = 1;
                for (int i = c_shape.size() - 1; i >= 0; --i) {
                    strides[i] = sizeof(primitive_data_type) * stride_items;
                    stride_items *= c_shape[i];
                }

                // this method transfers the ownership of the data to numpy.
                // Since the DenseVector does not allow moving out the data without getting a call
                // to destroy DenseVector underlying data at the destructor, we are creating a copy of the
                // data by moving it to a temp so the underlying data container within the TensorAdaptor is cleared.
                DenseVector<primitive_data_type> temp = rSelf.MoveData();
                primitive_data_type* array = new primitive_data_type[temp.size()];
                std::copy(temp.begin(), temp.end(), array);

                // now we add the release to clear the data.
                pybind11::capsule release(array, [](void* a) {
                    delete[] reinterpret_cast<primitive_data_type*>(a);
                });

                return numpy_array_type(
                    c_shape,
                    strides,
                    array,
                    release
                );
        })
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
        .def(pybind11::init<typename non_historical_variant_variable_ta_type::ContainerType::Pointer, typename non_historical_variant_variable_ta_type::VariableType>(), pybind11::arg("container"), pybind11::arg("variable"))
        .def(pybind11::init<typename non_historical_variant_variable_ta_type::ContainerType::Pointer, typename non_historical_variant_variable_ta_type::VariableType, const std::vector<int>&>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("shape"))
        ;

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
        .def(pybind11::init<typename node_historical_variant_variable_ta_type::ContainerType::Pointer, typename node_historical_variant_variable_ta_type::VariableType, const int>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("step_index") = 0)
        .def(pybind11::init<typename node_historical_variant_variable_ta_type::ContainerType::Pointer, typename node_historical_variant_variable_ta_type::VariableType, const std::vector<int>&, const int>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("shape"), pybind11::arg("step_index") = 0)
        ;

    // adding gauss point calculation tensor adaptors
    using condition_gp_variant_variable_ta_type = VariantVariableTensorAdaptor<ModelPart::ConditionsContainerType, GaussPointIO, const ProcessInfo&>;
    pybind11::class_<condition_gp_variant_variable_ta_type, condition_gp_variant_variable_ta_type::Pointer, condition_gp_variant_variable_ta_type::BaseType>(condition_tensor_adaptors_sub_module, "ConditionGaussPointVariableTensorAdaptor")
        .def(pybind11::init<typename condition_gp_variant_variable_ta_type::ContainerType::Pointer, typename condition_gp_variant_variable_ta_type::VariableType, const ProcessInfo&>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("process_info"))
        .def(pybind11::init<typename condition_gp_variant_variable_ta_type::ContainerType::Pointer, typename condition_gp_variant_variable_ta_type::VariableType, const std::vector<int>&, const ProcessInfo&>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("shape"), pybind11::arg("process_info"))
        ;

    using element_gp_variant_variable_ta_type = VariantVariableTensorAdaptor<ModelPart::ElementsContainerType, GaussPointIO, const ProcessInfo&>;
    pybind11::class_<element_gp_variant_variable_ta_type, element_gp_variant_variable_ta_type::Pointer, element_gp_variant_variable_ta_type::BaseType>(element_tensor_adaptors_sub_module, "ElementGaussPointVariableTensorAdaptor")
        .def(pybind11::init<typename element_gp_variant_variable_ta_type::ContainerType::Pointer, typename element_gp_variant_variable_ta_type::VariableType, const ProcessInfo&>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("process_info"))
        .def(pybind11::init<typename element_gp_variant_variable_ta_type::ContainerType::Pointer, typename element_gp_variant_variable_ta_type::VariableType, const std::vector<int>&, const ProcessInfo&>(), pybind11::arg("container"), pybind11::arg("variable"), pybind11::arg("shape"), pybind11::arg("process_info"))
        ;

}

} // namespace Kratos::Python.