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
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

// Project includes
#include "includes/model_part.h"
#include "includes/define_python.h"
#include "utilities/container_io_utils.h"
#include "tensor_adaptors/tensor_adaptor.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"
#include "numpy_utils.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "add_tensor_adaptors_to_python.h"

namespace Kratos::Python {

namespace Detail {

using PybindArrayType = std::variant<
                            pybind11::array_t<int, pybind11::array::c_style>,
                            pybind11::array_t<bool, pybind11::array::c_style>,
                            pybind11::array_t<double, pybind11::array::c_style>
                        >;


template<class TDataType>
class TensorAdaptorTrampoline final : public TensorAdaptor<TDataType> {
public:
    using BaseType = TensorAdaptor<TDataType>;

    using ContainerType = typename BaseType::ContainerType;

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

template<class TDataType>
PybindArrayType GetNumpyArray(TensorAdaptor<TDataType>& rTensorAdaptor)
{
    using numpy_array_type = pybind11::array_t<TDataType, pybind11::array::c_style>;

    const auto& r_shape = rTensorAdaptor.Shape();

    std::vector<std::size_t> c_shape(r_shape.size());
    std::copy(r_shape.begin(), r_shape.end(), c_shape.begin());
    std::vector<std::size_t> strides(c_shape.size());

    std::size_t stride_items = 1;
    for (int i = c_shape.size() - 1; i >= 0; --i) {
        strides[i] = sizeof(TDataType) * stride_items;
        stride_items *= c_shape[i];
    }

    if (rTensorAdaptor.Size() > 0) {
        // do nothing in the release of the numpy array since the ownership is not passed
        // the ownership is kept with the TensorAdaptor.
        pybind11::capsule release(rTensorAdaptor.ViewData().data(), [](void* a){});
        return numpy_array_type(pybind11::buffer_info(
            rTensorAdaptor.ViewData().data(),                   // Pointer to data
            sizeof(TDataType),                                  // Size of one item
            pybind11::format_descriptor<TDataType>::format(),   // Python format descriptor
            c_shape.size(),                                     // Number of dimensions
            c_shape,                                            // Shape of the array
            strides                                             // Strides
        ), release);
    } else {
        return numpy_array_type(c_shape, strides);
    }
}

template<class TDataType, class... TArgs>
void SetNumpyArray(
    TensorAdaptor<TDataType>& rTensorAdaptor,
    std::variant<pybind11::array_t<TArgs, pybind11::array::c_style> const *...> pArray)
{
    std::visit([&rTensorAdaptor](const auto pArray) {
        auto& r_array = *pArray;

        KRATOS_ERROR_IF(r_array.ndim() == 0)
            << "Passed data is not compatible [ array = "
            << r_array << ", Tensor adaptor = " << rTensorAdaptor << " ].\n";

        std::vector<int> shape(r_array.ndim());
        std::copy(r_array.shape(), r_array.shape() + r_array.ndim(), shape.begin());

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
        IndexPartition<IndexType>(rTensorAdaptor.ViewData().size()).for_each([&r_array, &rTensorAdaptor](const auto Index) {
            rTensorAdaptor.ViewData()[Index] = static_cast<TDataType>(r_array.data()[Index]);
        });
    }, pArray);
}

template<class TDataType>
void AddBaseTensorAdaptor(
    pybind11::module& rModule,
    const std::string& rName)
{
    // add the base tensor adaptor
    using tensor_adaptor = TensorAdaptor<TDataType>;
    pybind11::class_<tensor_adaptor, typename tensor_adaptor::Pointer>(rModule, rName.c_str())
        .def("CollectData", &tensor_adaptor::CollectData)
        .def("StoreData", &tensor_adaptor::StoreData)
        .def("GetContainer", &tensor_adaptor::GetContainer)
        .def("Shape", &tensor_adaptor::Shape)
        .def("GetDataShape", &tensor_adaptor::GetDataShape)
        .def("Size", &tensor_adaptor::Size)
        .def("__str__", PrintObject<tensor_adaptor>)
        .def("ViewData", &Detail::GetNumpyArray<TDataType>)
        .def_property("data",
                      &Detail::GetNumpyArray<TDataType>,
                      pybind11::cpp_function(
                        &Detail::SetNumpyArray<
                            TDataType,
                            bool,
                            std::uint8_t,
                            std::uint16_t,
                            std::uint32_t,
                            std::uint64_t,
                            std::int8_t,
                            std::int16_t,
                            std::int32_t,
                            std::int64_t,
                            float,
                            double,
                            long double>,
                        pybind11::arg("self"),
                        pybind11::arg("array").noconvert()
                      ))
    ;
}

} // namespace Detail

void AddTensorAdaptorsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto tensor_adaptor_sub_module = m.def_submodule("TensorAdaptors");
    Detail::AddBaseTensorAdaptor<bool>(tensor_adaptor_sub_module, "BoolTensorAdaptor");
    Detail::AddBaseTensorAdaptor<int>(tensor_adaptor_sub_module, "IntTensorAdaptor");
    Detail::AddBaseTensorAdaptor<double>(tensor_adaptor_sub_module, "DoubleTensorAdaptor");

    using historical_variable_tensor_adaptor = VariableTensorAdaptor<HistoricalIO, const int>;
    py::class_<historical_variable_tensor_adaptor, historical_variable_tensor_adaptor::Pointer, historical_variable_tensor_adaptor::BaseType>(tensor_adaptor_sub_module, "HistoricalVariableTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, historical_variable_tensor_adaptor::VariableType, const int>(), py::arg("container"), py::arg("variable"), py::arg("step_index") = 0)
        .def(py::init<ModelPart::NodesContainerType::Pointer, historical_variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&, const int>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"), py::arg("step_index") = 0)
        ;

    using gauss_point_variable_tensor_adaptor = VariableTensorAdaptor<GaussPointIO, const ProcessInfo&>;
    py::class_<gauss_point_variable_tensor_adaptor, gauss_point_variable_tensor_adaptor::Pointer, gauss_point_variable_tensor_adaptor::BaseType>(tensor_adaptor_sub_module, "GaussPointVariableTensorAdaptor")
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, gauss_point_variable_tensor_adaptor::VariableType, const ProcessInfo&>(), py::arg("container"), py::arg("variable"), py::arg("process_info"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, gauss_point_variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&, const ProcessInfo&>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"), py::arg("process_info"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, gauss_point_variable_tensor_adaptor::VariableType, const ProcessInfo&>(), py::arg("container"), py::arg("variable"), py::arg("process_info"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, gauss_point_variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&, const ProcessInfo&>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"), py::arg("process_info"))
        ;

    using variable_tensor_adaptor = VariableTensorAdaptor<NonHistoricalIO>;
    pybind11::class_<variable_tensor_adaptor, variable_tensor_adaptor::Pointer, variable_tensor_adaptor::BaseType>(tensor_adaptor_sub_module, "VariableTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, variable_tensor_adaptor::VariableType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::NodesContainerType::Pointer, variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, variable_tensor_adaptor::VariableType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, variable_tensor_adaptor::VariableType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"))
        .def(py::init<ModelPart::PropertiesContainerType::Pointer, variable_tensor_adaptor::VariableType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::PropertiesContainerType::Pointer, variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"))
        .def(py::init<ModelPart::GeometriesMapType::Pointer, variable_tensor_adaptor::VariableType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::GeometriesMapType::Pointer, variable_tensor_adaptor::VariableType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("return_data_shape"))
        ;

}

} // namespace Kratos::Python.