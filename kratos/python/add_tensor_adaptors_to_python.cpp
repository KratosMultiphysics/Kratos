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
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

// Project includes
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "utilities/container_io_utils.h"
#include "utilities/parallel_utilities.h"

// Tensor adaptors
#include "tensor_adaptors/tensor_adaptor.h"
#include "tensor_adaptors/combined_tensor_adaptor.h"
#include "tensor_adaptors/historical_variable_tensor_adaptor.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"
#include "tensor_adaptors/flags_tensor_adaptor.h"
#include "tensor_adaptors/equation_ids_tensor_adaptor.h"
#include "tensor_adaptors/gauss_point_variable_tensor_adaptor.h"
#include "tensor_adaptors/node_position_tensor_adaptor.h"

// Include base h
#include "add_tensor_adaptors_to_python.h"

namespace Kratos::Python {

namespace Detail {

template<template<class> class TTensorType, class TDataType>
pybind11::array_t<TDataType> GetPybindArray(TTensorType<TDataType>& rTensorAdaptor)
{
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
        return pybind11::array_t<TDataType>(pybind11::buffer_info(
            rTensorAdaptor.ViewData().data(),                   // Pointer to data
            sizeof(TDataType),                                  // Size of one item
            pybind11::format_descriptor<TDataType>::format(),   // Python format descriptor
            c_shape.size(),                                     // Number of dimensions
            c_shape,                                            // Shape of the array
            strides                                             // Strides
        ), release);
    } else {
        return pybind11::array_t<TDataType>(c_shape, strides);
    }
}

template<template<class> class TTensorType, class TDataType>
pybind11::array_t<TDataType> MovePybindArray(TensorAdaptor<TDataType>& rTensorAdaptor)
{
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
        auto data_span = rTensorAdaptor.MoveData();
        pybind11::capsule release(data_span.data(), [](void* a){ delete[] reinterpret_cast<TDataType*>(a); });
        return pybind11::array_t<TDataType>(pybind11::buffer_info(
            data_span.data(),                   // Pointer to data
            sizeof(TDataType),                                  // Size of one item
            pybind11::format_descriptor<TDataType>::format(),   // Python format descriptor
            c_shape.size(),                                     // Number of dimensions
            c_shape,                                            // Shape of the array
            strides                                             // Strides
        ), release);
    } else {
        return pybind11::array_t<TDataType>(c_shape, strides);
    }
}

template<template<class> class TTensorType, class TTensorStorageType, class TPybindArrayType>
bool AssignDataImpl(
    TTensorType<TTensorStorageType>& rTensor,
    const pybind11::array& rArray)
{
    if (pybind11::isinstance<pybind11::array_t<TPybindArrayType>>(rArray)) {

        KRATOS_ERROR_IF_NOT(rArray.flags() & pybind11::detail::npy_api::constants::NPY_ARRAY_C_CONTIGUOUS_)
            << "Only supports C-style (row-major) arrays from numpy.";

        const auto& casted_array = rArray.cast<pybind11::array_t<TPybindArrayType, pybind11::array::c_style>>();
        auto r_destination_span = rTensor.ViewData();
        const auto& r_origin_data = casted_array.data();
        IndexPartition<IndexType>(rTensor.ViewData().size()).for_each([&r_destination_span, &r_origin_data](const auto Index) {
            r_destination_span[Index] = static_cast<TTensorStorageType>(r_origin_data[Index]);
        });
        return true;
    } else {
        return false;
    }
}

template<template<class> class TTensorType, class TTensorStorageType, class... TPybindArrayType>
bool AssignData(
    TTensorType<TTensorStorageType>& rTensor,
    const pybind11::array& rArray)
{
    return (... || AssignDataImpl<TTensorType, TTensorStorageType, TPybindArrayType>(rTensor, rArray));
}

template<template<class> class TTensorType, class TTensorStorageType>
void SetPybindArray(
    TTensorType<TTensorStorageType>& rTensor,
    const pybind11::array&  rArray)
{
    KRATOS_ERROR_IF(rArray.ndim() == 0)
        << "Passed data is not compatible [ array = "
        << rArray << ", Tensor adaptor = " << rTensor << " ].\n";

    std::vector<unsigned int> shape(rArray.ndim());
    std::copy(rArray.shape(), rArray.shape() + rArray.ndim(), shape.begin());

    const auto& r_shape = rTensor.Shape();

    KRATOS_ERROR_IF_NOT(shape.size() == r_shape.size())
        << "Dimensions mismatch. [ Tensor dimensions = " << r_shape.size()
        << ", numpy array dimensions = " << shape.size()
        << ", Tensor adaptor = " << rTensor << " ].\n";

    for (unsigned int i = 0; i < shape.size(); ++i) {
        KRATOS_ERROR_IF_NOT(r_shape[i] == shape[i])
            << "Shape mismatch. [ Tensor shape = " << rTensor.Shape()
            << ", numpy array shape = " << shape
            << ", Tensor adaptor = " << rTensor << " ].\n";
    }

    if (!AssignData<
            TTensorType,
            TTensorStorageType,
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
            long double>(rTensor, rArray))
    {
        KRATOS_ERROR
            << "TensorsAdaptors cannot be assigned an numpy array with \""
            << rArray.dtype() << "\". They can be only set with numpy arrays having following dtypes:"
            << "\n\t numpy.bool"
            << "\n\t numpy.uint8"
            << "\n\t numpy.uint16"
            << "\n\t numpy.uint32"
            << "\n\t numpy.uint64"
            << "\n\t numpy.int8"
            << "\n\t numpy.int16"
            << "\n\t numpy.int32"
            << "\n\t numpy.int64"
            << "\n\t numpy.float32"
            << "\n\t numpy.float64"
            << "\n\t numpy.float128";
    }
}

template<class TDataType>
void AddBaseTensorAdaptor(
    pybind11::module& rModule,
    const std::string& rName)
{
    // add the base tensor adaptor
    using tensor_adaptor = TensorAdaptor<TDataType>;
    pybind11::class_<tensor_adaptor, typename tensor_adaptor::Pointer>(rModule, (rName + "Adaptor").c_str())
        .def(pybind11::init<const tensor_adaptor&, const bool>(), pybind11::arg("tensor_adaptor"), pybind11::arg("copy") = false)
        .def("Check", &tensor_adaptor::Check)
        .def("CollectData", &tensor_adaptor::CollectData)
        .def("StoreData", &tensor_adaptor::StoreData)
        .def("GetContainer", &tensor_adaptor::GetContainer)
        .def("HasContainer", &tensor_adaptor::HasContainer)
        .def("Shape", &tensor_adaptor::Shape)
        .def("DataShape", &tensor_adaptor::DataShape)
        .def("Size", &tensor_adaptor::Size)
        .def("__str__", PrintObject<tensor_adaptor>)
        .def("ViewData", &Detail::GetPybindArray<TensorAdaptor, TDataType>)
        .def("MoveData", &Detail::MovePybindArray<TensorAdaptor, TDataType>)
        .def("SetData", &Detail::SetPybindArray<TensorAdaptor, TDataType>, pybind11::arg("array").noconvert())
        .def_property("data",
            &Detail::GetPybindArray<TensorAdaptor, TDataType>,
            pybind11::cpp_function(
                &Detail::SetPybindArray<TensorAdaptor, TDataType>,
                pybind11::arg("self"),
                // no convert makes sure that the numpy arrays are
                // not converted, hence nothing will be copied. numpy
                // array will be passed as it is to the SetPybindArray
                // method.
                pybind11::arg("array").noconvert())
            )
    ;
}

template<class TDataType>
void AddCombinedTensorAdaptor(
    pybind11::module& rModule,
    const std::string& rName)
{
    using combined_ta_type = CombinedTensorAdaptor<TDataType>;
    pybind11::class_<combined_ta_type, typename combined_ta_type::Pointer, typename combined_ta_type::BaseType>(rModule, rName.c_str())
    .def(pybind11::init<const typename combined_ta_type::TensorAdaptorVectorType&, const bool>(), pybind11::arg("list_of_tensor_adaptors"), pybind11::arg("collect_and_store_recursively") = true) // reveling ctor
        .def(pybind11::init<const typename combined_ta_type::TensorAdaptorVectorType&, const unsigned int, const bool>(), pybind11::arg("list_of_tensor_adaptors"), pybind11::arg("axis"), pybind11::arg("collect_and_store_recursively") = true) // axis based ctor
        .def(pybind11::init<const combined_ta_type&, const bool, const bool>(), pybind11::arg("list_of_tensor_adaptors"), pybind11::arg("collect_and_store_recursively") = true, pybind11::arg("copy") = true)
        .def("GetTensorAdaptors", &combined_ta_type::GetTensorAdaptors)
        ;
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
        .def(py::init<ModelPart::NodesContainerType::Pointer, HistoricalVariableTensorAdaptor::VariablePointerType, const int>(), py::arg("container"), py::arg("variable"), py::arg("step_index") = 0)
        .def(py::init<ModelPart::NodesContainerType::Pointer, TensorAdaptorUtils::VariablePointerType, const std::vector<unsigned int>&, const int>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"), py::arg("step_index") = 0)
        .def(py::init<const HistoricalVariableTensorAdaptor::BaseType&, HistoricalVariableTensorAdaptor::VariablePointerType, const int, const bool>(), py::arg("tensor_adaptor"), py::arg("variable"), py::arg("step_index") = 0, py::arg("copy") = true)
        ;

    py::class_<VariableTensorAdaptor, VariableTensorAdaptor::Pointer, VariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "VariableTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"))
        .def(py::init<ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"))
        .def(py::init<ModelPart::GeometryContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::GeometryContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"))
        .def(py::init<ModelPart::MasterSlaveConstraintContainerType::Pointer, VariableTensorAdaptor::VariablePointerType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<ModelPart::MasterSlaveConstraintContainerType::Pointer, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"))
        .def(py::init<const VariableTensorAdaptor::BaseType&, VariableTensorAdaptor::VariablePointerType, const bool>(), py::arg("tensor_adaptor"), py::arg("variable"), py::arg("copy") = true)
        ;

    py::class_<GaussPointVariableTensorAdaptor, GaussPointVariableTensorAdaptor::Pointer, GaussPointVariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "GaussPointVariableTensorAdaptor")
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer>(), py::arg("container"), py::arg("variable"), py::arg("process_info"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer>(), py::arg("container"), py::arg("variable"), py::arg("process_info"))
        .def(py::init<const GaussPointVariableTensorAdaptor::BaseType&, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer, const bool>(), py::arg("tensor_adaptor"), py::arg("variable"), py::arg("process_info"), py::arg("copy") = true)
        ;

    py::class_<EquationIdsTensorAdaptor, EquationIdsTensorAdaptor::Pointer, EquationIdsTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "EquationIdsTensorAdaptor")
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, ProcessInfo::Pointer>(), py::arg("container"), py::arg("process_info"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, ProcessInfo::Pointer>(), py::arg("container"), py::arg("process_info"))
        .def(py::init<const EquationIdsTensorAdaptor::BaseType&, ProcessInfo::Pointer, const bool>(), py::arg("tensor_adaptor"), py::arg("process_info"), py::arg("copy") = true)
        ;

    py::class_<FlagsTensorAdaptor, FlagsTensorAdaptor::Pointer, FlagsTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "FlagsTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, const Flags&>(), py::arg("container"), py::arg("flag"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, const Flags&>(), py::arg("container"), py::arg("flag"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, const Flags&>(), py::arg("container"), py::arg("flag"))
        .def(py::init<const FlagsTensorAdaptor::BaseType&, const Flags&, const bool>(), py::arg("tensor_adaptor"), py::arg("flag"), py::arg("copy") = true)
        ;

    py::class_<NodePositionTensorAdaptor, NodePositionTensorAdaptor::Pointer, NodePositionTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "NodePositionTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, Globals::Configuration>(), py::arg("container"), py::arg("configuration"))
        .def(py::init<ModelPart::NodesContainerType::Pointer, Globals::Configuration, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("configuration"), py::arg("data_shape"))
        .def(py::init<const NodePositionTensorAdaptor::BaseType&, Globals::Configuration, const bool>(), py::arg("tensor_adaptor"), py::arg("configuration"), py::arg("copy") = true)
        ;
}

} // namespace Kratos::Python.