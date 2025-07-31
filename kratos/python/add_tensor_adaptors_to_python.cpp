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

template<class TDataType>
class TensorAdaptorTrampoline final : public TensorAdaptor<TDataType> {
public:
    using BaseType = TensorAdaptor<TDataType>;

    void Check() const override
    {
        PYBIND11_OVERRIDE_PURE(void,                            /*return type*/
                               BaseType,                        /*base type*/
                               Check                            /*function name*/
        );
    }

    void CollectData() override
    {
        PYBIND11_OVERRIDE_PURE(void,                            /*return type*/
                               BaseType,                        /*base type*/
                               CollectData                      /*function name*/
        );
    }

    void StoreData() override
    {
        PYBIND11_OVERRIDE_PURE(void,                            /*return type*/
                               BaseType,                        /*base type*/
                               StoreData                        /*function name*/
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
        const auto& casted_array = rArray.cast<pybind11::array_t<TPybindArrayType, pybind11::array::c_style>>();
        IndexPartition<IndexType>(rTensor.ViewData().size()).for_each([&casted_array, &rTensor](const auto Index) {
            rTensor.ViewData()[Index] = static_cast<TTensorStorageType>(casted_array.data()[Index]);
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
    // add the base TensorStorage
    using tensor_storage = TensorStorage<TDataType>;
    pybind11::class_<tensor_storage, typename tensor_storage::Pointer>(rModule, (rName + "Data").c_str())
        .def(pybind11::init<typename tensor_storage::ContainerPointerType, const DenseVector<unsigned int>&>(), pybind11::arg("container"), pybind11::arg("tensor_shape"))
        .def("Copy", &tensor_storage::Copy)
        .def("GetContainer", &tensor_storage::GetContainer)
        .def("Shape", &tensor_storage::Shape)
        .def("DataShape", &tensor_storage::DataShape)
        .def("Size", &tensor_storage::Size)
        .def("__str__", PrintObject<tensor_storage>)
        .def("ViewData", &Detail::GetPybindArray<TensorStorage, TDataType>)
        .def("MoveData", &Detail::MovePybindArray<TensorStorage, TDataType>)
        .def("SetData", &Detail::SetPybindArray<TensorStorage, TDataType>, pybind11::arg("array").noconvert())
        .def_property("data",
            &Detail::GetPybindArray<TensorStorage, TDataType>,
            pybind11::cpp_function(
                &Detail::SetPybindArray<TensorStorage, TDataType>,
                pybind11::arg("self"),
                // no convert makes sure that the numpy arrays are
                // not converted, hence nothing will be copied. numpy
                // array will be passed as it is to the SetPybindArray
                // method.
                pybind11::arg("array").noconvert())
            )
    ;

    // add the base tensor adaptor
    using tensor_adaptor = TensorAdaptor<TDataType>;
    pybind11::class_<tensor_adaptor, typename tensor_adaptor::Pointer>(rModule, (rName + "Adaptor").c_str())
        .def("Check", &tensor_adaptor::Check)
        .def("CollectData", &tensor_adaptor::CollectData)
        .def("StoreData", &tensor_adaptor::StoreData)
        .def("GetContainer", &tensor_adaptor::GetContainer)
        .def("GetStorage", pybind11::overload_cast<>(&tensor_adaptor::GetStorage))
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

} // namespace Detail

void AddTensorAdaptorsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto tensor_adaptor_sub_module = m.def_submodule("TensorAdaptors");
    Detail::AddBaseTensorAdaptor<bool>(tensor_adaptor_sub_module, "BoolTensor");
    Detail::AddBaseTensorAdaptor<int>(tensor_adaptor_sub_module, "IntTensor");
    Detail::AddBaseTensorAdaptor<double>(tensor_adaptor_sub_module, "DoubleTensor");

    py::class_<HistoricalVariableTensorAdaptor, HistoricalVariableTensorAdaptor::Pointer, HistoricalVariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "HistoricalVariableTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, HistoricalVariableTensorAdaptor::VariablePointerType, const int>(), py::arg("container"), py::arg("variable"), py::arg("step_index") = 0)
        .def(py::init<ModelPart::NodesContainerType::Pointer, TensorAdaptorUtils::VariablePointerType, const std::vector<unsigned int>&, const int>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"), py::arg("step_index") = 0)
        .def(py::init<TensorStorage<double>::Pointer, HistoricalVariableTensorAdaptor::VariablePointerType, const int>(), py::arg("tensor_storage"), py::arg("variable"), py::arg("step_index") = 0)
        ;

    pybind11::class_<VariableTensorAdaptor, VariableTensorAdaptor::Pointer, VariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "VariableTensorAdaptor")
        .def(py::init<VariableTensorAdaptor::ContainerPointerType, VariableTensorAdaptor::VariablePointerType>(), py::arg("container"), py::arg("variable"))
        .def(py::init<VariableTensorAdaptor::ContainerPointerType, VariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("variable"), py::arg("data_shape"))
        .def(py::init<TensorStorage<double>::Pointer, VariableTensorAdaptor::VariablePointerType>(), py::arg("tensor_storage"), py::arg("variable"))
        ;

    py::class_<GaussPointVariableTensorAdaptor, GaussPointVariableTensorAdaptor::Pointer, GaussPointVariableTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "GaussPointVariableTensorAdaptor")
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer>(), py::arg("container"), py::arg("variable"), py::arg("process_info"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, GaussPointVariableTensorAdaptor::VariablePointerType, ProcessInfo::Pointer>(), py::arg("container"), py::arg("variable"), py::arg("process_info"))
        ;

    pybind11::class_<EquationIdsTensorAdaptor, EquationIdsTensorAdaptor::Pointer, EquationIdsTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "EquationIdsTensorAdaptor")
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, ProcessInfo::Pointer>(), pybind11::arg("container"), pybind11::arg("process_info"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, ProcessInfo::Pointer>(), pybind11::arg("container"), pybind11::arg("process_info"))
        ;

    pybind11::class_<FlagsTensorAdaptor, FlagsTensorAdaptor::Pointer, FlagsTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "FlagsTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, const Flags&>(), pybind11::arg("container"), pybind11::arg("flag"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, const Flags&>(), pybind11::arg("container"), pybind11::arg("flag"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, const Flags&>(), pybind11::arg("container"), pybind11::arg("flag"))
        .def(py::init<TensorStorage<bool>::Pointer, const Flags&>(), pybind11::arg("tensor_storage"), pybind11::arg("flag"))
        ;

    pybind11::class_<NodePositionTensorAdaptor, NodePositionTensorAdaptor::Pointer, NodePositionTensorAdaptor::BaseType>(tensor_adaptor_sub_module, "NodePositionTensorAdaptor")
        .def(py::init<ModelPart::NodesContainerType::Pointer, Globals::Configuration>(), py::arg("container"), py::arg("configuration"))
        .def(py::init<ModelPart::NodesContainerType::Pointer, Globals::Configuration, const std::vector<unsigned int>&>(), py::arg("container"), py::arg("configuration"), py::arg("data_shape"))
        .def(py::init<TensorStorage<double>::Pointer, Globals::Configuration>(), py::arg("tensor_storage"), py::arg("configuration"))
        ;
}

} // namespace Kratos::Python.