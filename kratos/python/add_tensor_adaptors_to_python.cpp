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
#include "numpy_utils.h"
#include "utilities/parallel_utilities.h"

// Tensor adaptors
#include "tensor_adaptors/tensor_adaptor.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"

// Include base h
#include "add_tensor_adaptors_to_python.h"

namespace Kratos::Python {

namespace Detail {

using PybindArrayType = std::variant<
                            pybind11::array_t<int, pybind11::array::c_style>,
                            pybind11::array_t<bool, pybind11::array::c_style>,
                            pybind11::array_t<double, pybind11::array::c_style>
                        >;


class TensorAdaptorTrampoline final : public TensorAdaptor {
public:
    using BaseType = TensorAdaptor;

    void CollectData() override
    {
        PYBIND11_OVERRIDE_PURE(void,                /*return type*/
                               BaseType,            /*base type*/
                               CollectData          /*function name*/
        );
    }

    void StoreData() override
    {
        PYBIND11_OVERRIDE_PURE(void,                /*return type*/
                               BaseType,            /*base type*/
                               StoreData            /*function name*/
        );
    }

    ContainerType GetContainer() const override
    {
        PYBIND11_OVERRIDE_PURE(ContainerType,       /*return type*/
                               BaseType,            /*base type*/
                               GetContainer         /*function name*/
        );
    }

    std::string Info() const override
    {
        PYBIND11_OVERRIDE_PURE(std::string,         /*return type*/
                               BaseType,            /*base type*/
                               Info                 /*function name*/
        );
    }
}; // class ExpressionTrampoline


PybindArrayType GetNumpyArray(TensorAdaptor& rTensorAdaptor)
{
    auto view_data = rTensorAdaptor.ViewData();

    return std::visit([&rTensorAdaptor](auto& rTensorData) -> PybindArrayType {
        using value_type = typename std::remove_cv_t<std::decay_t<decltype(rTensorData)>>::value_type;
        using numpy_array_type = pybind11::array_t<value_type, pybind11::array::c_style>;

        const auto& r_shape = rTensorAdaptor.Shape();

        std::vector<std::size_t> c_shape(r_shape.size());
        std::copy(r_shape.begin(), r_shape.end(), c_shape.begin());
        std::vector<std::size_t> strides(c_shape.size());

        std::size_t stride_items = 1;
        for (int i = c_shape.size() - 1; i >= 0; --i) {
            strides[i] = sizeof(value_type) * stride_items;
            stride_items *= c_shape[i];
        }

        if (rTensorAdaptor.Size() == 0) {
            return numpy_array_type(c_shape);
        } else {
            // do nothing in the release of the numpy array since the ownership is not passed
            // the ownership is kept with the TensorAdaptor.
            pybind11::capsule release(rTensorData.data(), [](void* a) {});

            return numpy_array_type(
                c_shape,
                strides,
                rTensorData.data(),
                release
            );
        }
    }, view_data);
}

template<class... TArgs>
void SetNumpyArray(
    TensorAdaptor& rTensorAdaptor,
    std::variant<pybind11::array_t<TArgs, pybind11::array::c_style> const *...> pArray)
{
    auto view_data = rTensorAdaptor.ViewData();

    std::visit([&rTensorAdaptor](const auto pArray, auto& rTensorData) {
        using value_type = typename std::remove_cv_t<std::decay_t<decltype(rTensorData)>>::value_type;
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
        IndexPartition<IndexType>(rTensorData.size()).for_each([&r_array, &rTensorData](const auto Index) {
            rTensorData[Index] = static_cast<value_type>(r_array.data()[Index]);
        });
    }, pArray, view_data);
}

} // namespace Detail


void AddTensorAdaptorsToPython(pybind11::module& m)
{
    auto tensor_adaptor_sub_module = m.def_submodule("TensorAdaptors");

    // add the base tensor adaptor
    pybind11::class_<TensorAdaptor, TensorAdaptor::Pointer>(tensor_adaptor_sub_module, "TensorAdaptor")
        .def("CollectData", &TensorAdaptor::CollectData)
        .def("StoreData", &TensorAdaptor::StoreData)
        .def("GetContainer", &TensorAdaptor::GetContainer)
        .def("Shape", &TensorAdaptor::Shape)
        .def("Size", &TensorAdaptor::Size)
        .def("ViewData", &Detail::GetNumpyArray)
        .def_property("data",
                      &Detail::GetNumpyArray,
                      pybind11::cpp_function(
                        &Detail::SetNumpyArray<
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
        .def("MoveData", [](TensorAdaptor& rSelf){
                auto moved_data = rSelf.MoveData();

                return std::visit([&rSelf](auto& rTensorData) -> Detail::PybindArrayType {
                    using value_type = typename std::remove_cv_t<std::decay_t<decltype(rTensorData)>>::value_type;
                    using numpy_array_type = pybind11::array_t<value_type, pybind11::array::c_style>;

                    const auto& r_shape = rSelf.Shape();

                    std::vector<std::size_t> c_shape(r_shape.size());
                    std::copy(r_shape.begin(), r_shape.end(), c_shape.begin());
                    std::vector<std::size_t> strides(c_shape.size());

                    std::size_t stride_items = 1;
                    for (int i = c_shape.size() - 1; i >= 0; --i) {
                        strides[i] = sizeof(value_type) * stride_items;
                        stride_items *= c_shape[i];
                    }

                    // this method transfers the ownership of the data to numpy.
                    // Since the DenseVector does not allow moving out the data without getting a call
                    // to destroy DenseVector underlying data at the destructor, we are creating a copy of the
                    // data by moving it to a temp so the underlying data container within the TensorAdaptor is cleared.
                    value_type* array = new value_type[rTensorData.size()];
                    std::copy(rTensorData.begin(), rTensorData.end(), array);

                    // now we add the release to clear the data.

                    if (rSelf.Size() == 0) {
                        return std::move(numpy_array_type(c_shape));
                    } else {
                        pybind11::capsule release(array, [](void* a) {
                            delete[] reinterpret_cast<value_type*>(a);
                        });

                        return std::move(numpy_array_type(
                            c_shape,
                            strides,
                            array,
                            release
                        ));
                    }
                }, moved_data);
        })
        .def("__str__", PrintObject<TensorAdaptor>);
    ;

    pybind11::class_<VariableTensorAdaptor, VariableTensorAdaptor::Pointer, TensorAdaptor>(tensor_adaptor_sub_module, "VariableTensorAdaptor")
        // nodal historical
        .def(pybind11::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariableType, const int>(), pybind11::arg("nodes"), pybind11::arg("variable"), pybind11::arg("step_index"))
        .def(pybind11::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariableType, const int, const std::vector<int>&>(), pybind11::arg("nodes"), pybind11::arg("variable"), pybind11::arg("step_index"), pybind11::arg("shape"))
        // non historical
        .def(pybind11::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariableType>(), pybind11::arg("nodes"), pybind11::arg("variable"))
        .def(pybind11::init<ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&>(), pybind11::arg("nodes"), pybind11::arg("variable"), pybind11::arg("shape"))
        .def(pybind11::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType>(), pybind11::arg("conditions"), pybind11::arg("variable"))
        .def(pybind11::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&>(), pybind11::arg("conditions"), pybind11::arg("variable"), pybind11::arg("shape"))
        .def(pybind11::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType>(), pybind11::arg("elements"), pybind11::arg("variable"))
        .def(pybind11::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&>(), pybind11::arg("elements"), pybind11::arg("variable"), pybind11::arg("shape"))
        .def(pybind11::init<ModelPart::GeometriesMapType::Pointer, VariableTensorAdaptor::VariableType>(), pybind11::arg("geometries"), pybind11::arg("variable"))
        .def(pybind11::init<ModelPart::GeometriesMapType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&>(), pybind11::arg("geometries"), pybind11::arg("variable"), pybind11::arg("shape"))
        .def(pybind11::init<ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariableType>(), pybind11::arg("properties"), pybind11::arg("variable"))
        .def(pybind11::init<ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&>(), pybind11::arg("properties"), pybind11::arg("variable"), pybind11::arg("shape"))
        // gauss point io
        .def(pybind11::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&>(), pybind11::arg("conditions"), pybind11::arg("variable"), pybind11::arg("process_info"))
        .def(pybind11::init<ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&, const std::vector<int>&>(), pybind11::arg("conditions"), pybind11::arg("variable"), pybind11::arg("process_info"), pybind11::arg("shape"))
        .def(pybind11::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&>(), pybind11::arg("elements"), pybind11::arg("variable"), pybind11::arg("process_info"))
        .def(pybind11::init<ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&, const std::vector<int>&>(), pybind11::arg("elements"), pybind11::arg("variable"), pybind11::arg("process_info"), pybind11::arg("shape"))
        ;
}

} // namespace Kratos::Python.