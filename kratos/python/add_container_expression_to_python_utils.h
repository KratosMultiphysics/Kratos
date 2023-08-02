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
#include <numeric>

// External includes
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// Project includes
#include "expression/container_expression.h"
#include "expression/specialized_container_expression.h"
#include "expression/c_array_expression_io.h"

namespace Kratos::Python
{

template <class TDataType>
using ContiguousNumpyArray = pybind11::array_t<
    TDataType,
    /*column-major*/ pybind11::array::c_style
>;

template <class TDataType>
pybind11::array_t<TDataType> AllocateNumpyArray(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    const IndexType size = std::accumulate(rShape.begin(),
                                           rShape.end(),
                                           1UL,
                                           [](IndexType left, IndexType right) {return left * right;});

    TDataType* array = new TDataType[size * NumberOfEntities];
    pybind11::capsule release(array, [](void* a) {
        delete[] reinterpret_cast<TDataType*>(a);
    });

    // we allocate one additional dimension for the shape to hold the
    // number of entities. So if the shape within kratos is [2, 3] with 70 entities,
    // then the final shape of the numpy array will be [70,2,3] so it can keep the
    // existing shape preserved for each entitiy.
    std::vector<IndexType> c_shape(rShape.size() + 1);
    c_shape[0] = NumberOfEntities;
    std::copy(rShape.begin(), rShape.end(), c_shape.begin() + 1);

    // we have to allocate number of bytes to be skipped to reach next element
    // in each dimension. So this is calculated backwards.
    std::vector<IndexType> strides(c_shape.size());
    IndexType stride_items = 1;
    for (int i = c_shape.size() - 1; i >= 0; --i) {
        strides[i] = sizeof(double) * stride_items;
        stride_items *= c_shape[i];
    }

    return ContiguousNumpyArray<TDataType>(
        c_shape,
        strides,
        array,
        release
    );
}

template <class TDataType>
pybind11::array_t<TDataType> MakeNumpyArray(
    TDataType const* pBegin,
    TDataType const* pEnd,
    const std::vector<IndexType>& rShape)
{
    auto array = AllocateNumpyArray<TDataType>(rShape);
    KRATOS_ERROR_IF_NOT(std::distance(pBegin, pEnd) == array.size()) << "Size mismatch.";
    std::copy(pBegin,
              pEnd,
              array.mutable_data());
    return array;
}

template<class TContainerType>
void AddContainerExpressionToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_expression_holder_base = ContainerExpression<TContainerType>;
    py::class_<container_expression_holder_base, typename container_expression_holder_base::Pointer>(m, rName.c_str())
        .def(py::init<ModelPart&>(), py::arg("model_part"))
        .def("CopyFrom", &container_expression_holder_base::CopyFrom, py::arg("origin_container_expression"))
        .def("MoveFrom", [](container_expression_holder_base& rSelf, py::array_t<int>& rData){
            KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

            // dimension of the numpy array is always one dimension greater than the kratos stored dimension for each
            // entity. That is because, first dimension of the numpy array shows how many entities are there
            // in the numpy array to be read in. If the numpy array dimension is [45, 3, 4] then it shows
            // there are 45 entities each having matrices of shape [3, 4].
            std::vector<int> shape(rData.ndim() - 1);
            std::copy(rData.shape() + 1, rData.shape() + rData.ndim(), shape.begin());

            rSelf.MoveFrom(rData.mutable_data(),
                           rData.shape()[0],
                           shape.data(),
                           shape.size());
        }, py::arg("numpy_int_array").noconvert())
        .def("MoveFrom", [](container_expression_holder_base& rSelf, py::array_t<double>& rData){
            KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

            // dimension of the numpy array is always one dimension greater than the kratos stored dimension for each
            // entity. That is because, first dimension of the numpy array shows how many entities are there
            // in the numpy array to be read in. If the numpy array dimension is [45, 3, 4] then it shows
            // there are 45 entities each having matrices of shape [3, 4].
            std::vector<int> shape(rData.ndim() - 1);
            std::copy(rData.shape() + 1, rData.shape() + rData.ndim(), shape.begin());

            rSelf.MoveFrom(rData.mutable_data(),
                           rData.shape()[0],
                           shape.data(),
                           shape.size());
        }, py::arg("numpy_double_array").noconvert())
        .def("SetExpression", &container_expression_holder_base::SetExpression)
        .def("HasExpression", &container_expression_holder_base::HasExpression)
        .def("GetExpression", &container_expression_holder_base::pGetExpression)
        .def("GetModelPart", py::overload_cast<>(&container_expression_holder_base::GetModelPart), py::return_value_policy::reference)
        .def("GetContainer", py::overload_cast<>(&container_expression_holder_base::GetContainer), py::return_value_policy::reference)
        .def("GetItemShape", &container_expression_holder_base::GetItemShape)
        .def("GetItemComponentCount", &container_expression_holder_base::GetItemComponentCount)
        .def("Slice",
             &container_expression_holder_base::Slice,
             py::arg("offset"),
             py::arg("stride"))
        .def("Reshape",
             &container_expression_holder_base::Reshape,
             py::arg("new_shape"))
        .def("Comb",
             [](container_expression_holder_base& rSelf,
                const container_expression_holder_base& rOther)
                {return rSelf.Comb(rOther);},
             py::arg("other"))
        .def("Comb",
             [](container_expression_holder_base& rSelf,
                const std::vector<typename container_expression_holder_base::Pointer>& rOthers)
                {return rSelf.Comb(rOthers);},
             py::arg("others"))
        .def("Evaluate", [](const container_expression_holder_base& rSelf){
            const auto& r_shape = rSelf.GetItemShape();
            auto array = AllocateNumpyArray<double>(rSelf.GetContainer().size(), r_shape);
            CArrayExpressionIO::Write(rSelf, array.mutable_data(), array.size());
            return array;
        })
        .def("Clone", &container_expression_holder_base::Clone)
        .def("Scale", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther){auto copy = rSelf; copy.SetExpression(Scale(rSelf.pGetExpression(), rOther.pGetExpression())); return copy;})
        .def("__add__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf + rOther; })
        .def("__iadd__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf + rOther; return rSelf; })
        .def("__add__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf + Value; })
        .def("__iadd__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf + Value; return rSelf; })
        .def("__sub__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf - rOther; })
        .def("__isub__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf - rOther; return rSelf; })
        .def("__sub__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf - Value; })
        .def("__isub__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf - Value; return rSelf; })
        .def("__mul__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf * rOther; })
        .def("__imul__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf * rOther; return rSelf; })
        .def("__mul__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf * Value; })
        .def("__imul__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf * Value; return rSelf; })
        .def("__truediv__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf / rOther; })
        .def("__itruediv__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf / rOther; return rSelf; })
        .def("__truediv__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf / Value; })
        .def("__itruediv__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf / Value; return rSelf; })
        .def("__pow__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rInput) { container_expression_holder_base result(rSelf.GetModelPart()); result = Power(rSelf, rInput); return result; })
        .def("__ipow__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rInput) { rSelf = Power(rSelf, rInput); return rSelf; })
        .def("__pow__", [](container_expression_holder_base& rSelf, const double Value) { container_expression_holder_base result(rSelf.GetModelPart()); result = Power(rSelf, Value); return result; })
        .def("__ipow__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = Power(rSelf, Value); return rSelf; })
        .def("__neg__", [](container_expression_holder_base& rSelf) { return rSelf *= -1.0; })
        .def("PrintData", &container_expression_holder_base::PrintData)
        .def("__str__", &container_expression_holder_base::Info)
        ;
}

template<class TContainerType, class TContainerDataIOTag>
void AddSpecializedContainerExpressionToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_type = SpecializedContainerExpression<TContainerType, ContainerDataIO<TContainerDataIOTag>>;
    py::class_<container_type, typename container_type::Pointer, ContainerExpression<TContainerType>>(m, rName.c_str())
        .def(py::init<ModelPart&>(), py::arg("model_part"), py::doc("Creates a new container expression object with model_part."))
        .def(py::init<const container_type&>(), py::arg("other_container_expression_to_copy_from"), py::doc("Creates a new same type container expression object by copying data from other_container_expression_to_copy_from."))
        .def(py::init<const typename container_type::BaseType&>(), py::arg("other_container_expression_to_copy_from"), py::doc("Creates a new destination type container expression object by copying data from compatible other_container_expression_to_copy_from."))
        .def("Evaluate", [](const container_type& rSelf){
            const auto& r_shape = rSelf.GetItemShape();
            auto array = AllocateNumpyArray<double>(rSelf.GetContainer().size(), r_shape);
            CArrayExpressionIO::Write(rSelf, array.mutable_data(), array.size());
            return array;
        })
        .def("Evaluate", &container_type::template Evaluate<int>, py::arg("scalar_int_variable"))
        .def("Evaluate", &container_type::template Evaluate<double>, py::arg("scalar_double_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 4>>, py::arg("Array4_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 6>>, py::arg("Array6_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 9>>, py::arg("Array9_variable"))
        .def("Evaluate", &container_type::template Evaluate<Vector>, py::arg("Vector_variable"))
        .def("Evaluate", &container_type::template Evaluate<Matrix>, py::arg("Matrix_variable"))
        .def("Read", &container_type::template Read<int>, py::arg("scalar_int_variable"))
        .def("Read", &container_type::template Read<double>, py::arg("scalar_double_variable"))
        .def("Read", [](container_type& rSelf, const py::array_t<int>& rData){
            KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

            std::vector<int> shape(rData.ndim() - 1);
            std::copy(rData.shape() + 1, rData.shape() + rData.ndim(), shape.begin());

            rSelf.Read(rData.data(),
                       rData.shape()[0],
                       shape.data(),
                       shape.size());
        }, py::arg("numpy_int_array").noconvert())
        .def("Read", [](container_type& rSelf, const py::array_t<double>& rData){
            KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

            std::vector<int> shape(rData.ndim() - 1);
            std::copy(rData.shape() + 1, rData.shape() + rData.ndim(), shape.begin());

            rSelf.Read(rData.data(),
                       rData.shape()[0],
                       shape.data(),
                       shape.size());
        }, py::arg("numpy_double_array").noconvert())
        .def("Read", &container_type::template Read<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("Read", &container_type::template Read<array_1d<double, 4>>, py::arg("Array4_variable"))
        .def("Read", &container_type::template Read<array_1d<double, 6>>, py::arg("Array6_variable"))
        .def("Read", &container_type::template Read<array_1d<double, 9>>, py::arg("Array9_variable"))
        .def("Read", &container_type::template Read<Vector>, py::arg("Vector_variable"))
        .def("Read", &container_type::template Read<Matrix>, py::arg("Matrix_variable"))
        .def("SetData", &container_type::template SetData<int>, py::arg("scalar_int_value"))
        .def("SetData", &container_type::template SetData<double>, py::arg("scalar_double_value"))
        .def("SetData", &container_type::template SetData<array_1d<double, 3>>, py::arg("Array3_value"))
        .def("SetData", &container_type::template SetData<array_1d<double, 4>>, py::arg("Array4_value"))
        .def("SetData", &container_type::template SetData<array_1d<double, 6>>, py::arg("Array6_value"))
        .def("SetData", &container_type::template SetData<array_1d<double, 9>>, py::arg("Array9_value"))
        .def("SetData", &container_type::template SetData<Vector>, py::arg("Vector_value"))
        .def("SetData", &container_type::template SetData<Matrix>, py::arg("Matrix_value"))
        .def("SetZero", &container_type::template SetZero<double>, py::arg("scalar_variable"))
        .def("SetZero", &container_type::template SetZero<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("SetZero", &container_type::template SetZero<array_1d<double, 4>>, py::arg("Array4_variable"))
        .def("SetZero", &container_type::template SetZero<array_1d<double, 6>>, py::arg("Array6_variable"))
        .def("SetZero", &container_type::template SetZero<array_1d<double, 9>>, py::arg("Array9_variable"))
        .def("SetZero", &container_type::template SetZero<Vector>, py::arg("Vector_variable"))
        .def("SetZero", &container_type::template SetZero<Matrix>, py::arg("Matrix_variable"))
        .def("Clone", &container_type::Clone)
        .def("Slice", &container_type::Slice, py::arg("offset"), py::arg("stride"))
        .def("Reshape", [](const container_type& rSelf, const std::vector<IndexType>& rShape) { return rSelf.Reshape(rShape); }, py::arg("shape"))
        .def("Comb", [](const container_type& rSelf, const typename container_type::BaseType& rOther) { return rSelf.Comb(rOther); }, py::arg("other_container_expression_to_combine_with"))
        .def("Comb", [](const container_type& rSelf, const std::vector<typename container_type::BaseType::Pointer>& rListOfOthersContainerExpressions) { return rSelf.Comb(rListOfOthersContainerExpressions); }, py::arg("other_container_expressions_list_to_combine_with"))
        .def("Scale", [](const container_type& rSelf, const container_type& rOther){auto copy = rSelf; copy.SetExpression(Scale(rSelf.pGetExpression(), rOther.pGetExpression())); return copy;})
        .def("__add__", [](const container_type& rSelf, const container_type& rOther) { return rSelf + rOther; })
        .def("__iadd__", [](container_type& rSelf, const container_type& rOther) { rSelf = rSelf + rOther; return rSelf; })
        .def("__add__", [](const container_type& rSelf, const double Value) { return rSelf + Value; })
        .def("__iadd__", [](container_type& rSelf, const double Value) { rSelf = rSelf + Value; return rSelf; })
        .def("__sub__", [](const container_type& rSelf, const container_type& rOther) { return rSelf - rOther; })
        .def("__isub__", [](container_type& rSelf, const container_type& rOther) { rSelf = rSelf - rOther; return rSelf; })
        .def("__sub__", [](const container_type& rSelf, const double Value) { return rSelf - Value; })
        .def("__isub__", [](container_type& rSelf, const double Value) { rSelf = rSelf - Value; return rSelf; })
        .def("__mul__", [](const container_type& rSelf, const container_type& rOther) { return rSelf * rOther; })
        .def("__imul__", [](container_type& rSelf, const container_type& rOther) { rSelf = rSelf * rOther; return rSelf; })
        .def("__mul__", [](const container_type& rSelf, const double Value) { return rSelf * Value; })
        .def("__imul__", [](container_type& rSelf, const double Value) { rSelf = rSelf * Value; return rSelf; })
        .def("__truediv__", [](const container_type& rSelf, const container_type& rOther) { return rSelf / rOther; })
        .def("__itruediv__", [](container_type& rSelf, const container_type& rOther) { rSelf = rSelf / rOther; return rSelf; })
        .def("__truediv__", [](const container_type& rSelf, const double Value) { return rSelf / Value; })
        .def("__itruediv__", [](container_type& rSelf, const double Value) { rSelf = rSelf / Value; return rSelf; })
        .def("__pow__", [](container_type& rSelf, const container_type& rInput) { container_type result(rSelf.GetModelPart()); result = rSelf.Power(rInput); return result; })
        .def("__ipow__", [](container_type& rSelf, const container_type& rInput) { rSelf = rSelf.Power(rInput); return rSelf; })
        .def("__pow__", [](container_type& rSelf, const double Value) { container_type result(rSelf.GetModelPart()); result = rSelf.Power(Value); return result; })
        .def("__ipow__", [](container_type& rSelf, const double Value) { rSelf = rSelf.Power(Value); return rSelf; })
        .def("__neg__", [](container_type& rSelf) { return rSelf.operator*(-1.0); })
        ;
}

} // namespace Kratos::Python
