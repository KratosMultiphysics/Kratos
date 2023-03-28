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
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

// Project includes
#include "containers/container_variable_data/container_data_io.h"
#include "containers/container_variable_data/container_variable_data.h"
#include "containers/container_variable_data/specialized_container_variable_data.h"

// Include base h
#include "add_container_variable_data_to_python.h"

namespace Kratos::Python
{

template <class TDataType>
using ContiguousNumpyArray = pybind11::array_t<
    TDataType,
    /*column-major*/ pybind11::array::c_style | /*cast instead of looking for a suitable overload*/ pybind11::array::forcecast
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

    std::vector<IndexType> c_shape(rShape.size() + 1);
    c_shape[0] = NumberOfEntities;
    std::copy(rShape.begin(), rShape.end(), c_shape.begin() + 1);

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

#define KRATOS_FORBIDDEN_CAST(METHOD_NAME, SELF_TYPE, CONST, DATA_TYPE)                           \
    .def(METHOD_NAME, [](SELF_TYPE& rSelf, CONST py::array_t<DATA_TYPE>& rData) {                 \
        KRATOS_ERROR << "Unsupported numpy array is passed. Please change "                       \
                     << "it to dtype = numpy.float64. [ data_type = " << #DATA_TYPE << " ].\n"; })

template<class TContainerType>
void AddContainerVariableDataToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_variable_data_holder_base = ContainerVariableData<TContainerType>;
    py::class_<container_variable_data_holder_base, typename container_variable_data_holder_base::Pointer>(m, rName.c_str())
        .def("CopyFrom", &container_variable_data_holder_base::CopyFrom, py::arg("origin_container_data"))
        .def("MoveFrom", [](container_variable_data_holder_base& rSelf, py::array_t<double>& rData){
            KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

            std::vector<int> shape(rData.ndim() - 1);
            std::copy(rData.shape() + 1, rData.shape() + rData.ndim(), shape.begin());

            rSelf.MoveFrom(rData.mutable_data(),
                           rData.shape()[0],
                           shape.data(),
                           shape.size());

            std::for_each(rData.mutable_data(), rData.mutable_data() + rData.size(), [](double& rValue) { rValue += 1;});
        })
        KRATOS_FORBIDDEN_CAST("MoveFrom", container_variable_data_holder_base, , float)
        KRATOS_FORBIDDEN_CAST("MoveFrom", container_variable_data_holder_base, , long double)
        KRATOS_FORBIDDEN_CAST("MoveFrom", container_variable_data_holder_base, , int)
        KRATOS_FORBIDDEN_CAST("MoveFrom", container_variable_data_holder_base, , long int)
        KRATOS_FORBIDDEN_CAST("MoveFrom", container_variable_data_holder_base, , unsigned int)
        KRATOS_FORBIDDEN_CAST("MoveFrom", container_variable_data_holder_base, , long unsigned int)
        .def("Read", &container_variable_data_holder_base::Read, py::arg("starting_value"), py::arg("number_of_entities"), py::arg("starting_value_of_shape"), py::arg("shape_size"))
        .def("GetModelPart", py::overload_cast<>(&container_variable_data_holder_base::GetModelPart), py::return_value_policy::reference)
        .def("GetContainer", py::overload_cast<>(&container_variable_data_holder_base::GetContainer), py::return_value_policy::reference)
        .def("PrintData", &container_variable_data_holder_base::PrintData)
        .def("__str__", &container_variable_data_holder_base::Info)
        ;
}

template<class TContainerType, class TContainerDataIOTag>
void AddSpecializedContainerVariableDataToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_type = SpecializedContainerVariableData<TContainerType, ContainerDataIO<TContainerDataIOTag>>;
    py::class_<container_type, typename container_type::Pointer, ContainerVariableData<TContainerType>>(m, rName.c_str())
        .def(py::init<ModelPart&>(), py::arg("model_part"), py::doc("Creates a new container data object with model_part."))
        .def(py::init<const container_type&>(), py::arg("other_container_data_to_copy_from"), py::doc("Creates a new same type container data object by copying data from other_container_data_to_copy_from."))
        .def(py::init<const typename container_type::BaseType&>(), py::arg("other_container_data_to_copy_from"), py::doc("Creates a new destination type container data object by copying data from compatible other_container_data_to_copy_from."))
        .def("Evaluate", [](const container_type& rSelf){
            const auto& r_shape = rSelf.GetShape();
            auto array = AllocateNumpyArray<double>(rSelf.GetContainer().size(), r_shape);

            std::vector<int> shape(r_shape.size());
            std::transform(r_shape.begin(), r_shape.end(), shape.begin(), [](const IndexType Value) -> int { return Value; });

            rSelf.Evaluate(array.mutable_data(),
                           rSelf.GetContainer().size(),
                           shape.data(),
                           shape.size());

            return array;
        })
        .def("Evaluate", &container_type::template Evaluate<double>, py::arg("scalar_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 4>>, py::arg("Array4_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 6>>, py::arg("Array6_variable"))
        .def("Evaluate", &container_type::template Evaluate<array_1d<double, 9>>, py::arg("Array9_variable"))
        .def("Evaluate", &container_type::template Evaluate<Vector>, py::arg("Vector_variable"))
        .def("Evaluate", &container_type::template Evaluate<Matrix>, py::arg("Matrix_variable"))
        .def("Read", &container_type::template Read<double>, py::arg("scalar_variable"))
        .def("Read", [](container_type& rSelf, const py::array_t<double>& rData){
            KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

            std::vector<int> shape(rData.ndim() - 1);
            std::copy(rData.shape() + 1, rData.shape() + rData.ndim(), shape.begin());

            rSelf.Read(rData.data(),
                       rData.shape()[0],
                       shape.data(),
                       shape.size());
        })
        KRATOS_FORBIDDEN_CAST("Read", container_type, const, float)
        KRATOS_FORBIDDEN_CAST("Read", container_type, const, long double)
        KRATOS_FORBIDDEN_CAST("Read", container_type, const, int)
        KRATOS_FORBIDDEN_CAST("Read", container_type, const, long int)
        KRATOS_FORBIDDEN_CAST("Read", container_type, const, unsigned int)
        KRATOS_FORBIDDEN_CAST("Read", container_type, const, long unsigned int)
        .def("Read", &container_type::template Read<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("Read", &container_type::template Read<array_1d<double, 4>>, py::arg("Array4_variable"))
        .def("Read", &container_type::template Read<array_1d<double, 6>>, py::arg("Array6_variable"))
        .def("Read", &container_type::template Read<array_1d<double, 9>>, py::arg("Array9_variable"))
        .def("Read", &container_type::template Read<Vector>, py::arg("Vector_variable"))
        .def("Read", &container_type::template Read<Matrix>, py::arg("Matrix_variable"))
        .def("SetData", &container_type::template SetData<double>, py::arg("scalar_value"))
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
        .def(py::self +  py::self)
        .def(py::self += py::self)
        .def(py::self +  float())
        .def(py::self += float())
        .def(py::self -  py::self)
        .def(py::self -= py::self)
        .def(py::self -  float())
        .def(py::self -= float())
        .def(py::self *  py::self)
        .def(py::self *= py::self)
        .def(py::self *  float())
        .def(py::self *= float())
        .def(py::self /  py::self)
        .def(py::self /= py::self)
        .def(py::self /  float())
        .def(py::self /= float())
        .def("__pow__", [](container_type& rSelf, const container_type& rInput) { container_type result(rSelf.GetModelPart()); result = rSelf.Pow(rInput); return result; })
        .def("__ipow__", [](container_type& rSelf, const container_type& rInput) { rSelf = rSelf.Pow(rInput); return rSelf; })
        .def("__pow__", [](container_type& rSelf, const double Value) { container_type result(rSelf.GetModelPart()); result = rSelf.Pow(Value); return result; })
        .def("__ipow__", [](container_type& rSelf, const double Value) { rSelf = rSelf.Pow(Value); return rSelf; })
        .def("__neg__", [](container_type& rSelf) { return rSelf.operator*(-1.0); })
        ;
}

void  AddContainerVariableDataToPython(pybind11::module& m)
{
    auto sub_module = m.def_submodule("ContainerVariableData");

    AddContainerVariableDataToPython<ModelPart::NodesContainerType>(sub_module, "NodalVariableData");
    AddContainerVariableDataToPython<ModelPart::ConditionsContainerType>(sub_module, "ConditionVariableData");
    AddContainerVariableDataToPython<ModelPart::ElementsContainerType>(sub_module, "ElementVariableData");

    AddSpecializedContainerVariableDataToPython<ModelPart::NodesContainerType, ContainerDataIOTags::Historical>(sub_module, "HistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "NodalNonHistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ConditionNonHistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ElementNonHistoricalVariableData");
}

#undef KRATOS_FORBIDDEN_CAST

} // namespace Kratos::Python
