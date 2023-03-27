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
#include <pybind11/operators.h>

// Project includes

// Application includes
#include "containers/container_variable_data/container_data_io.h"
#include "containers/container_variable_data/container_variable_data.h"
#include "containers/container_variable_data/specialized_container_variable_data.h"

// Include base h
#include "add_container_variable_data_to_python.h"

namespace Kratos::Python
{

template<class TContainerType>
void AddContainerVariableDataToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_variable_data_holder_base = ContainerVariableData<TContainerType>;
    py::class_<container_variable_data_holder_base, typename container_variable_data_holder_base::Pointer>(m, rName.c_str())
        .def("CopyDataFrom", &container_variable_data_holder_base::CopyDataFrom, py::arg("origin_container_data"))
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
        .def("AssignData", &container_type::template AssignData<double>, py::arg("scalar_variable"))
        .def("AssignData", &container_type::template AssignData<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("AssignData", &container_type::template AssignData<array_1d<double, 4>>, py::arg("Array4_variable"))
        .def("AssignData", &container_type::template AssignData<array_1d<double, 6>>, py::arg("Array6_variable"))
        .def("AssignData", &container_type::template AssignData<array_1d<double, 9>>, py::arg("Array9_variable"))
        .def("AssignData", &container_type::template AssignData<Vector>, py::arg("Vector_variable"))
        .def("AssignData", &container_type::template AssignData<Matrix>, py::arg("Matrix_variable"))
        .def("ReadData", &container_type::template ReadData<double>, py::arg("scalar_variable"))
        .def("ReadData", &container_type::template ReadData<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("ReadData", &container_type::template ReadData<array_1d<double, 4>>, py::arg("Array4_variable"))
        .def("ReadData", &container_type::template ReadData<array_1d<double, 6>>, py::arg("Array6_variable"))
        .def("ReadData", &container_type::template ReadData<array_1d<double, 9>>, py::arg("Array9_variable"))
        .def("ReadData", &container_type::template ReadData<Vector>, py::arg("Vector_variable"))
        .def("ReadData", &container_type::template ReadData<Matrix>, py::arg("Matrix_variable"))
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

}
