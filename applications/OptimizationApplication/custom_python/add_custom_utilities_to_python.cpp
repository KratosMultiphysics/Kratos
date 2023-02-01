//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes
#include <pybind11/stl.h>
#include <pybind11/operators.h>

// Project includes
#include "includes/model_part.h"
#include "includes/data_communicator.h"

// Application includes
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/container_variable_data_holder.h"
#include "custom_utilities/container_variable_data_holder_utils.h"

// Include base h
#include "add_custom_utilities_to_python.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace Kratos {
namespace Python {

template<class ContainerVariableDataHolderType>
void AddContainerTypeToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_type = ContainerVariableDataHolder<ContainerVariableDataHolderType>;
    py::class_<container_type, typename container_type::Pointer, ContainerVariableDataHolderBase>(m, rName.c_str())
        .def(py::init<ModelPart&>(), py::arg("model_part"), py::doc("Creates a new container data object with model_part."))
        .def(py::init<const container_type&>(), py::arg("other_container_data_to_copy_from"), py::doc("Creates a new same type container data object by copying data from other_container_data_to_copy_from."))
        .def(py::init<const ContainerVariableDataHolderBase&>(), py::arg("other_container_data_to_copy_from"), py::doc("Creates a new destination type container data object by copying data from compatible other_container_data_to_copy_from."))
        .def("AssignDataToContainerVariable", &container_type::template AssignDataToContainerVariable<double>, py::arg("scalar_variable"))
        .def("AssignDataToContainerVariable", &container_type::template AssignDataToContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("ReadDataFromContainerVariable", &container_type::template ReadDataFromContainerVariable<double>, py::arg("scalar_variable"))
        .def("ReadDataFromContainerVariable", &container_type::template ReadDataFromContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("SetDataForContainerVariable", &container_type::template SetDataForContainerVariable<double>, py::arg("scalar_variable"), py::arg("scalar_value"))
        .def("SetDataForContainerVariable", &container_type::template SetDataForContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"), py::arg("Array3_value"))
        .def("SetDataForContainerVariableToZero", &container_type::template SetDataForContainerVariableToZero<double>, py::arg("scalar_variable"))
        .def("SetDataForContainerVariableToZero", &container_type::template SetDataForContainerVariableToZero<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("Clone", &container_type::Clone)
        .def("GetContainer", py::overload_cast<>(&container_type::GetContainer), py::return_value_policy::reference)
        .def(py::self +  py::self)
        .def(py::self += py::self)
        .def(py::self +  float())
        .def(py::self += float())
        .def(py::self -  py::self)
        .def(py::self -= py::self)
        .def(py::self -  float())
        .def(py::self -= float())
        .def(py::self *  float())
        .def(py::self *= float())
        .def(py::self /  float())
        .def(py::self /= float())
        .def("__pow__", &container_type::operator^)
        .def("__ipow__", &container_type::operator^=)
        .def("__neg__", [](container_type& rSelf) { return rSelf.operator*(-1.0); })
        ;
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ContainerVariableDataHolderBase, ContainerVariableDataHolderBase::Pointer>(m, "ContainerVariableDataHolderBase")
        .def("CopyDataFrom", &ContainerVariableDataHolderBase::CopyDataFrom, py::arg("origin_container_data"))
        .def("GetModelPart", py::overload_cast<>(&ContainerVariableDataHolderBase::GetModelPart))
        .def("IsCompatibleWithContainerVariableDataHolder", &ContainerVariableDataHolderBase::IsCompatibleWithContainerVariableDataHolder, py::arg("other_container_data"))
        .def("__str__", &ContainerVariableDataHolderBase::Info)
        ;

    AddContainerTypeToPython<HistoricalDataValueContainer>(m, "HistoricalContainerVariableDataHolder");
    AddContainerTypeToPython<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>(m, "NodalContainerVariableDataHolder");
    AddContainerTypeToPython<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>(m, "ConditionContainerVariableDataHolder");
    AddContainerTypeToPython<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>(m, "ElementContainerVariableDataHolder");
    AddContainerTypeToPython<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>(m, "ConditionPropertiesContainerVariableDataHolder");
    AddContainerTypeToPython<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>(m, "ElementPropertiesContainerVariableDataHolder");

    py::class_<ContainerVariableDataHolderUtils>(m, "ContainerVariableDataHolderUtils")
        .def_static("NormInf", &ContainerVariableDataHolderUtils::NormInf, py::arg("container_data"))
        .def_static("EntityMaxNormL2", &ContainerVariableDataHolderUtils::EntityMaxNormL2, py::arg("container_data"))
        .def_static("InnerProduct", &ContainerVariableDataHolderUtils::InnerProduct, py::arg("container_data_1"), py::arg("container_data_2"))
        ;

    py::class_<OptimizationUtils >(m, "OptimizationUtils")
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ConditionsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ElementsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ConditionsContainerType>)
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ElementsContainerType>)
        ;
}

}  // namespace Python.
} // Namespace Kratos

