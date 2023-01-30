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
#include "custom_utilities/solver_utilities/gradient_projection_solver_utils.h"
#include "custom_utilities/container_data.h"
#include "custom_utilities/container_data_utils.h"

// Include base h
#include "add_custom_utilities_to_python.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace Kratos {
namespace Python {

template<class ContainerDataType>
void AddContainerTypeToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_type = ContainerData<ContainerDataType>;
    py::class_<container_type, typename container_type::Pointer, ContainerDataBase>(m, rName.c_str())
        .def(py::init<ModelPart&>(), py::arg("model_part"))
        .def(py::init<const container_type&>(), py::arg("other_container_data_to_copy_from"))
        .def("AssignDataToContainerVariable", &container_type::template AssignDataToContainerVariable<double>, py::arg("scalar_variable"))
        .def("AssignDataToContainerVariable", &container_type::template AssignDataToContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("ReadDataFromContainerVariable", &container_type::template ReadDataFromContainerVariable<double>, py::arg("scalar_variable"))
        .def("ReadDataFromContainerVariable", &container_type::template ReadDataFromContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("SetDataForContainerVariable", &container_type::template SetDataForContainerVariable<double>, py::arg("scalar_variable"), py::arg("scalar_value"))
        .def("SetDataForContainerVariable", &container_type::template SetDataForContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"), py::arg("Array3_value"))
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

    py::class_<ContainerDataBase, ContainerDataBase::Pointer>(m, "ContainerDataBase")
        .def("CopyDataFrom", &ContainerDataBase::CopyDataFrom, py::arg("origin_container_data"))
        .def("GetModelPart", py::overload_cast<>(&ContainerDataBase::GetModelPart))
        .def("IsCompatibleWithContainerData", &ContainerDataBase::IsCompatibleWithContainerData, py::arg("other_container_data"))
        .def("__str__", &ContainerDataBase::Info)
        ;

    AddContainerTypeToPython<HistoricalDataValueContainer>(m, "HistoricalContainerData");
    AddContainerTypeToPython<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>(m, "NodalContainerData");
    AddContainerTypeToPython<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>(m, "ConditionContainerData");
    AddContainerTypeToPython<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>(m, "ElementContainerData");
    AddContainerTypeToPython<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>(m, "ConditionPropertiesContainerData");
    AddContainerTypeToPython<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>(m, "ElementPropertiesContainerData");

    py::class_<ContainerDataUtils >(m, "ContainerDataUtils")
        .def_static("NormInf", &ContainerDataUtils::NormInf)
        .def_static("InnerProduct", &ContainerDataUtils::InnerProduct)
        ;

    py::class_<OptimizationUtils >(m, "OptimizationUtils")
        .def_static("GetContainerIds", [](const ModelPart::NodesContainerType& rNodes) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationUtils::GetContainerIds(values, rNodes); return values;})
        .def_static("GetContainerIds", [](const ModelPart::ConditionsContainerType& rConditions) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationUtils::GetContainerIds(values, rConditions); return values;})
        .def_static("GetContainerIds", [](const ModelPart::ElementsContainerType& rElements) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationUtils::GetContainerIds(values, rElements); return values;})
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::NodesContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ConditionsContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ElementsContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationUtils::GetContainerPropertiesVariableToVector<ModelPart::ConditionsContainerType, double>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationUtils::GetContainerPropertiesVariableToVector<ModelPart::ElementsContainerType, double>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationUtils::GetContainerPropertiesVariableToVector<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationUtils::GetContainerPropertiesVariableToVector<ModelPart::ElementsContainerType, array_1d<double, 3>>)
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
        .def_static("CalculateVectorL2Norm", &OptimizationUtils::CalculateVectorL2Norm)
        .def_static("AssignVectorToContainerPropertiesVariable", &OptimizationUtils::AssignVectorToContainerPropertiesVariable<ModelPart::ConditionsContainerType, double>)
        .def_static("AssignVectorToContainerPropertiesVariable", &OptimizationUtils::AssignVectorToContainerPropertiesVariable<ModelPart::ElementsContainerType, double>)
        .def_static("AssignVectorToContainerPropertiesVariable", &OptimizationUtils::AssignVectorToContainerPropertiesVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("AssignVectorToContainerPropertiesVariable", &OptimizationUtils::AssignVectorToContainerPropertiesVariable<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        .def_static("AssignVectorToContainerVariable", &OptimizationUtils::AssignVectorToContainerVariable<ModelPart::NodesContainerType, double>)
        .def_static("AssignVectorToContainerVariable", &OptimizationUtils::AssignVectorToContainerVariable<ModelPart::ConditionsContainerType, double>)
        .def_static("AssignVectorToContainerVariable", &OptimizationUtils::AssignVectorToContainerVariable<ModelPart::ElementsContainerType, double>)
        .def_static("AssignVectorToContainerVariable", &OptimizationUtils::AssignVectorToContainerVariable<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("AssignVectorToContainerVariable", &OptimizationUtils::AssignVectorToContainerVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("AssignVectorToContainerVariable", &OptimizationUtils::AssignVectorToContainerVariable<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ConditionsContainerType>)
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ElementsContainerType>)
        .def_static("AssignVectorToContainerHistoricalVariable", &OptimizationUtils::AssignVectorToContainerHistoricalVariable<double>)
        .def_static("AssignVectorToContainerHistoricalVariable", &OptimizationUtils::AssignVectorToContainerHistoricalVariable<array_1d<double, 3>>)
        .def_static("GetHistoricalContainerVariableToVector", &OptimizationUtils::GetHistoricalContainerVariableToVector<double>)
        .def_static("GetHistoricalContainerVariableToVector", &OptimizationUtils::GetHistoricalContainerVariableToVector<array_1d<double, 3>>)
        .def_static("AddVectors", &OptimizationUtils::AddVectors)
        .def_static("SubstractVectors", &OptimizationUtils::SubstractVectors)
        .def_static("MultiplyVector", &OptimizationUtils::MultiplyVector)
        .def_static("DivideVector", &OptimizationUtils::DivideVector)
        .def_static("NormInf", &OptimizationUtils::NormInf)
        ;

    py::class_<GradientProjectionSolverUtils >(m, "GradientProjectionSolverUtils")
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection)
        .def_static("CalculateControlUpdate", &GradientProjectionSolverUtils::CalculateControlUpdate)
        ;
}

}  // namespace Python.
} // Namespace Kratos

