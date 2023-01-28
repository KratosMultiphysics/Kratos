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

// Include base h
#include "add_custom_utilities_to_python.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::enum_<ContainerData::ContainerDataType>(m,"ContainerDataType")
        .value("NodalHistorical", ContainerData::ContainerDataType::NodalHistorical)
        .value("NodalNonHistorical", ContainerData::ContainerDataType::NodalNonHistorical)
        .value("ConditionNonHistorical", ContainerData::ContainerDataType::ConditionNonHistorical)
        .value("ConditionProperties", ContainerData::ContainerDataType::ConditionProperties)
        .value("ElementNonHistorical", ContainerData::ContainerDataType::ElementNonHistorical)
        .value("ElementProperties", ContainerData::ContainerDataType::ElementProperties)
        ;

    py::class_<ContainerData, ContainerData::Pointer>(m, "ContainerData")
        .def(py::init<ModelPart&, const ContainerData::ContainerDataType&>())
        .def(py::init<const ContainerData&>())
        .def("AssignDataToContainerVariable", &ContainerData::AssignDataToContainerVariable<double>)
        .def("AssignDataToContainerVariable", &ContainerData::AssignDataToContainerVariable<array_1d<double, 3>>)
        .def("ReadDataFromContainerVariable", &ContainerData::ReadDataFromContainerVariable<double>)
        .def("ReadDataFromContainerVariable", &ContainerData::ReadDataFromContainerVariable<array_1d<double, 3>>)
        .def("SetDataForContainerVariable", &ContainerData::SetDataForContainerVariable<double>)
        .def("SetDataForContainerVariable", &ContainerData::SetDataForContainerVariable<array_1d<double, 3>>)
        .def("Clone", &ContainerData::Clone)
        .def("IsSameContainer", &ContainerData::IsSameContainer)
        .def("NormInf", &ContainerData::NormInf)
        .def("InnerProduct", &ContainerData::InnerProduct)
        .def("GetContainerDataType", &ContainerData::GetContainerDataType)
        .def("GetData", py::overload_cast<>(&ContainerData::GetData))
        .def("GetModelPart", py::overload_cast<>(&ContainerData::GetModelPart))
        .def("GetContainer", py::overload_cast<>(&ContainerData::GetContainer))
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(py::self * float())
        .def(py::self *= float())
        .def(py::self / float())
        .def(py::self /= float())
        .def("__str__", &ContainerData::Info)
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

