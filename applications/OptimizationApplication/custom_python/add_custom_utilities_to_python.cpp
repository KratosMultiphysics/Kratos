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

// Project includes
#include "includes/model_part.h"
#include "includes/data_communicator.h"

// Application includes
#include "custom_utilities/optimization_variable_utils.h"
#include "custom_utilities/solver_utilities/gradient_projection_solver_utilities.h"

// Include base h
#include "add_custom_utilities_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<OptimizationVariableUtils >(m, "OptimizationVariableUtils")
        .def_static("GetContainerIds", [](const ModelPart::NodesContainerType& rNodes) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationVariableUtils::GetContainerIds(rNodes, values); return values;})
        .def_static("GetContainerIds", [](const ModelPart::ConditionsContainerType& rConditions) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationVariableUtils::GetContainerIds(rConditions, values); return values;})
        .def_static("GetContainerIds", [](const ModelPart::ElementsContainerType& rElements) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationVariableUtils::GetContainerIds(rElements, values); return values;})
        .def_static("GetContainerVariableToVector", &OptimizationVariableUtils::GetContainerVariableToVector<ModelPart::NodesContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationVariableUtils::GetContainerVariableToVector<ModelPart::ConditionsContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationVariableUtils::GetContainerVariableToVector<ModelPart::ElementsContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationVariableUtils::GetContainerVariableToVector<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("GetContainerVariableToVector", &OptimizationVariableUtils::GetContainerVariableToVector<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("GetContainerVariableToVector", &OptimizationVariableUtils::GetContainerVariableToVector<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationVariableUtils::GetContainerPropertiesVariableToVector<ModelPart::ConditionsContainerType>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationVariableUtils::GetContainerPropertiesVariableToVector<ModelPart::ElementsContainerType>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ConditionsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationVariableUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ElementsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationVariableUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        ;

    py::class_<GradientProjectionSolverUtilities >(m, "GradientProjectionSolverUtilities")
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtilities::CalculateProjectedSearchDirectionAndCorrection<ModelPart::NodesContainerType, double>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtilities::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ConditionsContainerType, double>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtilities::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ElementsContainerType, double>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtilities::CalculateProjectedSearchDirectionAndCorrection<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtilities::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtilities::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtilities::CalculateControlChange<ModelPart::NodesContainerType, double>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtilities::CalculateControlChange<ModelPart::ConditionsContainerType, double>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtilities::CalculateControlChange<ModelPart::ElementsContainerType, double>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtilities::CalculateControlChange<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtilities::CalculateControlChange<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtilities::CalculateControlChange<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        ;
}

}  // namespace Python.
} // Namespace Kratos

