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
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/solver_utilities/gradient_projection_solver_utils.h"

// Include base h
#include "add_custom_utilities_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<OptimizationUtils >(m, "OptimizationUtils")
        .def_static("GetContainerIds", [](const ModelPart::NodesContainerType& rNodes) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationUtils::GetContainerIds(rNodes, values); return values;})
        .def_static("GetContainerIds", [](const ModelPart::ConditionsContainerType& rConditions) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationUtils::GetContainerIds(rConditions, values); return values;})
        .def_static("GetContainerIds", [](const ModelPart::ElementsContainerType& rElements) -> std::vector<IndexType> { std::vector<IndexType> values; OptimizationUtils::GetContainerIds(rElements, values); return values;})
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::NodesContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ConditionsContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ElementsContainerType, double>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("GetContainerVariableToVector", &OptimizationUtils::GetContainerVariableToVector<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationUtils::GetContainerPropertiesVariableToVector<ModelPart::ConditionsContainerType>)
        .def_static("GetContainerPropertiesVariableToVector", &OptimizationUtils::GetContainerPropertiesVariableToVector<ModelPart::ElementsContainerType>)
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
        ;

    py::class_<GradientProjectionSolverUtils >(m, "GradientProjectionSolverUtils")
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection<ModelPart::NodesContainerType, double>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ConditionsContainerType, double>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ElementsContainerType, double>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("CalculateProjectedSearchDirectionAndCorrection", &GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtils::CalculateControlChange<ModelPart::NodesContainerType, double>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtils::CalculateControlChange<ModelPart::ConditionsContainerType, double>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtils::CalculateControlChange<ModelPart::ElementsContainerType, double>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtils::CalculateControlChange<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtils::CalculateControlChange<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("CalculateControlChange", &GradientProjectionSolverUtils::CalculateControlChange<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        ;
}

}  // namespace Python.
} // Namespace Kratos

