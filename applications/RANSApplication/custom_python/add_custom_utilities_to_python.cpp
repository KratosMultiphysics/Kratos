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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_adjoint_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "custom_utilities/rans_nut_utility.h"
#include "custom_utilities/test_utilities.h"

namespace Kratos
{
namespace Python
{
void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RansNutUtility, RansNutUtility::Pointer>(m, "RansNutUtility")
        .def(py::init<ModelPart&, const double, const double, const int>())
        .def("Initialize", &RansNutUtility::Initialize)
        .def("InitializeCalculation", &RansNutUtility::InitializeCalculation)
        .def("CheckConvergence", &RansNutUtility::CheckConvergence)
        .def("UpdateTurbulentViscosity", &RansNutUtility::UpdateTurbulentViscosity)
        ;

    m.def_submodule("RansVariableUtilities")
        .def("ClipScalarVariable", &RansVariableUtilities::ClipScalarVariable)
        .def("GetMinimumScalarValue", &RansVariableUtilities::GetMinimumScalarValue)
        .def("GetMaximumScalarValue", &RansVariableUtilities::GetMaximumScalarValue)
        .def("AddAnalysisStep", &RansVariableUtilities::AddAnalysisStep)
        .def("IsAnalysisStepCompleted", &RansVariableUtilities::IsAnalysisStepCompleted)
        .def("AssignConditionVariableValuesToNodes", &RansVariableUtilities::AssignConditionVariableValuesToNodes<double>, py::arg("model_part"), py::arg("variable"), py::arg("flag"), py::arg("flag_value") = true)
        .def("AssignConditionVariableValuesToNodes", &RansVariableUtilities::AssignConditionVariableValuesToNodes<array_1d<double, 3>>, py::arg("model_part"), py::arg("variable"), py::arg("flag"), py::arg("flag_value") = true)
        .def("AssignBoundaryFlagsToGeometries", &RansVariableUtilities::AssignBoundaryFlagsToGeometries)
        .def("CopyNodalSolutionStepVariablesList", &RansVariableUtilities::CopyNodalSolutionStepVariablesList)
        .def("CalculateTransientVariableConvergence", &RansVariableUtilities::CalculateTransientVariableConvergence<double>)
        .def("CalculateTransientVariableConvergence", &RansVariableUtilities::CalculateTransientVariableConvergence<array_1d<double, 3>>)
        .def("InitializeContainerEntities", &RansVariableUtilities::InitializeContainerEntities<ModelPart::ConditionsContainerType>)
        .def("InitializeContainerEntities", &RansVariableUtilities::InitializeContainerEntities<ModelPart::ElementsContainerType>)
        .def("AssignMaximumVectorComponents", &RansVariableUtilities::AssignMaximumVectorComponents)
        .def("AssignMinimumVectorComponents", &RansVariableUtilities::AssignMinimumVectorComponents)
        .def("CalculateNodalNormal", &RansVariableUtilities::CalculateNodalNormal)
        .def("GetSolutionstepVariableNamesList", &RansVariableUtilities::GetSolutionstepVariableNamesList)
        ;

    m.def_submodule("RansCalculationUtilities")
        .def("CalculateLogarithmicYPlusLimit", &RansCalculationUtilities::CalculateLogarithmicYPlusLimit, py::arg("kappa"), py::arg("beta"), py::arg("max_iterations") = 20, py::arg("tolerance") = 1e-6)
        .def("CalculateWallHeight", &RansCalculationUtilities::CalculateWallHeight)
        .def("CalculateWallDistances", &RansCalculationUtilities::CalculateWallDistances, py::arg("model_part"), py::arg("structure_line_node_ids"))
        ;

    m.def_submodule("RansAdjointUtilities")
        .def("CopyAdjointSolutionToNonHistorical", &RansAdjointUtilities::CopyAdjointSolutionToNonHistorical)
        .def("RescaleAdjointSolution", &RansAdjointUtilities::RescaleAdjointSolution)
        .def("RescaleShapeSensitivity", &RansAdjointUtilities::RescaleShapeSensitivity)
        .def("CalculateTransientReponseFunctionInterpolationError", &RansAdjointUtilities::CalculateTransientReponseFunctionInterpolationError)
        ;
}

} // namespace Python.
} // Namespace Kratos
