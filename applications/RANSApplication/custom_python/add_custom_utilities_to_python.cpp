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

// Project includes

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rans_calculation_utilities.h"
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
        .def(py::init<ModelPart&, const Variable<double>&, const Variable<double>&, const double, const double, const int>())
        .def("Initialize", &RansNutUtility::Initialize)
        .def("InitializeCalculation", &RansNutUtility::InitializeCalculation)
        .def("CheckConvergence", &RansNutUtility::CheckConvergence)
        .def("UpdateTurbulenceData", &RansNutUtility::UpdateTurbulenceData)
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
        .def("SetElementConstitutiveLaws", &RansVariableUtilities::SetElementConstitutiveLaws)
        ;

    m.def_submodule("RansCalculationUtilities")
        .def("CalculateLogarithmicYPlusLimit", &RansCalculationUtilities::CalculateLogarithmicYPlusLimit, py::arg("kappa"), py::arg("beta"), py::arg("max_iterations") = 20, py::arg("tolerance") = 1e-6)
        .def("CalculateWallHeight", &RansVariableUtilities::CalculateWallHeight)
        ;
}

} // namespace Python.
} // Namespace Kratos
