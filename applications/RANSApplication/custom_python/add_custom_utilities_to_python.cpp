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
#include "custom_utilities/rans_variable_difference_norm_calculation_utility.h"
#include "custom_utilities/test_utilities.h"

namespace Kratos
{
namespace Python
{
void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using RansScalarVariableDifferenceNormCalculationUtilityType = RansVariableDifferenceNormsCalculationUtility<double>;
    py::class_<RansScalarVariableDifferenceNormCalculationUtilityType, RansScalarVariableDifferenceNormCalculationUtilityType::Pointer>(m, "ScalarVariableDifferenceNormCalculationUtility")
        .def(py::init<const ModelPart&, const Variable<double>&>())
        .def("InitializeCalculation", &RansScalarVariableDifferenceNormCalculationUtilityType::InitializeCalculation)
        .def("CalculateDifferenceNorm", &RansScalarVariableDifferenceNormCalculationUtilityType::CalculateDifferenceNorm);

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
        .def("CalculateLogarithmicYPlusLimit", &RansCalculationUtilities::CalculateLogarithmicYPlusLimit, py::arg("kappa"), py::arg("beta"), py::arg("max_iterations") = 20, py::arg("tolerance") = 1e-6);

    m.def_submodule("RansTestUtilities")
        .def("RandomFillNodalHistoricalVariable", &RansApplicationTestUtilities::RandomFillNodalHistoricalVariable<double>)
        .def("RandomFillNodalHistoricalVariable", &RansApplicationTestUtilities::RandomFillNodalHistoricalVariable<array_1d<double, 3>>)
        .def("RandomFillNodalNonHistoricalVariable", &RansApplicationTestUtilities::RandomFillContainerVariable<ModelPart::NodesContainerType, double>)
        .def("RandomFillNodalNonHistoricalVariable", &RansApplicationTestUtilities::RandomFillContainerVariable<ModelPart::NodesContainerType, array_1d<double, 3>>)
        .def("RandomFillConditionVariable", &RansApplicationTestUtilities::RandomFillContainerVariable<ModelPart::ConditionsContainerType, double>)
        .def("RandomFillConditionVariable", &RansApplicationTestUtilities::RandomFillContainerVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def("RandomFillElementVariable", &RansApplicationTestUtilities::RandomFillContainerVariable<ModelPart::ElementsContainerType, double>)
        .def("RandomFillElementVariable", &RansApplicationTestUtilities::RandomFillContainerVariable<ModelPart::ElementsContainerType, array_1d<double, 3>>)
        ;
}

} // namespace Python.
} // Namespace Kratos
