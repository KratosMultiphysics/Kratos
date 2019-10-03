// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

// Utilities
#include "custom_python/process_factory_utility.h"
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/active_set_utilities.h"
#include "custom_utilities/interface_preprocess.h"
#include "custom_utilities/self_contact_utilities.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    // Process Factory utility
    py::class_<ProcessFactoryUtility, typename ProcessFactoryUtility::Pointer>(m, "ProcessFactoryUtility")
    .def(py::init<>())
    .def(py::init<py::list&>())
    .def(py::init<py::object&>())
    .def("AddProcess",&ProcessFactoryUtility::AddProcess)
    .def("AddProcesses",&ProcessFactoryUtility::AddProcesses)
    .def("ExecuteMethod",&ProcessFactoryUtility::ExecuteMethod)
    .def("ExecuteInitialize",&ProcessFactoryUtility::ExecuteInitialize)
    .def("ExecuteBeforeSolutionLoop",&ProcessFactoryUtility::ExecuteBeforeSolutionLoop)
    .def("ExecuteInitializeSolutionStep",&ProcessFactoryUtility::ExecuteInitializeSolutionStep)
    .def("ExecuteFinalizeSolutionStep",&ProcessFactoryUtility::ExecuteFinalizeSolutionStep)
    .def("ExecuteBeforeOutputStep",&ProcessFactoryUtility::ExecuteBeforeOutputStep)
    .def("ExecuteAfterOutputStep",&ProcessFactoryUtility::ExecuteAfterOutputStep)
    .def("ExecuteFinalize",&ProcessFactoryUtility::ExecuteFinalize)
    .def("IsOutputStep",&ProcessFactoryUtility::IsOutputStep)
    .def("PrintOutput",&ProcessFactoryUtility::PrintOutput)
    .def("Clear",&ProcessFactoryUtility::Clear)
    ;

    // Contact utilities
    py::class_<ContactUtilities, typename ContactUtilities::Pointer>(m, "ContactUtilities")
    .def(py::init<>())
    .def("CalculateRelativeSizeMesh",&ContactUtilities::CalculateRelativeSizeMesh)
    .def("CalculateMaxNodalH",&ContactUtilities::CalculateMaxNodalH)
    .def("CalculateMeanNodalH",&ContactUtilities::CalculateMeanNodalH)
    .def("CalculateMinimalNodalH",&ContactUtilities::CalculateMinimalNodalH)
    .def("CheckActivity",&ContactUtilities::CheckActivity)
    .def("ComputeExplicitContributionConditions",&ContactUtilities::ComputeExplicitContributionConditions)
    .def("ActivateConditionWithActiveNodes",&ContactUtilities::ActivateConditionWithActiveNodes)
    ;

    // Active set utilities
    auto active_set_utilities = m.def_submodule("ActiveSetUtilities");
    active_set_utilities.def("ComputePenaltyFrictionlessActiveSet",&ActiveSetUtilities::ComputePenaltyFrictionlessActiveSet);
    active_set_utilities.def("ComputePenaltyFrictionalActiveSet",&ActiveSetUtilities::ComputePenaltyFrictionalActiveSet);

    // Self-contact utilities
    auto self_contact_utilities = m.def_submodule("SelfContactUtilities");
    self_contact_utilities.def("ComputeSelfContactPairing",&SelfContactUtilities::ComputeSelfContactPairing);
    self_contact_utilities.def("FullAssignmentOfPairs",&SelfContactUtilities::FullAssignmentOfPairs);
    self_contact_utilities.def("NotPredefinedMasterSlave",&SelfContactUtilities::NotPredefinedMasterSlave);

    // Interface preprocess
    py::class_<InterfacePreprocessCondition, typename InterfacePreprocessCondition::Pointer>(m, "InterfacePreprocessCondition")
    .def(py::init<ModelPart&>())
    .def("GenerateInterfacePart",&InterfacePreprocessCondition::GenerateInterfacePart)
    ;
}

}  // namespace Python.

} // Namespace Kratos

