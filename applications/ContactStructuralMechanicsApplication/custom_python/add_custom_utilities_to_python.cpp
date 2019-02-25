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

//Utilities
#include "custom_utilities/simple_contact_search.h"
#include "custom_utilities/advanced_contact_search.h"
#include "custom_utilities/process_factory_utility.h"
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/active_set_utilities.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    // Simple contact search
    py::class_<SimpleContactSearch<2, 2>, typename SimpleContactSearch<2, 2>::Pointer>(m, "SimpleContactSearch2D2N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearch<2, 2>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearch<2, 2>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearch<2, 2>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearch<2, 2>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearch<2, 2>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearch<2, 2>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearch<2, 2>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearch<2, 2>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearch<2, 2>::InvertSearch)
    ;
    py::class_<SimpleContactSearch<3, 3>, typename SimpleContactSearch<3, 3>::Pointer>(m, "SimpleContactSearch3D3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearch<3, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearch<3, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearch<3, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearch<3, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearch<3, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearch<3, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearch<3, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearch<3, 3>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearch<3, 3>::InvertSearch)
    ;
    py::class_<SimpleContactSearch<3, 4>, typename SimpleContactSearch<3, 4>::Pointer>(m, "SimpleContactSearch3D4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearch<3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearch<3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearch<3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearch<3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearch<3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearch<3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearch<3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearch<3, 4>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearch<3, 4>::InvertSearch)
    ;
    py::class_<SimpleContactSearch<3, 3, 4>, typename SimpleContactSearch<3, 3, 4>::Pointer>(m, "SimpleContactSearch3D3N4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearch<3, 3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearch<3, 3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearch<3, 3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearch<3, 3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearch<3, 3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearch<3, 3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearch<3, 3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearch<3, 3, 4>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearch<3, 3, 4>::InvertSearch)
    ;
    py::class_<SimpleContactSearch<3, 4, 3>, typename SimpleContactSearch<3, 4, 3>::Pointer>(m, "SimpleContactSearch3D4N3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearch<3, 4, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearch<3, 4, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearch<3, 4, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearch<3, 4, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearch<3, 4, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearch<3, 4, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearch<3, 4, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearch<3, 4, 3>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearch<3, 4, 3>::InvertSearch)
    ;

    // Advanced contact search
    py::class_<AdvancedContactSearch<2, 2>, typename AdvancedContactSearch<2, 2>::Pointer>(m, "AdvancedContactSearch2D2N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearch<2, 2>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearch<2, 2>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearch<2, 2>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearch<2, 2>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearch<2, 2>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearch<2, 2>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearch<2, 2>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearch<2, 2>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearch<2, 2>::InvertSearch)
    ;
    py::class_<AdvancedContactSearch<3, 3>, typename AdvancedContactSearch<3, 3>::Pointer>(m, "AdvancedContactSearch3D3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearch<3, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearch<3, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearch<3, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearch<3, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearch<3, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearch<3, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearch<3, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearch<3, 3>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearch<3, 3>::InvertSearch)
    ;
    py::class_<AdvancedContactSearch<3, 4>, typename AdvancedContactSearch<3, 4>::Pointer>(m, "AdvancedContactSearch3D4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearch<3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearch<3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearch<3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearch<3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearch<3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearch<3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearch<3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearch<3, 4>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearch<3, 4>::InvertSearch)
    ;
    py::class_<AdvancedContactSearch<3, 3, 4>, typename AdvancedContactSearch<3, 3, 4>::Pointer>(m, "AdvancedContactSearch3D3N4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearch<3, 3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearch<3, 3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearch<3, 3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearch<3, 3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearch<3, 3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearch<3, 3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearch<3, 3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearch<3, 3, 4>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearch<3, 3, 4>::InvertSearch)
    ;
    py::class_<AdvancedContactSearch<3, 4, 3>, typename AdvancedContactSearch<3, 4, 3>::Pointer>(m, "AdvancedContactSearch3D4N3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearch<3, 4, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearch<3, 4, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearch<3, 4, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearch<3, 4, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearch<3, 4, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearch<3, 4, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearch<3, 4, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearch<3, 4, 3>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearch<3, 4, 3>::InvertSearch)
    ;

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
    py::class_<ActiveSetUtilities, typename ActiveSetUtilities::Pointer>(m, "ActiveSetUtilities")
    .def(py::init<>())
    .def("ComputePenaltyFrictionlessActiveSet",&ActiveSetUtilities::ComputePenaltyFrictionlessActiveSet)
    .def("ComputePenaltyFrictionalActiveSet",&ActiveSetUtilities::ComputePenaltyFrictionalActiveSet)
    ;
}

}  // namespace Python.

} // Namespace Kratos

