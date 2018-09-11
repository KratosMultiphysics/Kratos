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
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Utilities
#include "custom_utilities/tree_contact_search.h"
#include "custom_utilities/process_factory_utility.h"
#include "custom_utilities/contact_utilities.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    // Tree contact search
    class_<TreeContactSearch<2, 2>, typename TreeContactSearch<2, 2>::Pointer>(m, "TreeContactSearch2D2N")
    .def(init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<2, 2>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<2, 2>::ClearMortarConditions)
    .def("CheckContactModelParts",&TreeContactSearch<2, 2>::CheckContactModelParts)
    .def("CreatePointListMortar",&TreeContactSearch<2, 2>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<2, 2>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<2, 2>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<2, 2>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<2, 2>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<2, 2>::InvertSearch)
    ;
    class_<TreeContactSearch<3, 3>, typename TreeContactSearch<3, 3>::Pointer>(m, "TreeContactSearch3D3N")
    .def(init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<3, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<3, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&TreeContactSearch<3, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&TreeContactSearch<3, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<3, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<3, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<3, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<3, 3>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<3, 3>::InvertSearch)
    ;
    class_<TreeContactSearch<3, 4>, typename TreeContactSearch<3, 4>::Pointer>(m, "TreeContactSearch3D4N")
    .def(init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&TreeContactSearch<3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&TreeContactSearch<3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<3, 4>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<3, 4>::InvertSearch)
    ;
    class_<TreeContactSearch<3, 3, 4>, typename TreeContactSearch<3, 3, 4>::Pointer>(m, "TreeContactSearch3D3N4N")
    .def(init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<3, 3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<3, 3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&TreeContactSearch<3, 3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&TreeContactSearch<3, 3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<3, 3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<3, 3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<3, 3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<3, 3, 4>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<3, 3, 4>::InvertSearch)
    ;
    class_<TreeContactSearch<3, 4, 3>, typename TreeContactSearch<3, 4, 3>::Pointer>(m, "TreeContactSearch3D4N3N")
    .def(init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<3, 4, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<3, 4, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&TreeContactSearch<3, 4, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&TreeContactSearch<3, 4, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<3, 4, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<3, 4, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<3, 4, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<3, 4, 3>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<3, 4, 3>::InvertSearch)
    ;

    // Process Factory utility
    class_<ProcessFactoryUtility, typename ProcessFactoryUtility::Pointer>(m, "ProcessFactoryUtility")
    .def(init< >())
    .def(init<  list&>())
    .def(init<  object&>())
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
    class_<ContactUtilities, typename ContactUtilities::Pointer>(m, "ContactUtilities")
    .def(init<>())
    .def("CalculateRelativeSizeMesh",&ContactUtilities::CalculateRelativeSizeMesh)
    .def("CalculateMaxNodalH",&ContactUtilities::CalculateMaxNodalH)
    .def("CalculateMeanNodalH",&ContactUtilities::CalculateMeanNodalH)
    .def("CalculateMinimalNodalH",&ContactUtilities::CalculateMinimalNodalH)
    ;
}

}  // namespace Python.

} // Namespace Kratos

