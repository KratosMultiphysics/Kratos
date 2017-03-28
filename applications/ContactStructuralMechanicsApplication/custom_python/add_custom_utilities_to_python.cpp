// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Utilities
#include "custom_utilities/tree_contact_search.h"
#include "custom_utilities/deprecated_tree_contact_search.h"
#include "custom_utilities/exact_mortar_segmentation_utility.h"

namespace Kratos
{
namespace Python
{
void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    
    class_<DeprecatedTreeContactSearch>("DeprecatedTreeContactSearch", init<ModelPart&, ModelPart&, const unsigned int>())
    .def("ModelPartSetter",&DeprecatedTreeContactSearch::ModelPartSetter)
    .def("InitializeMortarConditions",&DeprecatedTreeContactSearch::InitializeMortarConditions)
    .def("TotalClearMortarConditions",&DeprecatedTreeContactSearch::TotalClearMortarConditions)
    .def("PartialClearMortarConditions",&DeprecatedTreeContactSearch::PartialClearMortarConditions)
    .def("CreatePointListMortar",&DeprecatedTreeContactSearch::CreatePointListMortar)
    .def("UpdatePointListMortar",&DeprecatedTreeContactSearch::UpdatePointListMortar)
    .def("CreateMortarConditions",&DeprecatedTreeContactSearch::CreateMortarConditions)
    .def("UpdateMortarConditions",&DeprecatedTreeContactSearch::UpdateMortarConditions)
    .def("CheckMortarConditions",&DeprecatedTreeContactSearch::CheckMortarConditions)
    ;
    
    class_<TreeContactSearch>("TreeContactSearch", init<ModelPart&>())
    .def(init<ModelPart&, const unsigned int>())
    .def(init<ModelPart&, const unsigned int, const double>())
    .def(init<ModelPart&, const unsigned int, const double, const std::string>())
    .def("InitializeMortarConditions",&TreeContactSearch::InitializeMortarConditions)
    .def("TotalClearScalarMortarConditions",&TreeContactSearch::TotalClearScalarMortarConditions)
    .def("TotalClearComponentsMortarConditions",&TreeContactSearch::TotalClearComponentsMortarConditions)
    .def("TotalClearALMFrictionlessMortarConditions",&TreeContactSearch::TotalClearALMFrictionlessMortarConditions)
    .def("PartialClearScalarMortarConditions",&TreeContactSearch::PartialClearScalarMortarConditions)
    .def("PartialClearComponentsMortarConditions",&TreeContactSearch::PartialClearComponentsMortarConditions)
    .def("PartialClearALMFrictionlessMortarConditions",&TreeContactSearch::PartialClearALMFrictionlessMortarConditions)
    .def("CreatePointListMortar",&TreeContactSearch::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch::UpdateMortarConditions)
    .def("CheckMortarConditions",&TreeContactSearch::CheckMortarConditions)
    ;
    
    // Exact integration (for testing)
    class_<ExactMortarIntegrationUtility<2,2>>("ExactMortarIntegrationUtility2D2N", init<const unsigned int>())
    .def(init< >())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactIntegration)
    ;
    class_<ExactMortarIntegrationUtility<3,3>>("ExactMortarIntegrationUtility3D3N", init<const unsigned int>())
    .def(init< >())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactIntegration)
    ;
    class_<ExactMortarIntegrationUtility<3,4>>("ExactMortarIntegrationUtility3D4N", init<const unsigned int>())
    .def(init< >())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactIntegration)
    ;
  
}

}  // namespace Python.

} // Namespace Kratos

