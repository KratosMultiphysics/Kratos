// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
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
    
    class_<TreeContactSearch>("TreeContactSearch", init<ModelPart&, ModelPart&, const unsigned int>())
    .def("ModelPartSetter",&TreeContactSearch::ModelPartSetter)
    .def("InitializeNTNConditions",&TreeContactSearch::InitializeNTNConditions)
    .def("InitializeNTSConditions",&TreeContactSearch::InitializeNTSConditions)
    .def("InitializeMortarConditions",&TreeContactSearch::InitializeMortarConditions)
    .def("InitializeALMFrictionlessMortarConditions",&TreeContactSearch::InitializeALMFrictionlessMortarConditions)
    .def("TotalClearNTNConditions",&TreeContactSearch::TotalClearNTNConditions)
    .def("TotalClearNTSConditions",&TreeContactSearch::TotalClearNTSConditions)
    .def("TotalClearMortarConditions",&TreeContactSearch::TotalClearMortarConditions)
    .def("TotalClearALMFrictionlessMortarConditions",&TreeContactSearch::TotalClearALMFrictionlessMortarConditions)
    .def("PartialClearNTNConditions",&TreeContactSearch::PartialClearNTNConditions)
    .def("PartialClearNTSConditions",&TreeContactSearch::PartialClearNTSConditions)
    .def("PartialClearMortarConditions",&TreeContactSearch::PartialClearMortarConditions)
    .def("PartialClearALMFrictionlessMortarConditions",&TreeContactSearch::PartialClearALMFrictionlessMortarConditions)
    .def("CreatePointListNTN",&TreeContactSearch::CreatePointListNTN)
    .def("CreatePointListNTS",&TreeContactSearch::CreatePointListNTS)
    .def("CreatePointListMortar",&TreeContactSearch::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch::UpdatePointListMortar)
    .def("CreateNTNConditions",&TreeContactSearch::CreateNTNConditions)
    .def("CreateNTSConditions",&TreeContactSearch::CreateNTSConditions)
    .def("CreateMortarConditions",&TreeContactSearch::CreateMortarConditions)
    .def("UpdateNTNConditions",&TreeContactSearch::UpdateNTNConditions)
    .def("UpdateNTSConditions",&TreeContactSearch::UpdateNTSConditions)
    .def("UpdateMortarConditions",&TreeContactSearch::UpdateMortarConditions)
    .def("CheckMortarConditions",&TreeContactSearch::CheckMortarConditions)
    ;
    
    // Exact integration (for testing)
    class_<ExactMortarIntegrationUtility<2,2>>("ExactMortarIntegrationUtility2D2N", init<const GeometryType&, const array_1d<double, 3>&, const unsigned int>())
    .def(init<const GeometryType&, const array_1d<double, 3>&>())
    .def(init<const Condition::Pointer>())
    .def(init<const Condition::Pointer, const unsigned int>())
    .def("GetExactIntegration",&ExactMortarIntegrationUtility<2,2>::GetExactIntegration)
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactIntegration)
    ;
    class_<ExactMortarIntegrationUtility<3,3>>("ExactMortarIntegrationUtility3D3N", init<const GeometryType&, const array_1d<double, 3>&, const unsigned int>())
    .def(init<const GeometryType&, const array_1d<double, 3>&>())
    .def(init<const Condition::Pointer>())
    .def(init<const Condition::Pointer, const unsigned int>())
    .def("GetExactIntegration",&ExactMortarIntegrationUtility<3,3>::GetExactIntegration)
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactIntegration)
    ;
    class_<ExactMortarIntegrationUtility<3,4>>("ExactMortarIntegrationUtility3D4N", init<const GeometryType&, const array_1d<double, 3>&, const unsigned int>())
    .def(init<const GeometryType&, const array_1d<double, 3>&>())
    .def(init<const Condition::Pointer>())
    .def(init<const Condition::Pointer, const unsigned int>())
    .def("GetExactIntegration",&ExactMortarIntegrationUtility<3,4>::GetExactIntegration)
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactIntegration)
    ;
  
}

}  // namespace Python.

} // Namespace Kratos

