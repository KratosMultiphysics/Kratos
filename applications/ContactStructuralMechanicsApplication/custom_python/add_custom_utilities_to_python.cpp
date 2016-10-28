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
    .def("ClearNTNConditions",&TreeContactSearch::ClearNTNConditions)
    .def("ClearNTSConditions",&TreeContactSearch::ClearNTSConditions)
    .def("ClearMortarConditions",&TreeContactSearch::ClearMortarConditions)
    .def("CreatePointListNTN",&TreeContactSearch::CreatePointListNTN)
    .def("CreatePointListNTS",&TreeContactSearch::CreatePointListNTS)
    .def("CreatePointListMortar",&TreeContactSearch::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch::UpdatePointListMortar)
    .def("CreateNTNConditions",&TreeContactSearch::CreateNTNConditions)
    .def("CreateNTSConditions",&TreeContactSearch::CreateNTSConditions)
    .def("CreateMortarConditions",&TreeContactSearch::CreateMortarConditions)
    .def("CheckMortarConditions",&TreeContactSearch::CheckMortarConditions)
    ;
  
}

}  // namespace Python.

} // Namespace Kratos

