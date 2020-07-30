//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//

// System includes


// External includes


// Project includes
#include "custom_python/add_custom_strategies_to_python.h"


#include "spaces/ublas_space.h"

// Strategies
#include "custom_utilities/solver_settings.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/fs_strategy_for_chimera.h"

// Builder and solvers
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{


void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef FSStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef FSStrategyForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType > FSSTrategyForChimeraType;
    //*************************STRATEGY CLASSES***************************
    py::class_< FSSTrategyForChimeraType,
                typename FSSTrategyForChimeraType::Pointer,
                BaseSolvingStrategyType >
                (m,"FSStrategyForChimera")
                .def(py::init< ModelPart&, FractionalStepSettingsForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool >() )
                ;
    //*************************B&S CLASSES***************************
    py::class_< ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >,
     typename ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >(m,"ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera")
                .def(py::init<LinearSolverType::Pointer>());
}

}  // namespace Python.

} // Namespace Kratos

