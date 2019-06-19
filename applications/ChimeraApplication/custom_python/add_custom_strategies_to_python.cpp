//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//                   Rishith Ellath Meethal
//
// ==============================================================================

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

// Convergence criteria
#include "linear_solvers/linear_solver.h"


namespace Kratos
{

namespace Python
{
using namespace pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    //typedef FSStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
    class_< FSStrategyForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >,
                typename FSStrategyForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >::Pointer,
                BaseSolvingStrategyType >
                (m,"FSStrategyForChimera")
                .def(init< ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,bool,double,double,int,int,unsigned int,unsigned int,bool>())
                .def(init< ModelPart&, FractionalStepSettingsForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool >() )
                .def(init< ModelPart&, FractionalStepSettingsForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool, const Kratos::Variable<int>& >() )
                ;
    class_< ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >,
     typename ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >(m,"ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera")
                .def(init<LinearSolverType::Pointer>());
}

}  // namespace Python.

} // Namespace Kratos

