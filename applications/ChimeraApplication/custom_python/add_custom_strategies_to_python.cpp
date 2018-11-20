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
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/fs_strategy_for_chimera.h"

//#include "custom_strategies/custom_strategies/residual_based_arc_length_strategy.hpp"
//#include "custom_strategies/custom_strategies/eigensolver_strategy.hpp"

// Schemes
#include "solving_strategies/schemes/scheme.h"
//#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"
//#include "custom_strategies/custom_schemes/eigensolver_dynamic_scheme.hpp"

// Builder and solvers
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
//#include "custom_strategies/custom_builder_and_solver/residualbased_block_builder_and_solver_with_mpc_chimera.h"
#include "custom_strategies/custom_builder_and_solver/residualbased_block_builder_and_solver_with_constraints_for_chimera.h"


// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
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
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    // Custom convergence criterion types

    // Custom builder and solvers types

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
    class_< FSStrategyForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >,
                typename FSStrategyForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >::Pointer,
                BaseSolvingStrategyType >
                (m,"FSStrategyForChimera")
                .def(init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,bool,double,double,int,int,unsigned int,unsigned int,bool>())
                .def(init< ModelPart&, SolverSettingsForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool >() )
                .def(init< ModelPart&, SolverSettingsForChimera< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool, const Kratos::Variable<int>& >() )
                .def("CalculateReactions",&FSStrategyForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateReactions)
                .def("AddIterationStep",&FSStrategyForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::AddIterationStep)
                .def("ClearExtraIterationSteps",&FSStrategyForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::ClearExtraIterationSteps)
                ;

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************
    /* class_< ResidualBasedBlockBuilderAndSolverWithMpcChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                bases< ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >,
                boost::noncopyable >
                ("ResidualBasedBlockBuilderAndSolverWithMpcChimera", init<LinearSolverType::Pointer>()); */

/*     class_< ResidualBasedBlockBuilderAndSolverWithMpcChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >,
     typename ResidualBasedBlockBuilderAndSolverWithMpcChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >(m,"ResidualBasedBlockBuilderAndSolverWithMpcChimera")
                .def(init<LinearSolverType::Pointer>()); */

    class_< ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >,
     typename ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
                ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >(m,"ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera")
                .def(init<LinearSolverType::Pointer>());

    

}

}  // namespace Python.

} // Namespace Kratos

