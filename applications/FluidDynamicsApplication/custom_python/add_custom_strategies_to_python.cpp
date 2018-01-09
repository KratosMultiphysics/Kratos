//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_utilities/solver_settings.h"

#include "spaces/ublas_space.h"

// builder_and_solvers
#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_periodic.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/fs_strategy.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_strategies/strategies/residualbased_simple_steady_scheme.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bdf_scheme_turbulent_no_reaction.h"
#include "custom_strategies/strategies/gear_scheme.h"

// convergence criteria
#include "custom_strategies/convergence_criteria/vel_pr_criteria.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************

    class_< ResidualBasedBlockBuilderAndSolverPeriodic< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >,
            boost::noncopyable >
            ("ResidualBasedBlockBuilderAndSolverPeriodic", init<LinearSolverType::Pointer, const Variable<int>& >());

    class_< FSStrategy< SparseSpaceType,LocalSpaceType, LinearSolverType >, bases<BaseSolvingStrategyType>, boost::noncopyable >
            ("FSStrategy",init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,bool,double,double,int,int,unsigned int,unsigned int,bool>())
            .def(init< ModelPart&, SolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool >() )
            .def(init< ModelPart&, SolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool, const Kratos::Variable<int>& >() )
            .def("CalculateReactions",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateReactions)
            .def("AddIterationStep",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::AddIterationStep)
            .def("ClearExtraIterationSteps",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::ClearExtraIterationSteps)
            ;


    class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent< SparseSpaceType, LocalSpaceType >,
            bases< BaseSchemeType >,  boost::noncopyable >
            ("ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent",init<double,double,unsigned int,Process::Pointer >() )
            .def(init<double,double,unsigned int >())// constructor without a turbulence model
            .def(init<double,double,unsigned int,const Kratos::Variable<int>&>())// constructor without a turbulence model for periodic boundary conditions
            .def(init<double,double,unsigned int,Kratos::Variable<double>&>())// constructor with a non-default flag for slip conditions
            ;

    typedef ResidualBasedSimpleSteadyScheme< SparseSpaceType, LocalSpaceType > ResidualBasedSimpleSteadySchemeType;
    class_< ResidualBasedSimpleSteadySchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
        ("ResidualBasedSimpleSteadyScheme",init<double,double,unsigned int,Process::Pointer >() )
        .def(init<double,double,unsigned int >())// constructor without a turbulence model
        .def("GetVelocityRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::GetVelocityRelaxationFactor)
        .def("SetVelocityRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::SetVelocityRelaxationFactor)
        .def("GetPressureRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::GetPressureRelaxationFactor)
        .def("SetPressureRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::SetPressureRelaxationFactor)
        ;


    class_< ResidualBasedPredictorCorrectorBDFSchemeTurbulent< SparseSpaceType, LocalSpaceType >,
            bases< BaseSchemeType >,  boost::noncopyable >
            ("ResidualBasedPredictorCorrectorBDFSchemeTurbulent",init<unsigned int,Process::Pointer >() )
            .def(init<unsigned int >())// constructor without a turbulence model
            .def(init<unsigned int,Kratos::Variable<double>&>())// constructor with a non-default flag for slip conditions
            ;

	class_<ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction<SparseSpaceType, LocalSpaceType>,
		bases< ResidualBasedPredictorCorrectorBDFSchemeTurbulent<SparseSpaceType, LocalSpaceType> >, boost::noncopyable >
		("ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction", init<unsigned int, Process::Pointer >())
		.def(init<unsigned int >())// constructor without a turbulence model
		.def(init<unsigned int, Kratos::Variable<double>&>())// constructor with a non-default flag for slip conditions
		;
            
    class_< GearScheme< SparseSpaceType, LocalSpaceType >,
            bases< BaseSchemeType >,  boost::noncopyable >
            ("GearScheme",init<>()) // default constructor
            .def(init< Process::Pointer >()) // constructor passing a turbulence model
            ;

	// Convergence criteria
    class_< VelPrCriteria< SparseSpaceType, LocalSpaceType >,
            bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
            boost::noncopyable >
            ("VelPrCriteria", init< double, double, double, double>())
            .def("SetEchoLevel",&VelPrCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
            ;
}

}  // namespace Python.

} // Namespace Kratos

