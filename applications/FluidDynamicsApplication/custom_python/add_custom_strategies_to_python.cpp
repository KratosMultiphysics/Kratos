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
#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif // KRATOS_USE_AMATRIX

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "processes/process.h"
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
using namespace pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************

    class_< ResidualBasedBlockBuilderAndSolverPeriodic< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            typename ResidualBasedBlockBuilderAndSolverPeriodic< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >
            (m,"ResidualBasedBlockBuilderAndSolverPeriodic")
            .def(init<LinearSolverType::Pointer, const Variable<int>& >());

    class_< FSStrategy< SparseSpaceType,LocalSpaceType, LinearSolverType >,
            typename FSStrategy< SparseSpaceType,LocalSpaceType, LinearSolverType >::Pointer,
            BaseSolvingStrategyType >
            (m,"FSStrategy")
            .def(init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,bool,double,double,int,int,unsigned int,unsigned int,bool>())
            .def(init< ModelPart&, SolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool >() )
            .def(init< ModelPart&, SolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool, const Kratos::Variable<int>& >() )
            .def("CalculateReactions",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateReactions)
            .def("AddIterationStep",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::AddIterationStep)
            .def("ClearExtraIterationSteps",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::ClearExtraIterationSteps)
            ;


    class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent< SparseSpaceType, LocalSpaceType >,
            typename ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent< SparseSpaceType, LocalSpaceType >::Pointer,
            BaseSchemeType  >
            (m,"ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent")
            .def(init<double,double,unsigned int,Process::Pointer >() )
            .def(init<double,double,unsigned int >())// constructor without a turbulence model
            .def(init<double,unsigned int,const Kratos::Variable<int>&>())// constructor without a turbulence model for periodic boundary conditions
            ;

    typedef ResidualBasedSimpleSteadyScheme< SparseSpaceType, LocalSpaceType > ResidualBasedSimpleSteadySchemeType;
    class_< ResidualBasedSimpleSteadySchemeType,
            typename ResidualBasedSimpleSteadySchemeType::Pointer,
            BaseSchemeType >
        (m,"ResidualBasedSimpleSteadyScheme")
        .def(init<double,double,unsigned int,Process::Pointer >() )
        .def(init<double,double,unsigned int >())// constructor without a turbulence model
        .def("GetVelocityRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::GetVelocityRelaxationFactor)
        .def("SetVelocityRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::SetVelocityRelaxationFactor)
        .def("GetPressureRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::GetPressureRelaxationFactor)
        .def("SetPressureRelaxationFactor",&ResidualBasedSimpleSteadySchemeType::SetPressureRelaxationFactor)
        ;


    class_< ResidualBasedPredictorCorrectorBDFSchemeTurbulent< SparseSpaceType, LocalSpaceType >,
            typename ResidualBasedPredictorCorrectorBDFSchemeTurbulent< SparseSpaceType, LocalSpaceType >::Pointer,
            BaseSchemeType >
            (m,"ResidualBasedPredictorCorrectorBDFSchemeTurbulent")
            .def(init<unsigned int,Process::Pointer >() )
            .def(init<unsigned int >())// constructor without a turbulence model
            .def(init<unsigned int,Kratos::Variable<double>&>())// constructor with a non-default flag for slip conditions
            ;

    class_<ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction<SparseSpaceType, LocalSpaceType>,
            typename ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction<SparseSpaceType, LocalSpaceType>::Pointer,
		 ResidualBasedPredictorCorrectorBDFSchemeTurbulent<SparseSpaceType, LocalSpaceType>  >
		(m,"ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction")
                .def(init<unsigned int, Process::Pointer >())
		.def(init<unsigned int >())// constructor without a turbulence model
		.def(init<unsigned int, Kratos::Variable<double>&>())// constructor with a non-default flag for slip conditions
		;

    class_< GearScheme< SparseSpaceType, LocalSpaceType >,
            typename GearScheme< SparseSpaceType, LocalSpaceType >::Pointer,
            BaseSchemeType >
            (m,"GearScheme").
            def(init<>()) // default constructor
            .def(init< Process::Pointer >()) // constructor passing a turbulence model
            ;

	// Convergence criteria
    class_< VelPrCriteria< SparseSpaceType, LocalSpaceType >,
            typename VelPrCriteria< SparseSpaceType, LocalSpaceType >::Pointer,
            ConvergenceCriteria< SparseSpaceType, LocalSpaceType  >>
            (m,"VelPrCriteria")
            .def(init< double, double, double, double>())
            .def("SetEchoLevel",&VelPrCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel)
            ;
}

}  // namespace Python.

} // Namespace Kratos

