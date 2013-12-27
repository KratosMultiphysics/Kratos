//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/residual_based_newton_raphson_strategy.hpp"
#include "custom_strategies/component_wise_newton_raphson_strategy.hpp"
#include "custom_strategies/residual_based_newton_raphson_line_search_strategy.hpp"

//builders and solvers
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"
#include "custom_strategies/custom_builders_and_solvers/component_wise_builder_and_solver.hpp"
#include "custom_strategies/custom_builders_and_solvers/block_residual_based_builder_and_solver.hpp"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergence_criteria/displacement_convergence_criterion.hpp"
#include "custom_strategies/custom_convergence_criteria/component_wise_residual_convergence_criterion.hpp"

//schemes
#include "custom_strategies/custom_schemes/residual_based_static_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_newmark_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_bossak_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_contact_bossak_scheme.hpp"
#include "custom_strategies/custom_schemes/component_wise_bossak_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_rotation_bossak_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"

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

    //base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;


    //custom strategy types
    typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;
    typedef ComponentWiseNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ComponentWiseNewtonRaphsonStrategyType;
    typedef ResidualBasedNewtonRaphsonLineSearchStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonLineSearchStrategyType;

    //custom builder_and_solver types
    typedef ResidualBasedBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBuilderAndSolverType;
    typedef ComponentWiseBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ComponentWiseBuilderAndSolverType;
    typedef BlockResidualBasedBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BlockResidualBasedBuilderAndSolverType;

    //custom scheme types
    typedef ResidualBasedStaticScheme< SparseSpaceType, LocalSpaceType > ResidualBasedStaticSchemeType;
    typedef ResidualBasedNewmarkScheme< SparseSpaceType, LocalSpaceType > ResidualBasedNewmarkSchemeType;
    typedef ResidualBasedBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakSchemeType;
    typedef ResidualBasedContactBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedContactBossakSchemeType;    
    typedef ComponentWiseBossakScheme< SparseSpaceType, LocalSpaceType >  ComponentWiseBossakSchemeType;     
    typedef ResidualBasedRotationBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRotationBossakSchemeType;
    typedef ResidualBasedRelaxationScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRelaxationSchemeType;


    //custom convergence criterion types
    typedef DisplacementConvergenceCriterion< SparseSpaceType,  LocalSpaceType > DisplacementConvergenceCriterionType;
    typedef ComponentWiseResidualConvergenceCriterion< SparseSpaceType,  LocalSpaceType > ComponentWiseResidualConvergenceCriterionType;


    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    // Residual Based Builder and Solver
    class_< ResidualBasedBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
            (
	      "ResidualBasedBuilderAndSolver", init< LinearSolverType::Pointer > ()
            );

    // Component Wise Builder and Solver
    class_< ComponentWiseBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
            (
              "ComponentWiseBuilderAndSolver", init< LinearSolverType::Pointer > ()
            );

    // Block Residual Based Builder and Solver   
    class_< BlockResidualBasedBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > 
            (
              "BlockResidualBasedBuilderAndSolver", init< LinearSolverType::Pointer > ()
            );


    //********************************************************************
    //*************************SHCHEME CLASSES****************************
    //********************************************************************

    // Static Scheme Type
    class_< ResidualBasedStaticSchemeType,
	    bases< BaseSchemeType >, boost::noncopyable >
            (
	         "ResidualBasedStaticScheme", init< >() )
      
            .def("Initialize", &ResidualBasedStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

    // Residual Based Newmark Scheme Type
    class_< ResidualBasedNewmarkSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedNewmarkScheme", init< double >() )

            .def("Initialize", &ResidualBasedNewmarkScheme<SparseSpaceType, LocalSpaceType>::Initialize)

            ;

    // Residual Based Bossak Scheme Type
    class_< ResidualBasedBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedBossakScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

    // Residual Based Bossak Scheme Type
    class_< ResidualBasedContactBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedContactBossakScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedContactBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

    // Component Wise Bossak Scheme Type
    class_< ComponentWiseBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ComponentWiseBossakScheme", init< double , double >() )

            .def("Initialize", &ComponentWiseBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;


    // Residual Based Rotational Bossak Scheme Type
    class_< ResidualBasedRotationBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedRotationBossakScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedRotationBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;


    // Residual Based Relaxation Scheme Type
    class_< ResidualBasedRelaxationSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedRelaxationScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedRelaxationScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    // Displacement Convergence Criterion
    class_< DisplacementConvergenceCriterionType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
                "DisplacementConvergenceCriterion", init<double, double >()
            );


    // Component Wise Residual Convergence Criterion
    class_< ComponentWiseResidualConvergenceCriterionType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
                "ComponentWiseResidualConvergenceCriterion", init<double, double >()
            );



    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************


    // Residual Based Newton-Raphson Strategy
    class_< ResidualBasedNewtonRaphsonStrategyType, 
	    bases< BaseSolvingStrategyType >, boost::noncopyable >
            (
	     "ResidualBasedNewtonRaphsonStrategy",
	     init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool
	     >())
      
           .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
           .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
           .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
           .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
           .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
           .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
           .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
      ;

    // Component Wise Newton-Raphson Strategy
    class_< ComponentWiseNewtonRaphsonStrategyType, 
	    bases< BaseSolvingStrategyType >, boost::noncopyable >
            (
	     "ComponentWiseNewtonRaphsonStrategy",
	     init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool
	     >())
      
           .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
           .def("SetMaxIterationNumber", &ComponentWiseNewtonRaphsonStrategyType::SetMaxIterationNumber)
           .def("GetMaxIterationNumber", &ComponentWiseNewtonRaphsonStrategyType::GetMaxIterationNumber)
           .def("SetInitializePerformedFlag", &ComponentWiseNewtonRaphsonStrategyType::SetInitializePerformedFlag)
           .def("GetInitializePerformedFlag", &ComponentWiseNewtonRaphsonStrategyType::GetInitializePerformedFlag)
           .def("SetKeepSystemConstantDuringIterations", &ComponentWiseNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
           .def("GetKeepSystemConstantDuringIterations", &ComponentWiseNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
      ;
  
    // Residual Based Newton-Raphson Line Search Strategy
    class_< ResidualBasedNewtonRaphsonLineSearchStrategyType, 
	    bases< BaseSolvingStrategyType >, boost::noncopyable >
            (
	     "ResidualBasedNewtonRaphsonLineSearchStrategy",
	     init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool
	     >())
      
           .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
           .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetMaxIterationNumber)
           .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetMaxIterationNumber)
           .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetInitializePerformedFlag)
           .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetInitializePerformedFlag)
           .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetKeepSystemConstantDuringIterations)
           .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetKeepSystemConstantDuringIterations)
      ;
	   


}

}  // namespace Python.

} // Namespace Kratos

