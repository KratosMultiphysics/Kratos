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
#include "custom_strategies/explicit_strategy.hpp" 

//builders and solvers
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"
#include "custom_strategies/custom_builders_and_solvers/component_wise_builder_and_solver.hpp"

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
#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"
#include "custom_strategies/custom_schemes/explicit_central_differences_scheme.hpp" 

// modified schemes for the new custom operations on nodal variables
#include "custom_strategies/custom_schemes/residual_based_static_scheme_v2.hpp"
#include "custom_strategies/custom_schemes/residual_based_newmark_scheme_v2.hpp"
#include "custom_strategies/custom_schemes/residual_based_bossak_scheme_v2.hpp"


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
    typedef ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitStrategyType;

    //custom builder_and_solver types
    typedef ResidualBasedBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBuilderAndSolverType;
    typedef ComponentWiseBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ComponentWiseBuilderAndSolverType;

    //custom scheme types
    typedef ResidualBasedStaticScheme< SparseSpaceType, LocalSpaceType > ResidualBasedStaticSchemeType;
    typedef ResidualBasedNewmarkScheme< SparseSpaceType, LocalSpaceType > ResidualBasedNewmarkSchemeType;
    typedef ResidualBasedBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakSchemeType;
    typedef ResidualBasedContactBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedContactBossakSchemeType;    
    typedef ComponentWiseBossakScheme< SparseSpaceType, LocalSpaceType >  ComponentWiseBossakSchemeType;     
    typedef ResidualBasedRelaxationScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRelaxationSchemeType;
    typedef ExplicitCentralDifferencesScheme< SparseSpaceType, LocalSpaceType >  ExplicitCentralDifferencesSchemeType;

    //custom scheme types - modified
    typedef ResidualBasedStaticScheme_V2< SparseSpaceType, LocalSpaceType > ResidualBasedStaticSchemeType_V2;
    typedef ResidualBasedNewmarkScheme_V2< SparseSpaceType, LocalSpaceType > ResidualBasedNewmarkSchemeType_V2;
    typedef ResidualBasedBossakScheme_V2< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakSchemeType_V2;

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


    // Residual Based Relaxation Scheme Type
    class_< ResidualBasedRelaxationSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedRelaxationScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedRelaxationScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

    // Explicit scheme: Central differences 
    class_< ExplicitCentralDifferencesSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ExplicitCentralDifferencesScheme", init< const double, const double, const double, const bool >() )

            .def("Initialize", &ExplicitCentralDifferencesScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

  //********************************************************************
    //*************************SCHEME CLASSES*****************************
  //****MODIFIED TO TEST THE NEW CUSTOM OPERATION ON NODAL VARIABLES****
    //********************************************************************

  // Static Scheme Type
    class_< ResidualBasedStaticSchemeType_V2,
      bases< BaseSchemeType >, boost::noncopyable >
            (
           "ResidualBasedStaticScheme_V2", init< >() )
      
            .def("Initialize", &ResidualBasedStaticSchemeType_V2::Initialize)
            ;

    // Residual Based Newmark Scheme Type
    class_< ResidualBasedNewmarkSchemeType_V2,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedNewmarkScheme_V2", init< double >() )

            .def("Initialize", &ResidualBasedNewmarkSchemeType_V2::Initialize)

            ;

    // Residual Based Bossak Scheme Type
    class_< ResidualBasedBossakSchemeType_V2,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedBossakScheme_V2", init< double , double >() )

            .def("Initialize", &ResidualBasedBossakSchemeType_V2::Initialize)
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

    // Solid Mechanics Explicit Strategy
    class_< ExplicitStrategyType, 
	    bases< BaseSolvingStrategyType >, boost::noncopyable >
      (
       "ExplicitStrategy",
       init < ModelPart&, BaseSchemeType::Pointer,  LinearSolverType::Pointer, bool, bool, bool >())
      
         .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,  bool, bool, bool >())
      .def("SetInitializePerformedFlag", &ExplicitStrategyType::SetInitializePerformedFlag)
      .def("GetInitializePerformedFlag", &ExplicitStrategyType::GetInitializePerformedFlag)
      ;

    
    // Residual Based Newton-Raphson Strategy
    class_< ResidualBasedNewtonRaphsonStrategyType, 
	    bases< BaseSolvingStrategyType >, boost::noncopyable >
      (
       "ResidualBasedNewtonRaphsonStrategy",
       init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
      
      .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
      .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
      .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
      .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
      .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
      .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
      .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
      .def("SetFinalizeSolutionStepFlag", &ResidualBasedNewtonRaphsonStrategyType::SetFinalizeSolutionStepFlag)
      .def("GetFinalizeSolutionStepFlag", &ResidualBasedNewtonRaphsonStrategyType::GetFinalizeSolutionStepFlag)
      ;
    
    // Component Wise Newton-Raphson Strategy
    class_< ComponentWiseNewtonRaphsonStrategyType, 
	    bases< BaseSolvingStrategyType >, boost::noncopyable >
      (
       "ComponentWiseNewtonRaphsonStrategy",
       init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
      
      .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
      .def("SetMaxIterationNumber", &ComponentWiseNewtonRaphsonStrategyType::SetMaxIterationNumber)
      .def("GetMaxIterationNumber", &ComponentWiseNewtonRaphsonStrategyType::GetMaxIterationNumber)
      .def("SetInitializePerformedFlag", &ComponentWiseNewtonRaphsonStrategyType::SetInitializePerformedFlag)
      .def("GetInitializePerformedFlag", &ComponentWiseNewtonRaphsonStrategyType::GetInitializePerformedFlag)
      .def("SetKeepSystemConstantDuringIterations", &ComponentWiseNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
      .def("GetKeepSystemConstantDuringIterations", &ComponentWiseNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
      .def("SetFinalizeSolutionStepFlag", &ComponentWiseNewtonRaphsonStrategyType::SetFinalizeSolutionStepFlag)
      .def("GetFinalizeSolutionStepFlag", &ComponentWiseNewtonRaphsonStrategyType::GetFinalizeSolutionStepFlag)

      ;
  
    // Residual Based Newton-Raphson Line Search Strategy
    class_< ResidualBasedNewtonRaphsonLineSearchStrategyType, 
      bases< BaseSolvingStrategyType >, boost::noncopyable >
      (
       "ResidualBasedNewtonRaphsonLineSearchStrategy",
       init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
      
      .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
      .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetMaxIterationNumber)
      .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetMaxIterationNumber)
      .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetInitializePerformedFlag)
      .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetInitializePerformedFlag)
      .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetKeepSystemConstantDuringIterations)
      .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetKeepSystemConstantDuringIterations)
      .def("SetFinalizeSolutionStepFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::SetFinalizeSolutionStepFlag)
      .def("GetFinalizeSolutionStepFlag", &ResidualBasedNewtonRaphsonLineSearchStrategyType::GetFinalizeSolutionStepFlag)
      ;
     


}

}  // namespace Python.

} // Namespace Kratos

