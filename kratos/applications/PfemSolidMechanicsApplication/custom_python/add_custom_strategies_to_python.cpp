//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes 

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/residual_based_newton_raphson_line_search_implex_strategy.hpp"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//schemes
#include "custom_strategies/custom_schemes/residual_based_U_wP_static_scheme.hpp"

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
      typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

      //custom types
      typedef ResidualBasedUwPStaticScheme< SparseSpaceType, LocalSpaceType > ResidualBasedUwPStaticSchemeType;
      typedef ResidualBasedNewtonRaphsonLineSearchImplexStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonLineSearchImplexStrategyType;

      //********************************************************************
      //*************************STRATEGY CLASSES***************************
      //********************************************************************

    

      //********************************************************************
      //*************************BUILDER AND SOLVER*************************
      //********************************************************************
      


      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************
      // Static Scheme Type
      class_< ResidualBasedUwPStaticSchemeType,
        bases< BaseSchemeType >, boost::noncopyable >
              (
             "ResidualBasedUwPStaticScheme", init< >() )
        
              .def("Initialize", &ResidualBasedUwPStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)
              ;




      //********************************************************************
      //*******************CONVERGENCE CRITERIA CLASSES*********************
      //********************************************************************


    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Residual Based Newton-Raphson Line Search Strategy
    class_< ResidualBasedNewtonRaphsonLineSearchImplexStrategyType, 
      bases< BaseSolvingStrategyType >, boost::noncopyable >
      (
       "ResidualBasedNewtonRaphsonLineSearchImplexStrategy",
       init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaBaseType::Pointer, int, bool, bool, bool>())
      
      .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaBaseType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
      .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetMaxIterationNumber)
      .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetMaxIterationNumber)
      .def("SetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetInitializePerformedFlag)
      .def("GetInitializePerformedFlag", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetInitializePerformedFlag)
      .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::SetKeepSystemConstantDuringIterations)
      .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonLineSearchImplexStrategyType::GetKeepSystemConstantDuringIterations)
      ;
     



 

    }

  }  // namespace Python.

} // Namespace Kratos

