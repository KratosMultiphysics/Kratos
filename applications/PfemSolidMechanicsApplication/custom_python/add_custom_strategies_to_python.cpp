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

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//schemes

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
      typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
      typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;


      //custom types

      //********************************************************************
      //*************************STRATEGY CLASSES***************************
      //********************************************************************

      // class_< TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
      // 	      bases< BaseSolvingStrategyType >,  boost::noncopyable >
      // ("TestStrategy", 
      //  init<ModelPart&, LinearSolverType::Pointer, int, int, bool >() )
      // .def("MoveNodes",&TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
      // ;
    

      //********************************************************************
      //*************************BUILDER AND SOLVER*************************
      //********************************************************************
      
      // Residual Based Builder and Solver
      //typedef ResidualBasedBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBuilderAndSolverType;
      
      //class_< ResidualBasedBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedBuilderAndSolver", init< LinearSolverType::Pointer > ());

     
      // Block Residual Based Builder and Solver
      //typedef BlockResidualBasedBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BlockResidualBasedBuilderAndSolverType;
      //class_< BlockResidualBasedBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("BlockResidualBasedBuilderAndSolver", init< LinearSolverType::Pointer > ());


      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************

 
      // class_< ResidualBasedBossakSchemeType,
      // 	      bases< BaseSchemeType >,  boost::noncopyable >
      // 	(
      // 	 "ResidualBasedBossakScheme", init< double , double >() )

      // 	 .def("Initialize", &ResidualBasedBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      // 	;
    


      

      //********************************************************************
      //*******************CONVERGENCE CRITERIA CLASSES*********************
      //********************************************************************


     //  class_< DisplacementConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ,
     // 	      bases< ConvergenceCriteriaBaseType >, boost::noncopyable >
     // 	(
     // 	 "DisplacementConvergenceCriteria", init<double, double >()
     // 	 );


 

    }

  }  // namespace Python.

} // Namespace Kratos

