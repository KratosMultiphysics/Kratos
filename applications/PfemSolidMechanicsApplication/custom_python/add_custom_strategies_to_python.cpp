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
#include "custom_strategies/custom_schemes/residual_based_contact_bossak_scheme.hpp"

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
      typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

      //custom types

      //custom scheme types
      typedef ResidualBasedContactBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedContactBossakSchemeType;


      //********************************************************************
      //*************************STRATEGY CLASSES***************************
      //********************************************************************

    

      //********************************************************************
      //*************************BUILDER AND SOLVER*************************
      //********************************************************************
      


      //********************************************************************
      //*************************SHCHEME CLASSES****************************
      //********************************************************************

 
      // Residual Based Bossak Scheme Type
      class_< ResidualBasedContactBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedContactBossakScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedContactBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;



      //********************************************************************
      //*******************CONVERGENCE CRITERIA CLASSES*********************
      //********************************************************************


 

    }

  }  // namespace Python.

} // Namespace Kratos

