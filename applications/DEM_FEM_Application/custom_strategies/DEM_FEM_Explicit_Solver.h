//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Feng Chun $
//   Date:                $Date: 2013-9-16$
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_DEM_FEM_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_DEM_FEM_EXPLICIT_SOLVER_STRATEGY


// Project includes
#include "utilities/timer.h"


#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#define CUSTOMTIMER 1  // ACTIVATES AND DISABLES ::TIMER:::::

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "dem_fem__application.h"



#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "spatial_containers/spatial_search.h"



namespace Kratos
{
  
  template<
  class TSparseSpace,
  class TDenseSpace,
  class TLinearSolver>
  class DemFemExplicitSolverStrategy: public  SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
  {
      public:
      
	  KRATOS_CLASS_POINTER_DEFINITION(DemFemExplicitSolverStrategy);
        
      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>   BaseType;
	  
	  typedef typename BaseType::TDataType TDataType;	  
	  typedef TSparseSpace SparseSpaceType;
	  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;	  
	  typedef typename BaseType::TSchemeType TSchemeType;
	  
	  typedef ModelPart::NodesContainerType NodesArrayType;
	  typedef ModelPart::ElementsContainerType ElementsArrayType;
	  typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	  typedef ModelPart::ConditionsContainerType::ContainerType ConditionsContainerType;
	  
	   typedef Element::Pointer ElementPointer;
	   typedef typename std::vector<ElementPointer> ElementPointerVector;
	   typedef typename std::vector<ElementPointer>::iterator ElementPointerIterator;
	   
      typedef SpatialSearch::ResultElementsContainerType                ResultElementsContainerType;
      typedef SpatialSearch::VectorResultElementsContainerType          VectorResultElementsContainerType;
	  
	  typedef SpatialSearch::ResultConditionsContainerType              ResultConditionsContainerType;
      typedef SpatialSearch::VectorResultConditionsContainerType        VectorResultConditionsContainerType;

      typedef SpatialSearch::RadiusArrayType                            RadiusArrayType;
      typedef SpatialSearch::DistanceType                               DistanceType;
      typedef SpatialSearch::VectorDistanceType                         VectorDistanceType;	   
	  
      

      /// Default constructor.
      DemFemExplicitSolverStrategy(){}

      DemFemExplicitSolverStrategy(
                             ModelPart& fem_model_part,
							 ModelPart& dem_model_part,
							 const double fTimeStep,
							 const int   bDampOption,
							 const double DampRatio,
                             const bool MoveMeshFlag,
							 BaseType&   DEM_Strategy,
							 BaseType&   FEM_Strategy


      ): BaseType(dem_model_part, MoveMeshFlag)
      {
		 mpFem_model_part = &fem_model_part;
		 mpDem_model_part = &dem_model_part;
		 
		 mpDEM_Strategy   = &DEM_Strategy;
		 mpFEM_Strategy   = &FEM_Strategy;
		 
		 
		 mElementsAreInitialized = false;
		 
		 mStepFlag = 0;
      }
	  
	  
	  
      /// Destructor.
      virtual ~DemFemExplicitSolverStrategy()
      {
      }

      virtual void Initialize()
      {
               
		KRATOS_TRY
		
          
		KRATOS_CATCH("")
          

      }// Initialize()



      virtual double Solve()
      {
		  KRATOS_TRY
		  
          mpDEM_Strategy->Solve();
		  
		  
	      mStepFlag++;
        
          KRATOS_CATCH("") 
          
		  return 0.00;
          
      }//Solve()
	  
	  
	
	
	  ModelPart *mpFem_model_part;
	  ModelPart *mpDem_model_part;
	  BaseType * mpDEM_Strategy;
	  BaseType * mpFEM_Strategy;
	  
	  
	  bool mElementsAreInitialized;
	  
	  
	  int mStepFlag;


  }; // Class DemFemExplicitSolverStrategy


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined




