//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//




#if !defined(EXPLICIT_STRATEGY)
#define  EXPLICIT_STRATEGY


/* System includes */
#include <string>
#include <iostream> 
#include <algorithm>

/////////#define _OPENMP

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/element.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"


namespace Kratos
{
  
 template<
 class TSparseSpace,
 class TDenseSpace, 
 class TLinearSolver> 
 class ExplicitStrategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {

	  public:

	  KRATOS_CLASS_POINTER_DEFINITION(ExplicitStrategy);

	  typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

	  typedef typename BaseType::TDataType TDataType;
	  
	  typedef TSparseSpace SparseSpaceType;

	  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
	  
	  typedef typename BaseType::TSchemeType TSchemeType;

	  typedef typename BaseType::DofsArrayType DofsArrayType;

	  typedef typename Element::DofsVectorType DofsVectorType;

	  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

	  typedef typename BaseType::TSystemVectorType TSystemVectorType;

	  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

	  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

	  typedef ModelPart::NodesContainerType NodesArrayType;

	  typedef ModelPart::ElementsContainerType ElementsArrayType;

	  typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	  
	  typedef ModelPart::ConditionsContainerType::ContainerType ConditionsContainerType;
      
	  typedef ConditionsContainerType::iterator                 ConditionsContainerIterator;

	  typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

	  typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

          typedef ModelPart::PropertiesType PropertiesType;


          //typedef Element::Pointer ParticlePointer;
          //typedef typename std::vector<ParticlePointer> ParticlePointerVector;
          //typedef typename std::vector<ParticlePointer>::iterator ParticlePointerIterator;

          //typedef WeakPointerVector<Element > ParticleWeakVector;
          //typedef WeakPointerVector<Element >::iterator ParticleWeakIterator;
	  



	  ExplicitStrategy(
	                ModelPart& model_part, 
			const int        dimension,
			const bool       move_mesh_flag
			)
			
	  : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, move_mesh_flag)
	      {
			std::cout<< "*************************************"<< std::endl;
	        std::cout <<"*   EXPLICIT CALCULATIONS STRATEGY  *"<< std::endl;
            std::cout<< "*************************************"<< std::endl;
       
	      }

	  ~ExplicitStrategy () override {}
	           


//***************************************************************************
//***************************************************************************

void AssembleLoop()
{

	KRATOS_TRY
	ModelPart& r_model_part = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
	ElementsArrayType& pElements = r_model_part.Elements();
  

	typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() ;
	typename ElementsArrayType::iterator it_end   = pElements.ptr_end();
	for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
	{
	  //	Element::GeometryType& geom = it->GetGeometry();

	   //for (unsigned int i = 0; i < geom.size(); i++)
	   //			 geom(i)->SetLock();

		
		it->AddExplicitContribution(CurrentProcessInfo);

	   //for (unsigned int i = 0; i < geom.size(); i++)
	   //			 geom(i)->UnSetLock();

	}


	KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************
        
void NormalizeVariable(const Variable<array_1d<double, 3 > >& rRHSVariable, const Variable<double >& rNormalizationVariable)
{
      KRATOS_TRY
      
      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes(); 
	  //const double delta_t = CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      DenseVector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);

      #pragma omp parallel for 
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

      for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	  {
		   array_1d<double,3>& node_rhs_variable = (i)->FastGetSolutionStepValue(rRHSVariable);
		   double& normalization_variable = (i)->FastGetSolutionStepValue(rNormalizationVariable);
		   
		   node_rhs_variable /= normalization_variable;		   
	  }
	}
                             
     KRATOS_CATCH("")
}

void ExplicitUpdateLoop(const Variable<array_1d<double, 3 > >& rUpdateVariable, const Variable<array_1d<double, 3 > >& rRHSVariable, const double& factor)
{
      KRATOS_TRY
      
      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes(); 
	  //const double delta_t = CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      DenseVector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);

      #pragma omp parallel for 
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

      for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	  {
		   array_1d<double,3>& node_update_variable = (i)->FastGetSolutionStepValue(rUpdateVariable);
		   array_1d<double,3>& node_rhs_variable = (i)->FastGetSolutionStepValue(rRHSVariable);
		   noalias(node_update_variable) += factor* node_rhs_variable  ;
		   
	  }
	}
                             
     KRATOS_CATCH("")
}

inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, DenseVector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
  }
  
  
  
  
  
  
  
  
  //********************************************
  //********************************************
void InitializeSolutionStep() override
{
	KRATOS_TRY

	ModelPart& r_model_part = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
	ElementsArrayType& pElements = r_model_part.Elements();
  

	typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() ;
	typename ElementsArrayType::iterator it_end   = pElements.ptr_end();
	for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
	{
		it->InitializeSolutionStep(CurrentProcessInfo);
	}
	
	KRATOS_CATCH("")
}

void FinalizeSolutionStep() override
{
	KRATOS_TRY

	ModelPart& r_model_part = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
	ElementsArrayType& pElements = r_model_part.Elements();
  

	typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() ;
	typename ElementsArrayType::iterator it_end   = pElements.ptr_end();
	for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
	{
		it->FinalizeSolutionStep(CurrentProcessInfo);
	}
	
	KRATOS_CATCH("")
}

};

} /* namespace Kratos.*/
#endif /* EXPLICT_STRATEGY */

