/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
/* *********************************************************   
*          
*   Last Modified by:    $Author:  $
*   Date:                $Date: 2009-09-18 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/



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

	  virtual ~ExplicitStrategy () {}
	           


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

      vector<unsigned int> node_partition;
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

      vector<unsigned int> node_partition;
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

inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
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
void InitializeSolutionStep()
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

void FinalizeSolutionStep()
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

