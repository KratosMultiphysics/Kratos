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
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2009-09-18 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/



#if !defined(PFEM2_EXPLICIT_STRATEGY)
#define  KRATOS_PFEM2_EXPLICIT_STRATEGY


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
#include "solving_strategies/strategies/explicit_strategy.h"
#include "solving_strategies/schemes/scheme.h"
//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "pfem_2_application.h"
//#include "custom_utilities/neighbours_calculator.h"
//#include "custom_elements/2fluid_2d.h"
//#include "custom_elements/2fluid_3d.h"

namespace Kratos
{

 template<
 class TSparseSpace,
 class TDenseSpace,
 class TLinearSolver>
 class PFEM2_Explicit_Strategy : public ExplicitStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {

	  public:

	  KRATOS_CLASS_POINTER_DEFINITION(PFEM2_Explicit_Strategy);

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




	  PFEM2_Explicit_Strategy(
	                ModelPart& model_part,
			const int        dimension,
           //             const int        damp_type,
			//const double     damping_ratio,
             //           const bool       virtual_mass,
              //          const double     contact_stiffness_ratio,
               //         const double     max_delta_time,
			//const bool       CalculateReactions,
			//const bool       ComputeFemFemContact,
			const bool       MoveMeshFlag
			//typename         TLinearSolver::Pointer pNewLinearSolver,
			//typename         TSchemeType::Pointer pScheme,
			//typename         TBuilderAndSolverType::Pointer pNewBuilderAndSolver
			)

	  : ExplicitStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, dimension ,MoveMeshFlag)
	      {

	      }

	  virtual ~PFEM2_Explicit_Strategy () {}



  //********************************************
  //********************************************
  inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
  }

  //SPECIFIC FUNCTIONS FOR MY APPLICATION
void InitializeSolutionStep() override
{
	KRATOS_TRY

	ModelPart& r_model_part  = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	switch ( CurrentProcessInfo[FRACTIONAL_STEP] )
	{
		case 0:
		{
			SetToZeroVariablesInViscousIterations(CurrentProcessInfo);
			break;
		}
		case 3:
		{
			SetToZeroVariablesInPresureIterations(CurrentProcessInfo);
			break;
		}
		case 4:
		{
			SetToZeroVariablesForVolumetricStrain(CurrentProcessInfo);
			break;
		}
		case 5:
		{
			SetToZeroVariablesForPressure(CurrentProcessInfo);
			break;
		}
		case 6:
		{
			SetToZeroVariablesInPresureViscousCorrection(CurrentProcessInfo);
			break;
		}
		case 7:
		{
			SetToZeroMassAndArea(CurrentProcessInfo);
			break;
		}
		case 10:
		{
			SetToZeroVariablesInPresureProjection(CurrentProcessInfo);
			break;
		}



		default:
		{
			KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", CurrentProcessInfo[FRACTIONAL_STEP]);
		}
	}
	KRATOS_CATCH("")
}

void FinalizeSolutionStep() override
{
	KRATOS_TRY

	ModelPart& r_model_part  = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	switch ( CurrentProcessInfo[FRACTIONAL_STEP] )
	{
		case 0:
		{
			UpdateLoopForViscousIterationsWithNormalization(CurrentProcessInfo);
			break;
		}
		case 3:
		{
			UpdateLoopForPressureIterationsWithNormalization(CurrentProcessInfo);
			break;
		}
		case 4:
		{
			UpdateLoopForVolumetricStrain(CurrentProcessInfo);
			break;
		}
		case 5:
		{
			UpdateLoopForPressure(CurrentProcessInfo);
			break;
		}
		case 6:
		{
			UpdateLoopForPressureViscousCorrection(CurrentProcessInfo);
			break;
		}
		case 7:
		{
			UpdateLoopForMassAndArea(CurrentProcessInfo);
			break;
		}
		case 10:
		{
			NormalizePressureProjection(CurrentProcessInfo);
			break;
		}

		default:
		{
			KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", CurrentProcessInfo[FRACTIONAL_STEP]);
		}
	}
	KRATOS_CATCH("")
}




//VISCOUS ITERATIONS

void SetToZeroVariablesInViscousIterations(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  //const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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
			noalias(i->FastGetSolutionStepValue(RHS)) = ZeroVector(3);
			i->FastGetSolutionStepValue(NODAL_MASS)=0.0;
		}
	  }

	  KRATOS_CATCH("")
}

void UpdateLoopForViscousIterationsWithNormalization(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double delta_t = CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
	  //const double factor = delta_t;

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
		   array_1d<double,3>& rhs = (i->FastGetSolutionStepValue(RHS));
		   array_1d<double,3>& node_update_variable = (i->FastGetSolutionStepValue(VELOCITY));
		   //noalias(node_update_variable) += factor* acceleration  ;
		   noalias(node_update_variable) = rhs /(i->FastGetSolutionStepValue(NODAL_MASS))  ;

	  }
	}

     KRATOS_CATCH("")
}


//PRESSURE ITERATIONS

void SetToZeroVariablesInPresureIterations(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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

  			if (iteration_number == 1)
				noalias(i->FastGetSolutionStepValue(ACCELERATION)) = ZeroVector(3);
			else// if (iteration_number==2)
				noalias(i->FastGetSolutionStepValue(ACCELERATION)) =   i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO); //second order cos we are in the second (or higher) iteration.
			//else if (iteration_number==3)
			//	noalias(inode->FastGetSolutionStepValue(ACCELERATION)) = - inode->FastGetSolutionStepValue(ACCELERATION) + inode->FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
						noalias(i->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3);

			noalias(i->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3);
			noalias(i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO)) = ZeroVector(3);


			//noalias(in->GetSolutionStepValue(PRESSURE,1))=in->FastGetSolutionStepValue(PRESSURE);
			//noalias(in->FastGetSolutionStepValue(PRESSURE)) = 0.0;
			i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
			i->FastGetSolutionStepValue(NODAL_MASS)=0.0;
			//i->FastGetSolutionStepValue(MASS)=ZeroVector(3);

			//i->FastGetSolutionStepValue(PRESSURE)=0.0;
			i->FastGetSolutionStepValue(VOLUMETRIC_STRAIN)=0.0;

		}
	  }

	  KRATOS_CATCH("")
}


void SetToZeroVariablesForVolumetricStrain(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  //const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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
			i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
			i->FastGetSolutionStepValue(VOLUMETRIC_STRAIN)=0.0;
		}
	  }

	  KRATOS_CATCH("")
}

void SetToZeroVariablesForPressure(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  //const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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
			i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
			i->FastGetSolutionStepValue(PRESSURE)=0.0;
		}
	  }

	  KRATOS_CATCH("")
}

void UpdateLoopForPressureIterationsWithNormalization(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double factor =  CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
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
		   //normalizing variables:
		   array_1d<double,3>& press_proj_no_ro = (i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO));
		   array_1d<double,3>& press_proj_stabilization = (i->FastGetSolutionStepValue(PRESS_PROJ));
		   double& mass = (i->FastGetSolutionStepValue(NODAL_MASS)); //this already includes 1/delta_t
		   double& area = (i->FastGetSolutionStepValue(NODAL_AREA));
		   press_proj_no_ro /= mass; //so this is already (pres_proj / mass) * delta_t;
		   press_proj_stabilization /= area;



		   //updating acceleration
		   array_1d<double,3>& acceleration = (i->FastGetSolutionStepValue(ACCELERATION));
		   noalias(acceleration) -= (press_proj_no_ro);
		   //updating variable
		   array_1d<double,3>& velocity = (i->FastGetSolutionStepValue(VELOCITY));
		   noalias(velocity) += (acceleration)  ; //the nodal mass includes the delta_time, so this is actually acceleration*delta_t

		   i->FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE)=i->FastGetSolutionStepValue(PRESSURE);
	  }
	}

     KRATOS_CATCH("")
}


//to calculate only the pressure projection:

void SetToZeroVariablesInPresureProjection(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  //const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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

			noalias(i->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3);

			i->FastGetSolutionStepValue(NODAL_AREA)=0.0;

		}
	  }

	  KRATOS_CATCH("")
}


void NormalizePressureProjection(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double factor =  CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
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
		   //normalizing variables:
		   //array_1d<double,3>& press_proj_no_ro = (i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO));
		   array_1d<double,3>& press_proj_stabilization = (i->FastGetSolutionStepValue(PRESS_PROJ));
		   //double& mass = (i->FastGetSolutionStepValue(NODAL_MASS));
		   //array_1d<double,3>& vectorial_mass=(i->FastGetSolutionStepValue(MASS));
		   double& area = (i->FastGetSolutionStepValue(NODAL_AREA));
		   //press_proj_no_ro /= mass;
		   //press_proj_no_ro(0) /= vectorial_mass(0);
		   //press_proj_no_ro(1) /= vectorial_mass(1);
		   //press_proj_no_ro(2) /= vectorial_mass(2)+1e-20;
		   press_proj_stabilization /= area;

		   //(i->FastGetSolutionStepValue(VOLUMETRIC_STRAIN))=i->FastGetSolutionStepValue(VOLUMETRIC_STRAIN)/area;
		   //(i->FastGetSolutionStepValue(PRESSURE))=i->FastGetSolutionStepValue(PRESSURE)/area;
		   //(i->FastGetSolutionStepValue(ELASTIC_PRESSURE))=i->FastGetSolutionStepValue(ELASTIC_PRESSURE)/area;

	  }
	}

     KRATOS_CATCH("")
}

void UpdateLoopForVolumetricStrain(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double factor =  CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
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
		   double& area = (i->FastGetSolutionStepValue(NODAL_AREA));
		   (i->FastGetSolutionStepValue(VOLUMETRIC_STRAIN))=i->FastGetSolutionStepValue(VOLUMETRIC_STRAIN)/area;
	  }
	}

     KRATOS_CATCH("")
}

void UpdateLoopForPressure(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double factor =  CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
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
		   double& area = (i->FastGetSolutionStepValue(NODAL_AREA));
		   (i->FastGetSolutionStepValue(PRESSURE))=i->FastGetSolutionStepValue(PRESSURE)/area;
	  }
	}

     KRATOS_CATCH("")
}




//PRESSURE VISCOUS CORRECTION

void SetToZeroVariablesInPresureViscousCorrection(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  //const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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
			i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
			i->FastGetSolutionStepValue(NODAL_MASS)=0.0;
		}
	  }

	  KRATOS_CATCH("")
}

void UpdateLoopForPressureViscousCorrection(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double factor =  CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
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
			i->FastGetSolutionStepValue(PRESSUREAUX) = i->FastGetSolutionStepValue(NODAL_AREA)/i->FastGetSolutionStepValue(NODAL_MASS);
			i->FastGetSolutionStepValue(EXTERNAL_PRESSURE) = i->FastGetSolutionStepValue(PRESSUREAUX) + i->FastGetSolutionStepValue(PRESSURE);
	  }
	}

     KRATOS_CATCH("")
}


void SetToZeroMassAndArea(ProcessInfo& CurrentProcessInfo)
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double factor =  CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
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
			i->FastGetSolutionStepValue(NODAL_MASS) = 0.0;
			i->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
	  }
	}

     KRATOS_CATCH("")
}

void UpdateLoopForMassAndArea(ProcessInfo& CurrentProcessInfo)
{
	KRATOS_TRY
	KRATOS_CATCH("")
}


};


 template<
 class TSparseSpace,
 class TDenseSpace,
 class TLinearSolver>
 class Fluid_Phase_PFEM2_Explicit_Strategy : public ExplicitStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {

	  public:

	  KRATOS_CLASS_POINTER_DEFINITION(Fluid_Phase_PFEM2_Explicit_Strategy);

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




	  Fluid_Phase_PFEM2_Explicit_Strategy(
	                ModelPart& model_part,
			const int        dimension,
           //             const int        damp_type,
			//const double     damping_ratio,
             //           const bool       virtual_mass,
              //          const double     contact_stiffness_ratio,
               //         const double     max_delta_time,
			//const bool       CalculateReactions,
			//const bool       ComputeFemFemContact,
			const bool       MoveMeshFlag
			//typename         TLinearSolver::Pointer pNewLinearSolver,
			//typename         TSchemeType::Pointer pScheme,
			//typename         TBuilderAndSolverType::Pointer pNewBuilderAndSolver
			)

	  : ExplicitStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, dimension ,MoveMeshFlag)
	      {

	      }

	  virtual ~Fluid_Phase_PFEM2_Explicit_Strategy () {}



  //********************************************
  //********************************************
  inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
  }

  //SPECIFIC FUNCTIONS FOR MY APPLICATION
void InitializeSolutionStep() override
{
	KRATOS_TRY

	ModelPart& r_model_part  = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	switch ( CurrentProcessInfo[FRACTIONAL_STEP] )
	{
		case 0:
		{
			SetToZeroVariablesInViscousIterations();
			break;
		}
		case 3:
		{
			SetToZeroVariablesInPresureIterations();
			break;
		}

		default:
		{
			KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", CurrentProcessInfo[FRACTIONAL_STEP]);
		}
	}
	KRATOS_CATCH("")
}

void FinalizeSolutionStep() override
{
	KRATOS_TRY

	ModelPart& r_model_part  = BaseType::GetModelPart();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	switch ( CurrentProcessInfo[FRACTIONAL_STEP] )
	{
		case 0:
		{
			UpdateLoopForViscousIterationsWithNormalization();
			break;
		}
		case 3:
		{
			UpdateLoopForPressureIterationsWithNormalization();
			break;
		}

		default:
		{
			KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", CurrentProcessInfo[FRACTIONAL_STEP]);
		}
	}
	KRATOS_CATCH("")
}




//VISCOUS ITERATIONS

void SetToZeroVariablesInViscousIterations()
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
	  //ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  //const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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
			noalias(i->FastGetSolutionStepValue(RHS)) = ZeroVector(3);
			i->FastGetSolutionStepValue(NODAL_MASS)=0.0;
		}
	  }

	  KRATOS_CATCH("")
}

void UpdateLoopForViscousIterationsWithNormalization()
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
	  //ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes   = r_model_part.Nodes();

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
		   array_1d<double,3>& rhs = (i->FastGetSolutionStepValue(RHS));
		   array_1d<double,3>& node_update_variable = (i->FastGetSolutionStepValue(WATER_VELOCITY));
		   noalias(node_update_variable) = rhs/(i->FastGetSolutionStepValue(NODAL_MASS)) ;
		   if (i->IsFixed(WATER_VELOCITY_X)==true)
				noalias(node_update_variable) = (i->FastGetSolutionStepValue(WATER_VELOCITY,1));

	  }
	}

     KRATOS_CATCH("")
}


//PRESSURE ITERATIONS

void SetToZeroVariablesInPresureIterations()
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
	  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes   = r_model_part.Nodes();

      //const double delta_t = CurrentProcessInfo[DELTA_TIME];
	  const int iteration_number = CurrentProcessInfo[NL_ITERATION_NUMBER];

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

  			if (iteration_number == 1)
				noalias(i->FastGetSolutionStepValue(ACCELERATION)) = ZeroVector(3);
			else// if (iteration_number==2)
				noalias(i->FastGetSolutionStepValue(ACCELERATION)) =   i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO); //second order cos we are in the second (or higher) iteration.
			//else if (iteration_number==3)
			//	noalias(inode->FastGetSolutionStepValue(ACCELERATION)) = - inode->FastGetSolutionStepValue(ACCELERATION) + inode->FastGetSolutionStepValue(PRESS_PROJ_NO_RO);

			//noalias(i->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3);
			noalias(i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO)) = ZeroVector(3);


			//noalias(in->GetSolutionStepValue(PRESSURE,1))=in->FastGetSolutionStepValue(PRESSURE);
			//noalias(in->FastGetSolutionStepValue(PRESSURE)) = 0.0;
			//i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
			i->FastGetSolutionStepValue(NODAL_MASS)=0.0;
			//i->FastGetSolutionStepValue(MASS)=ZeroVector(3);
		}
	  }

	  KRATOS_CATCH("")
}

void UpdateLoopForPressureIterationsWithNormalization()
{
      KRATOS_TRY

      ModelPart& r_model_part  = BaseType::GetModelPart();
	  //ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes   = r_model_part.Nodes();
	  //const double factor =  CurrentProcessInfo.GetValue(DELTA_TIME); //included in factor
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
		   //normalizing variables:
		   array_1d<double,3>& press_proj_no_ro = (i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO));
		   //array_1d<double,3>& press_proj_stabilization = (i->FastGetSolutionStepValue(PRESS_PROJ));
		   double& mass = (i->FastGetSolutionStepValue(NODAL_MASS));
		   //array_1d<double,3>& vectorial_mass=(i->FastGetSolutionStepValue(MASS));
		   //double& area = (i->FastGetSolutionStepValue(NODAL_AREA));
		   press_proj_no_ro /= mass;
		   //press_proj_no_ro(0) /= vectorial_mass(0);
		   //press_proj_no_ro(1) /= vectorial_mass(1);
		   //press_proj_no_ro(2) /= vectorial_mass(2)+1e-20;
		   //press_proj_stabilization /= area;

		   //updating acceleration
		   array_1d<double,3>& acceleration = (i->FastGetSolutionStepValue(ACCELERATION));
		   noalias(acceleration) -= (press_proj_no_ro);
		   //updating variable
		   array_1d<double,3>& velocity = (i->FastGetSolutionStepValue(WATER_VELOCITY));
		   noalias(velocity) += (acceleration)  ;

		   i->FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE)=i->FastGetSolutionStepValue(WATER_PRESSURE);
		   //i->FastGetSolutionStepValue(FRACT_VEL)=i->FastGetSolutionStepValue(VELOCITY);
		   if (i->IsFixed(WATER_VELOCITY_X)==true)
				noalias(velocity) = (i->FastGetSolutionStepValue(WATER_VELOCITY,1));
	  }
	}

     KRATOS_CATCH("")
}


};

} /* namespace Kratos.*/
#endif /* KRATOS_RESIDUALBASED_CENTRAL_DIFERENCES_STRATEGY */
