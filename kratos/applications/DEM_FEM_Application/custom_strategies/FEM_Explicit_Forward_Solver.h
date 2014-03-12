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


/// WARNING = Los desplazameitos son obtenidos en le paso
/// n+1; las velocidades  y aceleraciones son obtenidas para el paso n.

#if !defined(KRATOS_FEM_EXPLICIT_FORWARD_STRATEGY)
#define  KRATOS_FEM_EXPLICIT_FORWARD_STRATEGY


/* System includes */
#include <limits>
#include<iostream>
#include<iomanip>

/////////#define _OPENMP

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "spatial_containers/spatial_containers.h"
#include "dem_fem__application.h"


namespace Kratos
{
  
 template<
 class TSparseSpace,
 class TDenseSpace, 
 class TLinearSolver> 
 class FemExplicitForwardStrategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {

	  public:

	  KRATOS_CLASS_POINTER_DEFINITION(FemExplicitForwardStrategy);

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


          typedef Element::Pointer ParticlePointer;
          typedef typename std::vector<ParticlePointer> ParticlePointerVector;
          typedef typename std::vector<ParticlePointer>::iterator ParticlePointerIterator;

          typedef WeakPointerVector<Element > ParticleWeakVector;
          typedef WeakPointerVector<Element >::iterator ParticleWeakIterator;
		  



	  FemExplicitForwardStrategy
	       (
			ModelPart&       model_part, 
			const int        dimension,
			const int        damp_type,
			const double     damping_ratio,
			const bool       If_Virtual_Mass,
			const double     max_delta_time,
			const bool       CalculateReactions,
			const bool       MoveMeshFlag
			)
			
	  : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
	   {
			std::cout<< "******************************************************"<< std::endl;
	        std::cout <<"    EXPLICIT FORWARD SOLVER ANALYSIS FOR FEM"<< std::endl;
		    std::cout <<"     TIME INTEGRATION METHOD  =  Forward Differences    "<< std::endl;
			std::cout <<"IMPLEMENTED BY = Chun feng Based on the Version of Nelson"<< std::endl;
			std::cout<< "******************************************************"<< std::endl;
               
		
	        mdimension                 = dimension; 
			mCalculateReactionsFlag    = CalculateReactions;                
		    mInitializeWasPerformed    = false;
		    mComputeTime               = false;

		    mtimestep                  = max_delta_time;             

                
 			mdamp_type               = damp_type; 	
		    mdamping_ratio       = damping_ratio;
		    malpha_damp          = 0.00;
		    mbeta_damp           = 0.00; 

			mbIfvirtual_mass            = If_Virtual_Mass;
			mnumber_step             = 0;
			
		
	      }

	  virtual ~FemExplicitForwardStrategy () {}
	           
  double Solve()
  {

	KRATOS_TRY

	//ModelPart& r_model_part              = BaseType::GetModelPart();

	//ProcessInfo& CurrentProcessInfo      = r_model_part.GetProcessInfo();

	if(mInitializeWasPerformed == false)
	{
		 Initialize();
	}


   if(mIfHaveInitialzeSolutionStep == false)
	{
		InitializeSolutionStep();
		mIfHaveInitialzeSolutionStep = true;
	}



	////////Cfeng: Calcualte the deformation force, and add them to the node
	GetForce();


	///Cfeng: cal damping forces
	ComputeDampingForces();


	////////Cfeng:Calculate the movementment of the node;
	Calculate_Node_Movement_By_Newton_Law();

	///Cfeng:Cal reaction
	if(mCalculateReactionsFlag)
	{
		CalculateReaction();
	}


   /// Actualizacion de los desplazamientos
   if(BaseType::MoveMeshFlag() == true)
	{
		BaseType::MoveMesh();
	}



	///Finalize solution step
	FinalizeSolutionStep();



	return 0.00;
	KRATOS_CATCH("")
  }


void Initialize()
{

    KRATOS_TRY
	

    if(mInitializeWasPerformed == false)
	{
		
		InitializeElements();
	
		InitializeConditions();
		
		CalculateVirtualMass();
		

		mInitializeWasPerformed   = true;	
	}

    KRATOS_CATCH("")
}



void ComputeDampingForces()
{
    if(mdamp_type == 1)
    {
        ComputeNonViscousDampingForces();
    }
    else if(mdamp_type == 2)
    {
        ComputeViscousDampingForces();
    }
}

void ComputeViscousDampingForces()
{
}


void ComputeNonViscousDampingForces()
{   
	KRATOS_TRY
	ModelPart& r_model_part          = BaseType::GetModelPart();
	NodesArrayType& pNodes           = r_model_part.Nodes(); 
	
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

		for(ModelPart::NodeIterator it=i_begin; it!= i_end; ++it)
		{
			array_1d<double,3> &rhs        = it->FastGetSolutionStepValue(RHS);

			const array_1d<double,3> Vel = it->FastGetSolutionStepValue(VELOCITY);
			
			for (int iDof = 0; iDof < 3; iDof++)
			{
				if (Vel[iDof] > 0.0)
				{
					rhs[iDof] = rhs[iDof] - mdamping_ratio * fabs(rhs[iDof]);
				}
				else
				{
					rhs[iDof] = rhs[iDof] + mdamping_ratio * fabs(rhs[iDof]);
				}
			}
		}
	}
	  
	KRATOS_CATCH("")
}
 



void Calculate_Node_Movement_By_Newton_Law()
{

    KRATOS_TRY

    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    NodesArrayType& pNodes           = r_model_part.Nodes();


    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> node_partition;
    CreatePartition(number_of_threads, pNodes.size(), node_partition);

   // #pragma omp parallel for private(mid_pos_velocity, mid_neg_velocity)

    double UnBalForce    = 0.0;
	double ReactionForce = 0.0;


    #pragma omp parallel for reduction (+:UnBalForce, ReactionForce)
    for(int k=0; k<number_of_threads; k++)
	{
		typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
		typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

		for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
		{
			array_1d<double,3>& acceleration = (i->FastGetSolutionStepValue(ACCELERATION));
			array_1d<double,3>& velocity     = (i->FastGetSolutionStepValue(VELOCITY));
			array_1d<double,3>& displacement = i->FastGetSolutionStepValue(DISPLACEMENT);
			array_1d<double,3>& nodeForce    = i->FastGetSolutionStepValue(RHS);
			

	         double nodal_mass         =  i->FastGetSolutionStepValue(NODAL_MASS);
			 
			 if(mbIfvirtual_mass == true)
			 {
				 nodal_mass         =  i->FastGetSolutionStepValue(VIRTUAL_NODAL_MASS);
			 }

            acceleration = nodeForce / nodal_mass;

            array_1d<double,3> velocity_old = velocity;
            array_1d<double,3> velocity_new = velocity_old + acceleration * mtimestep;



            ////Cfeng : fixed node should set acc-vel-dis zero
            if( (i->pGetDof(VELOCITY_X))->IsFixed() == true)
	        {
                acceleration[0] = 0.0;
                displacement[0] += velocity[0] * mtimestep;

                ReactionForce = ReactionForce + fabs (nodeForce[0]);
            }
            else
            {
               velocity    [0]  = velocity_new[0];
               displacement[0] += velocity[0] * mtimestep;

               UnBalForce = UnBalForce + fabs(nodeForce[0]);
            }


            if( (i->pGetDof(VELOCITY_Y))->IsFixed() == true)
	       {
                acceleration[1] = 0.0;
                displacement[1] += velocity[1] * mtimestep;

                ReactionForce = ReactionForce + fabs (nodeForce[1]);
            }
            else
            {
               velocity    [1]  = velocity_new[1];
               displacement[1] += velocity[1] * mtimestep;

               UnBalForce = UnBalForce + fabs(nodeForce[1]);
            }


            if( (i->pGetDof(VELOCITY_Z))->IsFixed() == true)
	       {
                acceleration[2] = 0.0;
                displacement[2] += velocity[2] * mtimestep;

                ReactionForce = ReactionForce + fabs (nodeForce[2]);
            }
            else
            {
                velocity    [2]  = velocity_new[2];
                displacement[2] += velocity[2] * mtimestep;

                UnBalForce = UnBalForce + fabs(nodeForce[2]);
            }

        }
    }


    mUnbalRatioCal = UnBalForce / ReactionForce;

    if(ReactionForce < 1.0e-6)
    {
        mUnbalRatioCal = UnBalForce;
    }


    CurrentProcessInfo[DEM_FEM_CONVERGENCE_RATIO] = mUnbalRatioCal;

    KRATOS_CATCH("")
}


void CalculateVirtualMass()
{
    KRATOS_TRY

    if(mbIfvirtual_mass == true)
    {   
		ModelPart& r_model_part          = BaseType::GetModelPart();
		ElementsArrayType& pElements     = r_model_part.Elements();
		//ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
		NodesArrayType& pNodes           = r_model_part.Nodes();
		
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
				i->FastGetSolutionStepValue(VIRTUAL_NODAL_MASS) = 0.0;
			}
		}



		vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
            {
				
				double Young  = it->GetProperties()[YOUNG_MODULUS];        
				double Length = it->GetGeometry().Length();
				double Volume = 0.0;
				double VirtualMass = 0.0;

				Element::GeometryType& geom = it->GetGeometry();

	
				if (it->GetGeometry().Dimension() == 2 && geom.size() > 2)
				{
				   Volume = it->GetGeometry().Area();

				   VirtualMass = Young / (Length * Length) * Volume;
				} 
				else if (it->GetGeometry().Dimension() == 3 && geom.size() > 3 )
				{
				   Volume = it->GetGeometry().Volume();

				   VirtualMass = Young / (Length * Length) * Volume;
				}

				for (unsigned int i = 0; i <geom.size(); i++)
				{
				  double& mass = geom(i)->FastGetSolutionStepValue(VIRTUAL_NODAL_MASS);
				  geom(i)->SetLock();
				  mass  = mass + VirtualMass / (double)geom.size();
				  geom(i)->UnSetLock();
				}
            }
        }
    }


    KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************
 void InitializeElements()
    {
        KRATOS_TRY
        ModelPart& r_model_part          = BaseType::GetModelPart();
        ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements     = r_model_part.Elements();

        Matrix MassMatrix;
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        unsigned int index = 0;

        #pragma omp parallel for private(index, MassMatrix)
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
            {
				
                Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento
                (it)->Initialize();
                (it)->GetValue(IS_INACTIVE) = false;
                (it)->CalculateMassMatrix(MassMatrix, CurrentProcessInfo);
                const unsigned int& dim   = geom.WorkingSpaceDimension();
                index = 0;
                for (unsigned int i = 0; i <geom.size(); i++)
                {
                    double& mass = geom(i)->FastGetSolutionStepValue(NODAL_MASS);
                    geom(i)->SetLock();
                    index = i*dim;
                    mass  = mass + MassMatrix(index,index);
                    geom(i)->UnSetLock();
					
					
                }
            }
        }
        KRATOS_CATCH("")
    }

//***************************************************************************
//***************************************************************************

///WARNING = Falta colocar el contacto
void InitializeConditions()
{
  
       KRATOS_TRY
       
      ModelPart& r_model_part          = BaseType::GetModelPart();  
      //ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();   
      ConditionsArrayType& pConditions = r_model_part.Conditions();
      
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> condition_partition;
      CreatePartition(number_of_threads, pConditions.size(), condition_partition);
      
      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
	    (it) -> Initialize();
	  }

      }
      
      KRATOS_CATCH("")
      
}
  
//***************************************************************************
//***************************************************************************

void InitializeSolutionStep()
{
    KRATOS_TRY
    ModelPart& r_model_part          = BaseType::GetModelPart();  
    ElementsArrayType& pElements     = r_model_part.Elements();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    ConditionsArrayType& pConditions = r_model_part.Conditions();

    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> element_partition;
    vector<unsigned int> condition_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);
    CreatePartition(number_of_threads, pConditions.size(), condition_partition);

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
      typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
      typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	    (it) -> InitializeSolutionStep(CurrentProcessInfo);
	  }
    }

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
      {
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
	    (it) -> InitializeSolutionStep(CurrentProcessInfo);
	  }

    }


    
    KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************



//***************************************************************************
//***************************************************************************

    void GetForce()
    {
        KRATOS_TRY

        /// Set to zero de RHS
        SetToZeroRHS();

        /// Compute the global external nodal force.
        Calculate_Conditions_RHS_and_Add();

        /// Compute the stress and body force of the element. ( No lineal analysis)
        Calculate_Elements_RHS_and_Add();


        KRATOS_CATCH("")

    }

//***************************************************************************
//***************************************************************************

    void SetToZeroRHS()

    {
        KRATOS_TRY

        ModelPart& r_model_part  = BaseType::GetModelPart();
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
                array_1d<double,3>& normal    = (i->FastGetSolutionStepValue(NORMAL));
                array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(RHS));
                noalias(normal)               = ZeroVector(3);
                noalias(node_rhs)             = ZeroVector(3);
            }
        }

        KRATOS_CATCH("")
    }



//***************************************************************************
//***************************************************************************

    void Calculate_Conditions_RHS_and_Add()
    {
        KRATOS_TRY

        ModelPart& r_model_part          = BaseType::GetModelPart();
        ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        Vector rhs_cond;

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif
        vector<unsigned int> condition_partition;
        CreatePartition(number_of_threads, pConditions.size(), condition_partition);
        unsigned int index;
        #pragma omp parallel for private (index, rhs_cond)
        for(int k=0; k<number_of_threads; k++)
        {
            typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
            typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

            for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
            {
                Condition::GeometryType& geom = it->GetGeometry();
                it->CalculateRightHandSide(rhs_cond,CurrentProcessInfo);
                const unsigned int& dim = geom.WorkingSpaceDimension();
                for (unsigned int i = 0; i <geom.size(); i++)
                {
                    index = i*dim;
                    array_1d<double,3>& node_rhs = geom(i)->FastGetSolutionStepValue(RHS);
                    for(unsigned int kk=0; kk<dim; kk++)
                    {
                        geom(i)->SetLock();
                        node_rhs[kk] = node_rhs[kk] + rhs_cond[index+kk];
                        geom(i)->UnSetLock();
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }



//***************************************************************************
//***************************************************************************

    void Calculate_Elements_RHS_and_Add()
    {

        KRATOS_TRY
        ModelPart& r_model_part = BaseType::GetModelPart();
        ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements     = r_model_part.Elements();



#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        Vector rhs_elem;
        unsigned int index;
        #pragma omp parallel for private (index, rhs_elem)
        for(int k=0; k<number_of_threads; k++)
        {
            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
            {
                Element::GeometryType& geom = it->GetGeometry();
                const unsigned int& dim = it->GetGeometry().WorkingSpaceDimension();
                it->CalculateRightHandSide(rhs_elem, CurrentProcessInfo);
                for (unsigned int i = 0; i <geom.size(); i++)
                {
                    index = i*dim;
                    array_1d<double,3>& node_rhs = geom(i)->FastGetSolutionStepValue(RHS);
                    for(unsigned int kk=0; kk<dim; kk++)
                    {
                        geom(i)->SetLock();
                        node_rhs[kk] += /*node_rhs[kk]+*/ rhs_elem[index+kk];
                        geom(i)->UnSetLock();
                    }

                }
            }
        }

        KRATOS_CATCH("")
    }



//***************************************************************************
//***************************************************************************




//***************************************************************************
//***************************************************************************





void FinalizeSolutionStep()
    {
    KRATOS_TRY
    ModelPart& r_model_part = BaseType::GetModelPart();  
    ElementsArrayType& pElements = r_model_part.Elements();
    ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
    ConditionsArrayType& pConditions = r_model_part.Conditions();

    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> element_partition;
    vector<unsigned int> condition_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);
    CreatePartition(number_of_threads, pConditions.size(), condition_partition);
   

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
      typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
      typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	    (it) -> FinalizeSolutionStep(CurrentProcessInfo);
	  }
    }

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
      {
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
	    (it) -> FinalizeSolutionStep(CurrentProcessInfo);
	  }

    }
    
    KRATOS_CATCH("")

}



void CalculateReaction() 
{
    ModelPart& r_model_part = BaseType::GetModelPart();   
    NodesArrayType& pNodes  = r_model_part.Nodes(); 
    
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
         //for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; ++i)
         {
	        array_1d<double,3>& reaction                 = (i->FastGetSolutionStepValue(REACTION));
			double mass                                  = (i->FastGetSolutionStepValue(NODAL_MASS));
			
			if(mbIfvirtual_mass == true)
			 {
				 mass         =  i->FastGetSolutionStepValue(VIRTUAL_NODAL_MASS);
			 }
			
            const array_1d<double,3>& rhs                = (i->FastGetSolutionStepValue(RHS));
            const array_1d<double,3>& acceleration       = (i->FastGetSolutionStepValue(ACCELERATION));
            const array_1d<double,3>& velocity           = (i->FastGetSolutionStepValue(VELOCITY));
            //array_1d<double,3> dif                 =  rhs; 

            if( (i->pGetDof(VELOCITY_X))->IsFixed() == true)
 		    {
				reaction[0] =  -rhs[0] + mass * acceleration[0] + mass * malpha_damp * velocity[0];
			}
				
                
            if( i->pGetDof(VELOCITY_Y)->IsFixed() == true )
		    {
				reaction[1] = -rhs[1] + mass * acceleration[1] + mass * malpha_damp * velocity[1];
			}
          
            if( i->HasDofFor(VELOCITY_Z))
	        {
		        if( i->pGetDof(VELOCITY_Z)->IsFixed() == true )
		        {
					reaction[2] = -rhs[2] + mass * acceleration[2] + mass * malpha_damp * velocity[2];}
                }
            }
      } 
 }


private:


unsigned int    mdimension;


bool   mCalculateReactionsFlag;   
bool   mInitializeWasPerformed;
bool   mComputeTime;
double mdamping_ratio;
double malpha_damp;
double mbeta_damp; 
double mtimestep;  /// la suma de los delta time

int    mdamp_type;
bool   mbIfvirtual_mass;

int    mnumber_step;

double mUnbalRatioCal;
bool   mIfHaveCalVirtualMass;


bool   mIfHaveInitialzeSolutionStep;



//******************************************************************************************
//******************************************************************************************
inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
  }


inline double Truncar_Delta_Time(double& num)
    {
      bool trunc = false;
      double num_trucado = num; 
      unsigned long int a     = 1;
      unsigned long int i     = 10;
      if(num!=0.00){
      while(trunc==false)
	{
	  num_trucado = num_trucado*i;
          a = a*i;
          if(num_trucado >= 1){
              num_trucado = static_cast<long unsigned int>(num_trucado);  
              num         = num_trucado/a;
              trunc       = true;
	   }
 	  }
	}
 
        return num;
 
    }

};

} /* namespace Kratos.*/
#endif /* KRATOS_FEM_EXPLICIT_FORWARD_STRATEGY */


