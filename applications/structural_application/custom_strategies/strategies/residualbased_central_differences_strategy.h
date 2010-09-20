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


#if !defined(KRATOS_RESIDUALBASED_CENTRAL_DIFERENCES_STRATEGY)
#define  KRATOS_RESIDUALBASED_CENTRAL_DIFERENCES_STRATEGY


/* System includes */


/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

namespace Kratos
{


template<class TSparseSpace,
class TDenseSpace, 
class TLinearSolver>
 
 class ResidualBasedCentralDiferencesStrategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {

	  public:

	  KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedCentralDiferencesStrategy);

	  typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

	  typedef typename BaseType::TDataType TDataType;

	  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

	  typedef typename BaseType::DofsArrayType DofsArrayType;

	  typedef typename Element::DofsVectorType DofsVectorType;

	  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

	  typedef typename BaseType::TSystemVectorType TSystemVectorType;

	  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

	  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

	  typedef ModelPart::NodesContainerType NodesArrayType;

	  typedef ModelPart::ElementsContainerType ElementsArrayType;

	  typedef ModelPart::ConditionsContainerType ConditionsArrayType;



	  typedef TSparseSpace SparseSpaceType;

	  typedef typename BaseType::TSchemeType TSchemeType;
	
	  typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
	  typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


	  ResidualBasedCentralDiferencesStrategy (ModelPart& model_part,  
			double alpha_damp,
                        double fraction_delta_time,
                        double max_delta_time,
			bool   CalculateReactions      = true,
			bool   MoveMeshFlag            = true
			)
	  : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
	      {

                			
			
			
			
		malpha_damp              = alpha_damp;
                mfraction_delta_time     = fraction_delta_time; 
                mmax_delta_time          = max_delta_time;
                mdelta_time              = 1.00;
                mCalculateReactionsFlag  = CalculateReactions; 
                mElementsAreInitialized  = false;
                mOldDisplacementComputed = false;
                mCalculateCriticalTime   = false;  
			
                 //set flags to start correcty the calculations
		 mSolutionStepIsInitialized = false;
		 mInitializeWasPerformed    = false;
                std::cout <<"USING THE CENTRAL DIFFERENCES  TIME INTEGRATION"<< std::endl;
	      }

	  virtual ~ResidualBasedCentralDiferencesStrategy (){}



double Solve() 
    { 
       ModelPart& r_model_part          = BaseType::GetModelPart();
       ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
      
      //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
      //if the operations needed were already performed this does nothing
      if(mInitializeWasPerformed == false)
           {Initialize();}

       //OPERATIONS THAT SHOULD BE DONE ONCE  
      // Orden es crucial
       if(mCalculateCriticalTime==false)
          {
           Compute_Critical_Time();
           }
       else
         {
            double step                    = CurrentProcessInfo[TIME_STEPS];   
            double delta_time_used         = mfraction_delta_time*mdelta_time;
            CurrentProcessInfo[DELTA_TIME] = delta_time_used; // reduzco el valor critico del tiempo en 75%
            r_model_part.CloneTimeStep(delta_time_used*step);
	    
	    std::cout<<"------------------------------------------------------------------"<<std::endl; 
	    std::cout<< "Factor Delta Critical Time   = "<< mfraction_delta_time << std::endl;
	    std::cout<< "Delta Critical Time Computed = "<< mdelta_time << std::endl; 
	    std::cout<< "Delta Time Used              = "<< delta_time_used << std::endl;
	    std::cout<< "Current Time                 = "<< time <<std::endl;
	    std::cout<< "Analysis Time Step           = "<< step <<std::endl;
            std::cout<<"------------------------------------------------------------------"<<std::endl;
	    
         }


     //Compute_Critical_Time();

      if(mOldDisplacementComputed==false)
         {ComputeOldDisplacement();}

      GetForce();

      GetNextDisplacement();

      Finalize();

      Update();

      if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

      //if(mCalculateReactionsFlag==true){CalculateReaction( ); }


       std::cout<<"FINISHED SOLVE"<<std::endl;
       return 0;     
    }

//***************************************************************************
//***************************************************************************


void Initialize()
{

KRATOS_TRY
//Initialize The Elements - OPERATIONS TO BE DONE ONCE
{InitializeElements();}
mInitializeWasPerformed  = true;
mElementsAreInitialized  = true;

KRATOS_CATCH("")
   

}


//***************************************************************************
//***************************************************************************

void InitializeElements()
{
      KRATOS_TRY
      // Accediendo a los elmentos
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
      #ifdef _OPENMP
      double start_prod = omp_get_wtime();   
      #endif
      std::cout<< "Initializing Elements " << std::endl;
      //KRATOS_WATCH(number_of_threads );
      //KRATOS_WATCH(element_partition );

      #pragma omp parallel for private(MassMatrix)
      for(int k=0; k<number_of_threads; k++)
      {

        typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
      {
	
        Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento
        (it)->Initialize(); 
        (it)->MassMatrix(MassMatrix, CurrentProcessInfo);

      unsigned int dim = geom.WorkingSpaceDimension();
      for (unsigned int i = 0; i <geom.size(); i++)
          {
	    geom[i].SetLock();
	    int index = i*dim;
	    geom[i].FastGetSolutionStepValue(NODAL_MASS) += MassMatrix(index,index);
	    geom[i].UnSetLock();
         }
       }
      }

     r_model_part.GetCommunicator().AssembleCurrentData(NODAL_MASS); 
     //mmax_delta_time = mdelta_time;
     #ifdef _OPENMP
     double stop_prod = omp_get_wtime();
     std::cout << "Time Initializing Elements  = " << stop_prod - start_prod << std::endl;
     #endif
     KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************


void Compute_Critical_Time()
{
    KRATOS_TRY
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    ElementsArrayType& pElements     = r_model_part.Elements(); 
    double delta_time_a    = 1.00;
    double delta_time_used = 1.00; 
    double step            = CurrentProcessInfo[TIME_STEPS];  
    double time = 0.00;

	
    #ifdef _OPENMP
       int number_of_threads = omp_get_max_threads();
    #else
       int number_of_threads = 1;
    #endif

    vector<unsigned int> element_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);
    //KRATOS_WATCH( element_partition );

    #ifdef _OPENMP
      double start_prod = omp_get_wtime();   
    #endif

    Vector dts(number_of_threads);
    for(int i = 0; i < number_of_threads; ++i)
    {
       dts[i] = 1.00;
    } 

    #pragma omp parallel for private(delta_time_a) shared(dts) 
    for(int k=0; k<number_of_threads; k++)
    {
    typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
    typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

    #ifdef _OPENMP
    int thread_id = omp_get_thread_num();
    #else
    int thread_id = 0;
    #endif

    for (ElementsArrayType::iterator it=it_begin; it!= it_end; ++it)
    {
       it->Calculate(DELTA_TIME, delta_time_a, CurrentProcessInfo);

    // WARNING = Los threads Id siempre empiezan por cero. 
       if(delta_time_a < dts[thread_id])
        {
    //KRATOS_WATCH(delta_time_a)
         dts[thread_id] = delta_time_a;
        } 
    }
    }

    std::cout<<"------------------------------------------------------------------"<<std::endl; 
    mdelta_time  =  (*std::min_element(dts.begin(), dts.end()));
    mdelta_time = Truncar_Delta_Time(mdelta_time);
    if(mdelta_time>mmax_delta_time) {mdelta_time = mmax_delta_time;}
    delta_time_used = mfraction_delta_time * mdelta_time;
    CurrentProcessInfo[DELTA_TIME] = delta_time_used; // reduzco el valor critico del tiempo en 75%
    r_model_part.CloneTimeStep( delta_time_used * step );
    time = CurrentProcessInfo[TIME];
    mCalculateCriticalTime = true;
    std::cout<< "Factor Delta Critical Time   = "<< mfraction_delta_time << std::endl;
    std::cout<< "Delta Critical Time Computed = "<< mdelta_time << std::endl; 
    std::cout<< "Delta Time Used              = "<< delta_time_used << std::endl;
    std::cout<< "Current Time                 = "<< time <<std::endl;
    std::cout<< "Analysis Time Step           = "<< step <<std::endl;

    #ifdef _OPENMP
    double stop_prod = omp_get_wtime();
    std::cout << "Time for calculating Delta Critical Time  = " << stop_prod - start_prod << std::endl;
    #endif
    std::cout<<"------------------------------------------------------------------"<<std::endl;
     
    KRATOS_CATCH("")

}


//***************************************************************************
//***************************************************************************
void ComputeOldDisplacement()
{
      KRATOS_TRY

      ModelPart& r_model_part         = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes          = r_model_part.Nodes(); 

      double DeltaTime     = CurrentProcessInfo[DELTA_TIME]; 
      double steps         = CurrentProcessInfo[TIME_STEPS];
      
      if(DeltaTime == 0)
	  KRATOS_ERROR(std::logic_error,"Detected delta_time = 0. Please check if the time step is created correctly for the current model part","");
    
      std::cout<< "Computing Old Displacements" << DeltaTime <<std::endl;
      std::cout<< "Delta Time    = "<< DeltaTime <<std::endl;
      std::cout<< "Analysis Step = "<< steps <<std::endl;     

      
	  #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif


      vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);
      
      #ifdef _OPENMP
      double start_prod = omp_get_wtime();  
      #endif

      array_1d<double,3> OldDisplacement;
      #pragma omp parallel for private(OldDisplacement) shared(DeltaTime)
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

     // debiera calcularse la aceleracion inicial segun libro de Chopra capitulo 5, tabla 5.3.1 
                 for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
		 //for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
		 //		i != r_model_part.NodesEnd() ; ++i)
		    { 
		      const array_1d<double,3>& CurrentDisplacement  =  (i->FastGetSolutionStepValue(DISPLACEMENT,1));
		      const array_1d<double,3>& CurrentVelocity      =  (i->FastGetSolutionStepValue(VELOCITY,1));
		      const array_1d<double,3>& CurrentAcceleration  =  (i->FastGetSolutionStepValue(ACCELERATION,1));

                      // if (i->Id()==392)
                      //    {KRATOS_WATCH(i->FastGetSolutionStepValue(DISPLACEMENT)) 
		      //     KRATOS_WATCH(i->FastGetSolutionStepValue(DISPLACEMENT,1))
                      //     KRATOS_WATCH(i->FastGetSolutionStepValue(DISPLACEMENT,2))  
                      //    }
		    // D_1 7.145 Libro de ""Estructuras Sometiadas a Acciones Sismicas"   
                     //KRATOS_WATCH(CurrentDisplacement)
		     OldDisplacement =  0.5*DeltaTime*DeltaTime*CurrentAcceleration - DeltaTime*CurrentVelocity + CurrentDisplacement;
		     i->FastGetSolutionStepValue(DISPLACEMENT,2) = OldDisplacement;
		    }
	}

	 #ifdef _OPENMP 
         double stop_prod = omp_get_wtime();
         std::cout << "Time for calculating Old Displacement  = " << stop_prod - start_prod << std::endl;
         #endif
         mOldDisplacementComputed  = true;
         KRATOS_CATCH("")
}


//***************************************************************************
//***************************************************************************

void GetForce()
{
      KRATOS_TRY
      // Set to zro de RHS
      Set_to_Zero_RHS();

      // Compute the stress and body force of the element.
      Calculate_Elements_RHS_and_Add();

      // Compute the global external nodal force.
      Calculate_Conditions_RHS_and_Add();
       //r_model_part.GetCommunicator().AssembleCurrentData(RHS);

      KRATOS_CATCH("")

}

//***************************************************************************
//***************************************************************************

void Set_to_Zero_RHS()

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
      //KRATOS_WATCH( number_of_threads );
      //KRATOS_WATCH( node_partition );

      #ifdef _OPENMP
      double start_prod = omp_get_wtime();  
      #endif

      #pragma omp parallel for 
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

     // debiera calcularse la aceleracion inicial segun libro de Chopra capitulo 5, tabla 5.3.1 
          for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      

           //for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	  //			i != r_model_part.NodesEnd() ; ++i)
	  {
	    //KRATOS_WATCH(i->Id())
	    //KRATOS_WATCH((i->FastGetSolutionStepValue(RHS)));
            array_1d<double,3>& reaction  = (i->FastGetSolutionStepValue(REACTION));
	    array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(RHS)); 
	    node_rhs[0] = 0.00; reaction[0] = 0.00;
	    node_rhs[1] = 0.00; reaction[1] = 0.00;   
	    node_rhs[2] = 0.00; reaction[2] = 0.00;
            
            //KRATOS_WATCH((i->FastGetSolutionStepValue(RHS)));
	  }
	}
               
         #ifdef _OPENMP
         double stop_prod = omp_get_wtime();
         std::cout << "Time for calculating Set to Zero RHS  = " << stop_prod - start_prod << std::endl;
         #endif

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

      #ifdef _OPENMP
      double start_prod = omp_get_wtime();
      #endif
  
      std::cout<< "Calculating RHS Conditions " << std::endl;
      //KRATOS_WATCH(number_of_threads );
      //KRATOS_WATCH(condition_partition );

     #pragma omp parallel for private (rhs_cond)
      for(int k=0; k<number_of_threads; k++)
      {    
	 // #ifdef _OPENMP
         // int thread_id = omp_get_thread_num();
         // #else
         // int thread_id = 1.00;
         // #endif
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)

      
      //for (ConditionsArrayType::iterator it=pConditions.begin(); it!=pConditions.end(); ++it)
	{      
	  Condition::GeometryType& geom = it->GetGeometry();
	  (it)->CalculateRightHandSide(rhs_cond,CurrentProcessInfo);
	    
	  unsigned int dim = geom.WorkingSpaceDimension();

//             Accediendo a los nodos por cada elemento   
//             Calculando masas y RHS poor nodos.    
          for (unsigned int i = 0; i <geom.size(); i++)
            {
	    int index = i*dim;
	    array_1d<double,3>& node_rhs = geom[i].FastGetSolutionStepValue(RHS);
 	    for(unsigned int kk=0; kk<dim; kk++)
	    { geom[i].SetLock();
	     node_rhs[kk] += rhs_cond[index+kk];
             geom[i].UnSetLock();
	     }
             }
        }
      }
 
        #ifdef _OPENMP
	double stop_prod = omp_get_wtime();
	std::cout << "Time for calculating Calculate_Conditions_RHS_and_Add  = " << stop_prod - start_prod << std::endl;
        #endif
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

      Vector rhs_elem;

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      //int thread_id = omp_get_thread_num();
      #else
      int number_of_threads = 1;
      //int thread_id = 1.00;
      #endif

      vector<unsigned int> element_partition;
      CreatePartition(number_of_threads, pElements.size(), element_partition);

      #ifdef _OPENMP
      double start_prod = omp_get_wtime();
      #endif  

      std::cout<< "Calculating RHS Elements " << std::endl;
      //KRATOS_WATCH(number_of_threads );
      //KRATOS_WATCH(element_partition );

     #pragma omp parallel for private (rhs_elem)
      for(int k=0; k<number_of_threads; k++)
      {

	  //#ifdef _OPENMP
          //int thread_id = omp_get_thread_num();
          //#else
          //int thread_id = 1.00;
          //#endif
	
        typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
  
         for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
     
        // Accediendo a los elmentos
       //for (ElementsArrayType::iterator it=pElements.begin(); it!=pElements.end(); ++it)
         {

	  Element::GeometryType& geom = it->GetGeometry();
	  unsigned int dim = (it)->GetGeometry().WorkingSpaceDimension();
	  (it)->CalculateRightHandSide(rhs_elem, CurrentProcessInfo);
    
	// Accediendo a los nodos por cada elemento   
	  // Calculando masas y RHS poor nodos.    
	  for (unsigned int i = 0; i <geom.size(); i++)
	      {
		unsigned int index = i*dim;
		array_1d<double,3>& node_rhs = geom[i].FastGetSolutionStepValue(RHS);
		
		for(unsigned int kk=0; kk<dim; kk++)
		    { geom[i].SetLock();      
		      node_rhs[kk] += rhs_elem[index+kk];
                      geom[i].UnSetLock();
		    }
	      }
	}
    }
    
         #ifdef _OPENMP
	 double stop_prod = omp_get_wtime();
         std::cout << "Time for calculating Calculate_Elements_RHS_and_Add  = " << stop_prod - start_prod << std::endl;
         #endif

     KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************

void GetNextDisplacement()
{

      KRATOS_TRY
      ModelPart& r_model_part = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes          = r_model_part.Nodes(); 
      //Element::EquationIdVectorType EquationId;


      double DeltaTime    = CurrentProcessInfo[DELTA_TIME];
      double factor_a     = (1.00/(DeltaTime*DeltaTime) + malpha_damp/(2.00*DeltaTime));

      array_1d<double,3> Final_Force;
      array_1d<double,3> NewDisp; 

      
	  #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

      vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);
      //KRATOS_WATCH( number_of_threads );
      //KRATOS_WATCH( node_partition );

      #ifdef _OPENMP
      double start_prod = omp_get_wtime();
      #endif
       #pragma omp parallel for  private(Final_Force, NewDisp) shared(factor_a)
       for(int k=0; k<number_of_threads; k++)
 	{
 	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
 	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
 
      // debiera calcularse la aceleracion inicial segun libro de Chopra capitulo 5, tabla 5.3.1 
           for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)     

//	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
//				i != r_model_part.NodesEnd() ; ++i)
	  {  

	     array_1d<double,3>&  NewDisplacement        = (i->FastGetSolutionStepValue(DISPLACEMENT)); 
	     array_1d<double,3>&  node_rhs               = (i->FastGetSolutionStepValue(RHS)); 
	     double Mass_Damping                         = (factor_a)*((i)->FastGetSolutionStepValue(NODAL_MASS));
	     //KRATOS_WATCH(factor_a)
	     //KRATOS_WATCH(DeltaTime)	    

	     // se le    suma la contribucion de los parametros de 7.152. Libro canet 
	     Calculate_Final_Force_Contribution(i ,Final_Force);
             //KRATOS_WATCH(node_rhs)
    

             /*
	     if(i->IsFixed(DISPLACEMENT_X) != true )
		{(NewDisp[0]) = node_rhs[0]/Mass_Damping;}
	     else NewDisp[0] = 0;

	     if (i->IsFixed(DISPLACEMENT_Y) != true  )
		{(NewDisp[1]) = node_rhs[1]/Mass_Damping;}
	     else NewDisp[1] = 0;

	     if( i->HasDofFor(DISPLACEMENT_Z)) {
		if( i->IsFixed(DISPLACEMENT_Z) != true  )
		    {(NewDisp[2]) = node_rhs[2]/Mass_Damping;}
		else NewDisp[2] = 0;
                     
		}
                */ 

             //   if (i->Id()==134)
 	     //   KRATOS_WATCH(i->FastGetSolutionStepValue(DISPLACEMENT_X,1))                  
  
                if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false )
		{NewDisp[0] = node_rhs[0]/Mass_Damping; NewDisplacement[0] = NewDisp[0];}
                //else
                //{ NewDisp[0] = i->FastGetSolutionStepValue(DISPLACEMENT_X); NewDisplacement[0] = NewDisp[0];}
		
                 if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
		{(NewDisp[1]) = node_rhs[1]/Mass_Damping; NewDisplacement[1] = NewDisp[1]; }
                //else
                // { NewDisp[1] = i->FastGetSolutionStepValue(DISPLACEMENT_Y); NewDisplacement[1] = NewDisp[1];}
		
               if( i->HasDofFor(DISPLACEMENT_Z)){
		if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
		{(NewDisp[2]) = node_rhs[2]/Mass_Damping; NewDisplacement[2] = NewDisp[2]; }
                //else
                // { NewDisp[2] = i->FastGetSolutionStepValue(DISPLACEMENT_Z,1); NewDisplacement[2] = NewDisp[2];}




                }
		//KRATOS_WATCH(i->Id())
                //KRATOS_WATCH(i->FastGetSolutionStepValue(NODAL_MASS))
                //KRATOS_WATCH(i->FastGetSolutionStepValue(RHS))
	        //KRATOS_WATCH(i->FastGetSolutionStepValue(DISPLACEMENT,0))
		//KRATOS_WATCH(i->FastGetSolutionStepValue(DISPLACEMENT,1))
                //KRATOS_WATCH(i->FastGetSolutionStepValue(DISPLACEMENT,2))
	        //noalias(NewDisplacement) = NewDisp; 

	
	  }
    }
       
         #ifdef _OPENMP
	 double stop_prod = omp_get_wtime();
         std::cout << "Time for calculating GetNextDisplacement  = " << stop_prod - start_prod << std::endl;
         #endif
         KRATOS_CATCH("")
	    
}

//***************************************************************************
//***************************************************************************

void Calculate_Final_Force_Contribution(
ModelPart::NodeIterator& i, array_1d<double,3>& Final_Force)  

{

	ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

	double DeltaTime    = CurrentProcessInfo[DELTA_TIME];  
	double factor_b     = (2.00/(DeltaTime*DeltaTime)); 
        double factor_c     = (-1.00/(DeltaTime*DeltaTime) + malpha_damp/(2.00*DeltaTime));
        double nodal_mass   = ((i)->FastGetSolutionStepValue(NODAL_MASS)); 


        array_1d<double,3>& node_rhs = (i->FastGetSolutionStepValue(RHS));    
	const array_1d<double,3>& CurrentDisplacement =  (i->FastGetSolutionStepValue(DISPLACEMENT,1));
        const array_1d<double,3>& OldDisplacement     =  (i->FastGetSolutionStepValue(DISPLACEMENT,2));

	Final_Force = nodal_mass*(factor_b*CurrentDisplacement + factor_c*OldDisplacement);
	node_rhs += Final_Force;  
}

//***************************************************************************
//***************************************************************************

void Update()
{

      ModelPart& r_model_part = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes          = r_model_part.Nodes(); 
      double DeltaTime    = CurrentProcessInfo[DELTA_TIME]; 

      array_1d<double, 3> OldAcc;
      array_1d<double, 3> OldVel;
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);
      //KRATOS_WATCH( number_of_threads );
      //KRATOS_WATCH( node_partition );

      #ifdef _OPENMP
      double start_prod = omp_get_wtime();
      #endif
      
      #pragma omp parallel for  private(OldAcc, OldVel) shared(DeltaTime)
      for(int k=0; k<number_of_threads; k++)
 	{
 	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
 	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
          for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)   
 //      for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
//				i != r_model_part.NodesEnd() ; ++i)
	  {
            // Las aceleraciones y las velocidades corresponden al paso anterior. No obstante se calculan como paso actual ya que no afecta en las ecuaciones
	    array_1d<double,3>& Oldacceleration    = (i->FastGetSolutionStepValue(ACCELERATION)); // WARNING = Estaba sin 1
            array_1d<double,3>& Oldvelocity        = (i->FastGetSolutionStepValue(VELOCITY));     // WARNING = Estaba sin 1
	    const array_1d<double,3>& currentdisp  = (i->FastGetSolutionStepValue(DISPLACEMENT)); 
            const array_1d<double,3>& oldDispl_1   = (i->FastGetSolutionStepValue(DISPLACEMENT,1));
            const array_1d<double,3>& oldDispl_2   = (i->FastGetSolutionStepValue(DISPLACEMENT,2));
         

            noalias(OldAcc) = (1.00/(DeltaTime*DeltaTime))*(currentdisp- 2.00*oldDispl_1 + oldDispl_2);
            noalias(OldVel) = (1.00/(2.00*DeltaTime))*(currentdisp - oldDispl_2); 

            noalias(Oldacceleration) =  OldAcc;
            noalias(Oldvelocity)     =  OldVel;

	  }
	  
	  
      }
      
      	 #ifdef _OPENMP
	 double stop_prod = omp_get_wtime();
         std::cout << "Time for Update= " << stop_prod - start_prod << std::endl;
         #endif
}


void Finalize()
    {
    KRATOS_TRY
    //finalizes solution step for all of the elements
    ModelPart& r_model_part = BaseType::GetModelPart();  
    ElementsArrayType& pElements = r_model_part.Elements();
    ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
    ConditionsArrayType& pConditions = r_model_part.Conditions();

    //mOldDisplacementComputed = false;  
          #ifdef _OPENMP
          int number_of_threads = omp_get_max_threads();
          #else
          int number_of_threads = 1;
          #endif

    vector<unsigned int> element_partition;
    vector<unsigned int> condition_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);
    CreatePartition(number_of_threads, pConditions.size(), condition_partition);
    #ifdef _OPENMP
    double start_prod = omp_get_wtime();
    #endif  
    //KRATOS_WATCH(number_of_threads );

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
         #ifdef _OPENMP
	 double stop_prod = omp_get_wtime();
         std::cout << "Time for Finalize = " << stop_prod - start_prod << std::endl;
         #endif
         KRATOS_CATCH("")
}



///* WARNING = check with rossi
void CalculateReaction( ) 
{
 
 ModelPart& r_model_part = BaseType::GetModelPart();   
 //ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
 //GetForce()
 for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
     {

            //KRATOS_WATCH(i->Id()) 
            double& mass                           = (i->FastGetSolutionStepValue(NODAL_MASS));
            array_1d<double,3>& reaction           = (i->FastGetSolutionStepValue(REACTION));
            array_1d<double,3>& rhs                = (i->FastGetSolutionStepValue(RHS,1));
            array_1d<double,3>& rhs_1              = (i->FastGetSolutionStepValue(RHS,2));
            array_1d<double,3>& acceleration       = (i->FastGetSolutionStepValue(ACCELERATION));
            array_1d<double,3>& velocity           = (i->FastGetSolutionStepValue(VELOCITY));
            array_1d<double,3> dif                =  rhs - rhs_1;
            //KRATOS_WATCH(dif)
            //KRATOS_WATCH(mass*acceleration)       
            //noalias(reaction)                      = -(mass*acceleration + mass*malpha_damp * velocity);  
             

            if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == true)
 		{reaction[0] = -dif[0] - mass * acceleration[0] - mass * malpha_damp * velocity[0];}
                
            if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == true )
		{reaction[1] = -dif[1] - mass * acceleration[1] - mass * malpha_damp * velocity[1];}
          
            if( i->HasDofFor(DISPLACEMENT_Z)){
		if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == true )
		{reaction[2] = -dif[2] - mass * acceleration[2] - mass * malpha_damp * velocity[2];}
                }
}
}

private:


bool   mInitialCalculations;
bool   mElementsAreInitialized;
bool   mOldDisplacementComputed;
bool   mCalculateReactionsFlag;  
bool   mReformDofSetAtEachStep;  
bool   mInitializeWasPerformed;
bool   mCalculateCriticalTime;
bool   mSolutionStepIsInitialized;
double malpha_damp;
double mfraction_delta_time;
double mdelta_time;
double mmax_delta_time;




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
      long double num_trucado = num; 
      unsigned long int a     = 1;
      unsigned long int i     = 10;
      while(trunc==false)
	{
	  num_trucado = num_trucado*i;
          a = a*i;
          if(num_trucado >= 100){
              num_trucado = static_cast<long unsigned int>(num_trucado);  
             
	      num         = num_trucado/a;
              trunc       = true;}
	}

        //KRATOS_WATCH(num)     
        return num;
 
    }



};  

} /* namespace Kratos.*/
#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME  defined */



