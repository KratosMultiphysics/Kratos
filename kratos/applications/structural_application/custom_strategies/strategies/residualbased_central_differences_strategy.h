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

#if !defined(KRATOS_RESIDUALBASED_CENTRAL_DIFERENCES_STRATEGY)
#define  KRATOS_RESIDUALBASED_CENTRAL_DIFERENCES_STRATEGY


/* System includes */
#include <limits>
#include<iostream>
#include<iomanip>



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


/* For contact analysis */
#include "custom_strategies/schemes/forward_increment_lagrange_multiplier_scheme.h"
#include "custom_utilities/boundary_conditions_and_contact_utilities.h"


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


	  ResidualBasedCentralDiferencesStrategy(ModelPart& model_part, int dimension,
			double alpha_damp,      /// para calcular la matriz de amortiguamiento
                        double fraction_delta_time,
                        double max_delta_time,
			bool   CalculateReactions,
			bool   ComputeContactConditions,
			bool   MoveMeshFlag
			)
	  : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
	      {


		std::cout<<std::fixed<<std::scientific<<std::setprecision(4);  
	        mdimension                 = dimension; 
		malpha_damp                = alpha_damp;
                mfraction_delta_time       = fraction_delta_time; 
                mmax_delta_time            = max_delta_time;
                mCalculateReactionsFlag    = CalculateReactions; 
	        mComputeContactConditions  = ComputeContactConditions;
                mElementsAreInitialized    = false;
		mConditionsAreInitialized  = false;
                mInitialConditions         = false;
                mCalculateOldTime          = false;  
	        mSolutionStepIsInitialized = false;
		mSolutionStepIsInitialized = false;
		mInitializeWasPerformed    = false;

		mtimestep                  = 0.00;
		
		if(mComputeContactConditions==true)
	  	    mpBCCU_Pointer =  typename BoundaryConditionsAndContactUtilities::Pointer (new BoundaryConditionsAndContactUtilities(model_part, mdimension) );
		
		mpLagrangianMultiplier = typename ForwardIncrementLagrangeMultiplierScheme::Pointer (new ForwardIncrementLagrangeMultiplierScheme(model_part, mdimension) ); 
		
                std::cout <<"TIME INTEGRATION METHOD  =  CENTRAL DIFFERENCES "<< std::endl;
	      }

	  virtual ~ResidualBasedCentralDiferencesStrategy () {}
	           
		   



double Solve() 
      { 

	#ifdef _OPENMP
	double start_prod = omp_get_wtime();   
	#endif
	std::cout<<"INITIALIZE SOLVE"<<std::endl;
	
	ModelPart& r_model_part              = BaseType::GetModelPart();
	ConditionsArrayType& pConditions     = r_model_part.Conditions();

	//ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();

	///Initializa los elementos y condiciones
	if(mInitializeWasPerformed == false){
	 Initialize();
	 minitial_conditions_size = pConditions.size();
	}

	///Initialize solution step
	if(mSolutionStepIsInitialized==false)
	   InitializeSolutionStep();

	///Computa el tiempo critico 
	ComputeCriticalTime();

	if(mInitialConditions==false){
	   ComputeInitialConditions();
	   GetForce();
	}

        /// Predict displacements
	/// Computa el nuevo desplazamiento a partir de los pasos anteriores
	ComputeIntermedialVelocityAndNewDisplacement();
	
	/// Actualizacion de los desplazamientos
	if(BaseType::MoveMeshFlag() == true) 
	   BaseType::MoveMesh();
   
	/// Compute contact force and displacements correcions
	/// realizando el bounding box
	/// verifico si hay condiciones de contacto y corrijo el desplazamiento
	if(mComputeContactConditions==true )
	  { 
                         /// Borra las antiguas condiciones de conatcto
	    mpBCCU_Pointer->CreateBoundariesAndLinkingConditions();      /// Crea las nuevas condiciones de contacto
	    double dist = CheckObjectContact(minitial_conditions_size);  /// crea los iteradores del contacto
	    
	    if (dist!=0)
	       Update();  
	  } 
	
	
	/// Computa los valores del paso de las velocidades y aceleraciones
	ComputeOldVelocitiesAndAccelerations();
	
	///Computing The new internal and external forces  for n+1 step
	GetForce();
	
	FinalizeSolutionStep();
	
	CheckConditionsStatus(minitial_conditions_size);

	      
	#ifdef _OPENMP
	double stop_prod = omp_get_wtime();
	std::cout << "TIME SOLVING  = " << stop_prod - start_prod << std::endl;
	#endif
	std::cout << "FINISHED SOLVE"<<std::endl;
	return 0;     
      }



// Borra las condiciones de contacto anteriores sin borrar las otras condiciones de las estructura
void CheckConditionsStatus(const unsigned int& initial_conditions_size)
{
      
      ModelPart& r_model_part              = BaseType::GetModelPart();
      ConditionsContainerType& pConditions = r_model_part.ConditionsArray();
   
      //WARNING = SOLO BORRAR LAS CONDIIONES DE CONTACTO 
      if(initial_conditions_size==0){
      // No bounday conditions in model part
      pConditions.clear();  
      }
      else{
      ConditionsContainerIterator  end_previos  = pConditions.begin() + initial_conditions_size; //- conditions;  
      ConditionsContainerIterator  end_actual   = pConditions.end();
      pConditions.erase(end_previos, end_actual);
      
      }
}


double CheckObjectContact(const unsigned int& initial_conditions_size)
{
  
      ModelPart& r_model_part                     = BaseType::GetModelPart();
      ConditionsContainerType& pConditions        = r_model_part.ConditionsArray();
      m_end_previos  = pConditions.begin() + initial_conditions_size;  
      m_end_actual   = pConditions.end();      
      return (std::distance(m_end_previos, m_end_actual)); 
}

/// comuta el D(n+1) a partir de V(n+1/2) 
void ComputeIntermedialVelocityAndNewDisplacement()
{
  
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    NodesArrayType& pNodes           = r_model_part.Nodes(); 
      
    const double current_delta_time  = CurrentProcessInfo[DELTA_TIME]; 
    const double mid_delta_time      = 0.50*(molddelta_time +  current_delta_time);

       
    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> node_partition;
    CreatePartition(number_of_threads, pNodes.size(), node_partition);

    array_1d<double,3> mid_neg_velocity; /// V(n-1/2)
    array_1d<double,3> mid_pos_velocity; /// V(n+1/2)
    //array_1d<double,3> new_displacement; /// V(n+1/2)
    
    #pragma omp parallel for private(mid_pos_velocity, mid_neg_velocity) 
    for(int k=0; k<number_of_threads; k++)
      {
	typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      

	{
	   
	  array_1d<double,3>& actual_displacement   = i->FastGetSolutionStepValue(DISPLACEMENT);   /// Estamos en paso  T(n+1)
	  array_1d<double,3>& current_displacement  = i->FastGetSolutionStepValue(DISPLACEMENT,1);   /// U(n)
	  array_1d<double,3>& old_displacement      = i->FastGetSolutionStepValue(DISPLACEMENT,2);   /// U(n-1)
	  array_1d<double,3>& current_Rhs           = i->FastGetSolutionStepValue(RHS);              /// Fext(n) - Fint(n)  
	  
	  
	  const double& nodal_mass    =  i->FastGetSolutionStepValue(NODAL_MASS); 
	  const double nodal_damping  =  malpha_damp * nodal_mass; 
	  const double factor_a       =  2.00 *  nodal_mass + nodal_damping *  current_delta_time;    
	  const double factor_b       =  2.00 *  nodal_mass - nodal_damping *  molddelta_time;  
	  const double factor_c       =  2.00 *  mid_delta_time;
	  
	  noalias(mid_neg_velocity)   = (1.00/molddelta_time) * (current_displacement - old_displacement);
	  noalias(mid_pos_velocity)   = (factor_b/factor_a) * mid_neg_velocity + (factor_c/factor_a) * current_Rhs ;

	  /// Update to the new displacement 
	  ///current_displacement = current_displacement + mid_pos_velocity * current_delta_time;
	  /// Una aproximacion de la velocidad
	  if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false )
	    actual_displacement[0]  = current_displacement[0] + mid_pos_velocity[0] * current_delta_time;

	  
	  if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
	    actual_displacement[1]     = current_displacement[1] + mid_pos_velocity[1] * current_delta_time;

	  if( i->HasDofFor(DISPLACEMENT_Z)){
	    if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
	       actual_displacement[2]     = current_displacement[2] + mid_pos_velocity[2] * current_delta_time; 
	}

      }   
    }
    
    
    std::cout<< "  OLD TIMESTEP                 = "<< molddelta_time     << "  SECONDS" << std::endl; 
    std::cout<< "  CURRENT TIMESTEP             = "<< current_delta_time << "  SECONDS" << std::endl; 
    std::cout<< "  AVERAGE TIMESTEP             = "<< mid_delta_time     << "  SECONDS" << std::endl; 
      
}

void Initialize()
{

    KRATOS_TRY

    /// Inicializando los elemtos
    if(mElementsAreInitialized == false)
	InitializeElements();

    /// Inicializando las condiciones
    if(mConditionsAreInitialized == false)
	InitializeConditions();

    mInitializeWasPerformed   = true;


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
	    unsigned int dim   = geom.WorkingSpaceDimension();
	    unsigned int index = 0;
	    for (unsigned int i = 0; i <geom.size(); i++)
	     {
		geom[i].SetLock();
		index = i*dim;
		geom[i].FastGetSolutionStepValue(NODAL_MASS) += MassMatrix(index,index);
		geom[i].UnSetLock();
	     }
	 }
      }
     r_model_part.GetCommunicator().AssembleCurrentData(NODAL_MASS); 
     mElementsAreInitialized   = true;
     KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************

///WARNING = Falta colocar el contacto
void InitializeConditions()
{
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
    
      
      mConditionsAreInitialized = true;
      
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

void ComputeCriticalTime()
{
    KRATOS_TRY
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();

    double delta_time_computed = 1.00;
    double delta_time_used     = 1.00; 
    double time                = 0.00;
    const int step             = CurrentProcessInfo[TIME_STEPS];  

    
    delta_time_computed = ComputeTime();
    delta_time_computed = Truncar_Delta_Time(delta_time_computed);
    if(delta_time_computed>mmax_delta_time) {delta_time_computed = mmax_delta_time;}
    delta_time_used = mfraction_delta_time * delta_time_computed;
    CurrentProcessInfo[DELTA_TIME] = delta_time_used; // reduzco el valor critico del tiempo en 75%
    
    mtimestep += delta_time_used; 
    r_model_part.CloneTimeStep(mtimestep);
    time = CurrentProcessInfo[TIME];
    
    if(mCalculateOldTime==false) {
        molddelta_time    = delta_time_used;
	mCalculateOldTime = true;
    }
    
    std::cout<< "  FACTOR DELTA CRITICAL TIME   = "<< mfraction_delta_time   << "        "  << std::endl;
    std::cout<< "  DELTA CRITICAL TIME COMPUTED = "<< delta_time_computed    << "  SECONDS" << std::endl; 
    std::cout<< "  DELTA TIME USED              = "<< delta_time_used        << "  SECONDS" << std::endl;
    std::cout<< "  CURRENT TIME                 = "<< time                   << "  SECONDS" << std::endl;
    std::cout<< "  TIME STEPS                   = "<< step                   << "         " << std::endl;
  
    KRATOS_CATCH("")

}

//***************************************************************************
//***************************************************************************

double ComputeTime()
{
  
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    ElementsArrayType& pElements     = r_model_part.Elements(); 
    
    int number_of_threads  =  0;
    //int thread_id          =  0;
    double delta_time_a    =  1.00;
    //double step            =  CurrentProcessInfo[TIME_STEPS];  

    #ifdef _OPENMP
       number_of_threads = omp_get_max_threads();
       //thread_id         = omp_get_thread_num();
    #else
       number_of_threads = 1;
       //thread_id = 0;
    #endif
    
    Vector values;
    vector<unsigned int> element_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);

    Vector dts(number_of_threads);
    for(int i = 0; i < number_of_threads; ++i)
    {
       dts[i] = 1.00;
    } 

    #pragma omp parallel for private(delta_time_a, values) shared(dts) 
    for(int k=0; k<number_of_threads; k++)
    {
    typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
    typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

    for (ElementsArrayType::iterator it=it_begin; it!= it_end; ++it)
       {
	  it-> Calculate(DELTA_TIME, delta_time_a, CurrentProcessInfo);
	  
	  // WARNING = Los threads Id siempre empiezan por cero. 
	  if(delta_time_a < dts[k])
	  dts[k] = delta_time_a;
        
       }
    }
    
   return  (*std::min_element(dts.begin(), dts.end())); 
  
}


//***************************************************************************
//***************************************************************************

void ComputeInitialConditions()
{
      KRATOS_TRY

      ModelPart& r_model_part         = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes          = r_model_part.Nodes(); 

      const double DeltaTime          = CurrentProcessInfo[DELTA_TIME]; 
      
      if(DeltaTime == 0)
	  KRATOS_ERROR(std::logic_error,"Detected delta_time = 0. Please check if the time step is created correctly for the current model part","");
    
      
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);
      //array_1d<double,3> OldDisplacement;
      #pragma omp parallel for //private(OldDisplacement) 
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

          /// debiera calcularse la aceleracion inicial segun libro de Chopra capitulo 5, tabla 5.3.1 
	  ///NOTA antes de hacer update (DISPLACEMENT) == (DISPLACEMENT,1);
	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
	  { 
	  const array_1d<double,3>& CurrentDisplacement  =  (i->FastGetSolutionStepValue(DISPLACEMENT));
	  const array_1d<double,3>& CurrentVelocity      =  (i->FastGetSolutionStepValue(VELOCITY));
	  const array_1d<double,3>& CurrentAcceleration  =  (i->FastGetSolutionStepValue(ACCELERATION));
	  array_1d<double,3>&       OldDisplacement      =  (i->FastGetSolutionStepValue(DISPLACEMENT,2));
	  
	  /// D_1 7.145 Libro de ""Estructuras Sometiadas a Acciones Sismicas"   
	  noalias(OldDisplacement) =  0.5*DeltaTime*DeltaTime*CurrentAcceleration - DeltaTime*CurrentVelocity + CurrentDisplacement;
	  //i->FastGetSolutionStepValue(DISPLACEMENT,2) = OldDisplacement; /// corresponde al time step = -1  
	  }
	}

         mInitialConditions  = true;
         KRATOS_CATCH("")
}


//***************************************************************************
//***************************************************************************

void GetForce()
{
      KRATOS_TRY
      
      ModelPart& r_model_part  = BaseType::GetModelPart();   

      /// Set to zero de RHS
      SetToZeroRHS();
      
      /// Compute the global external nodal force.
      Calculate_Conditions_RHS_and_Add();
      
      /// Compute the stress and body force of the element.
      Calculate_Elements_RHS_and_Add();
      
      r_model_part.GetCommunicator().AssembleCurrentData(RHS);
      
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
            array_1d<double,3>& reaction  = (i->FastGetSolutionStepValue(REACTION));
	    array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(RHS)); 
	    
	    noalias(reaction) = ZeroVector(3);
	    noalias(node_rhs) = ZeroVector(3);
	    
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
      unsigned int dim; 
      unsigned int index; 
     #pragma omp parallel for private (dim, index, rhs_cond)
      for(int k=0; k<number_of_threads; k++)
      {    
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
        for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
	{      
	  
	  Condition::GeometryType& geom = it->GetGeometry();  
	  it->CalculateRightHandSide(rhs_cond,CurrentProcessInfo);    
	  dim = geom.WorkingSpaceDimension();  
	  for (unsigned int i = 0; i <geom.size(); i++)
            {
	     
	    index = i*dim;
	    array_1d<double,3>& node_rhs = geom[i].FastGetSolutionStepValue(RHS);

	  
 	    for(unsigned int kk=0; kk<dim; kk++)
	    { geom[i].SetLock();
	     node_rhs[kk] += rhs_cond[index+kk];
             geom[i].UnSetLock();
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

      Vector rhs_elem;

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> element_partition;
      CreatePartition(number_of_threads, pElements.size(), element_partition);

      unsigned int dim;
      unsigned int index; 
      #pragma omp parallel for private (dim, index, rhs_elem)
      for(int k=0; k<number_of_threads; k++)
      {
        typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
         for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
         {
	  Element::GeometryType& geom = it->GetGeometry();
          dim = it->GetGeometry().WorkingSpaceDimension();
	  it->CalculateRightHandSide(rhs_elem, CurrentProcessInfo); 
	  for (unsigned int i = 0; i <geom.size(); i++)
	      {
		index = i*dim;
		array_1d<double,3>& node_rhs = geom[i].FastGetSolutionStepValue(RHS);
		for(unsigned int kk=0; kk<dim; kk++)
		    { geom[i].SetLock();      
		      node_rhs[kk] += rhs_elem[index+kk];
                      geom[i].UnSetLock();
		    }

	      }
	}
    }
    

     KRATOS_CATCH("")
}


//***************************************************************************
//***************************************************************************

void Calculate_Final_Force_Contribution(ModelPart::NodeIterator& i)  

{
        array_1d<double,3> Final_Force; 
        const double nodal_damping   =  malpha_damp * ((i)->FastGetSolutionStepValue(NODAL_MASS)); 

        array_1d<double,3>& node_rhs = (i->FastGetSolutionStepValue(RHS));    
	const array_1d<double,3>& velocity =  (i->FastGetSolutionStepValue(VELOCITY));
	
	noalias(Final_Force) = nodal_damping * velocity;    
	node_rhs -= Final_Force;  

}

//***************************************************************************
//***************************************************************************

/// Encuentra las velocidades y aceleraciones del paso n.
void ComputeOldVelocitiesAndAccelerations()
{ 
  
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    NodesArrayType& pNodes           = r_model_part.Nodes(); 
      
    const double current_delta_time  = CurrentProcessInfo[DELTA_TIME]; 
    const double mid_delta_time      = 0.50*(molddelta_time   +  current_delta_time);
    const double inv_sum_delta_time  = 1.00 / (molddelta_time + current_delta_time);

       
    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> node_partition;
    CreatePartition(number_of_threads, pNodes.size(), node_partition);

    array_1d<double,3> mid_neg_velocity; /// V(n-1/2)
    array_1d<double,3> mid_pos_velocity; /// V(n+1/2)
    array_1d<double,3> mid_neg_velocity_old; /// V(n-1/2)
    array_1d<double,3> mid_pos_velocity_old; /// V(n+1/2)
    
    #pragma omp parallel for private (mid_pos_velocity, mid_neg_velocity,  mid_neg_velocity_old,  mid_pos_velocity_old) 
    for(int k=0; k<number_of_threads; k++)
      {
	typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      

	{
	   
	  array_1d<double,3>& displacement          = i->FastGetSolutionStepValue(DISPLACEMENT);   /// Estamos en paso  T(n+1)
	  array_1d<double,3>& current_displacement  = i->FastGetSolutionStepValue(DISPLACEMENT,1);   /// U(n)
	  array_1d<double,3>& olddisplacement       = i->FastGetSolutionStepValue(DISPLACEMENT,2);   /// U(n-1)
	  array_1d<double,3>& veryolddisplacement   = i->FastGetSolutionStepValue(DISPLACEMENT,3);   /// U(n-1)
	  
	  noalias(mid_neg_velocity)       = (1.00/molddelta_time)     * (current_displacement   - olddisplacement);
	  noalias(mid_pos_velocity)       = (1.00/current_delta_time) * (displacement           - current_displacement);
	  noalias(mid_neg_velocity_old)   = (1.00/molddelta_time)     * (olddisplacement        - veryolddisplacement);
	  noalias(mid_pos_velocity_old)   = (1.00/current_delta_time) * (current_displacement   - olddisplacement);

            
	    array_1d<double,3>& velocity           = (i->FastGetSolutionStepValue(VELOCITY));  
 	    array_1d<double,3>& acceleration       = (i->FastGetSolutionStepValue(ACCELERATION));   
	    
	    /// X 
	    if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false )
	    {
	       
	       velocity[0]      =  inv_sum_delta_time * (displacement[0] - olddisplacement[0]);
	       acceleration[0]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[0] - mid_neg_velocity[0]);  
	       
	   }
            else
	    {
	       velocity[0]         =  inv_sum_delta_time * (current_displacement[0] - veryolddisplacement[0]);
	       acceleration[0]     =  (1.00 / mid_delta_time ) * (mid_pos_velocity_old[0] - mid_neg_velocity_old[0]);  
	    }
            
            /// Y
 	    if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
 	    {
	       velocity[1]      =  inv_sum_delta_time * (displacement[1] - olddisplacement[1]);
	       acceleration[1]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[1] - mid_neg_velocity[1]);  
 	    }
             else
 	    {
 	       velocity[1]         =  inv_sum_delta_time * (current_displacement[1] - veryolddisplacement[1]);
 	       acceleration[1]     =  (1.00 / mid_delta_time ) * (mid_pos_velocity_old[1] - mid_neg_velocity_old[1]);  
      
 	    }

            /// Z 
	    if( i->HasDofFor(DISPLACEMENT_Z)){
 	    if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
 	    {
	      velocity[2]      =  inv_sum_delta_time * (displacement[2] - olddisplacement[2]);
	      acceleration[2]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[2] - mid_neg_velocity[2]);  
 	    }
 	    else
 	    {
 	       velocity[2]         =  inv_sum_delta_time * (current_displacement[2] - veryolddisplacement[2]);
 	       acceleration[2]     =  (1.00 / mid_delta_time ) * (mid_pos_velocity_old[2] - mid_neg_velocity_old[2]);  
 	      
 	    }  
	  }  
      }    
  }
}
 

 
/// Computes the contact force and displacement corrections 
void Update()
{
   KRATOS_TRY
   mpLagrangianMultiplier->CalculateContactForceAndDisplacementCorrections(m_end_previos,m_end_actual); 
   KRATOS_CATCH("")
   return;
}

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

     molddelta_time = CurrentProcessInfo[DELTA_TIME]; 
    

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



///* WARNING = check with rossi
void CalculateReaction() 
{
 /*
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
*/
}

private:


unsigned int    mdimension;
unsigned int    minitial_conditions_size;  
bool   mInitialCalculations;
bool   mElementsAreInitialized;
bool   mConditionsAreInitialized;
bool   mInitialConditions;
bool   mCalculateReactionsFlag;  
bool   mReformDofSetAtEachStep;  
bool   mInitializeWasPerformed;
bool   mComputeContactConditions;
bool   mCalculateOldTime;
bool   mSolutionStepIsInitialized;
double malpha_damp;
double mfraction_delta_time;
double mmax_delta_time;
double molddelta_time;
double mtimestep;  /// la suma de los delta time


/// Contact Condition Iterations
BoundaryConditionsAndContactUtilities::Pointer mpBCCU_Pointer;
ForwardIncrementLagrangeMultiplierScheme::Pointer mpLagrangianMultiplier;
ConditionsContainerIterator m_end_previos;   
ConditionsContainerIterator m_end_actual;  


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
 
        return num;
 
    }



};  

} /* namespace Kratos.*/
#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME  defined */



