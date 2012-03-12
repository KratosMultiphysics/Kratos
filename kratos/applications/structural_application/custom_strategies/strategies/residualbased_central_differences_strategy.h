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
#include "custom_utilities/joint.h"
#include "custom_utilities/disconnect_utility.h"

namespace Kratos
{
  
 enum Constraint_Enforcement{Penalty_Methods, Lagrange_Multiplier_Methods};  
 template<
 class TSparseSpace,
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
	  
	  typedef BoundaryConditionsAndContactUtilities         BoundaryAndContactType;
	  
	  typedef ForwardIncrementLagrangeMultiplierScheme      LagrangeMultiplierType;

	   typedef Joint<4> Joint2D;
	   
	  
	  /// Notes:
	  /*
	  
	   a) Penalty methods corresponf to elastic colision and therefore the kinectic energy is preservated.
	      It need some damping to disipated kinetic energy.
	      
	   b) Using lagrangian multiplier aproach, the kinetic energy is disipated during the contact 
	   iteraction and therefore no contact damping is required.
	  */
	  
	  
	     ///@param 
	     ///CE                   Constraint_Enforcement&  
	     ///alpha_damp,          Para calcular la matriz de amortiguamiento proporcional a la masa
	     ///fraction_delta_time  Fracion del timepo critico (0-1)
	     ///max_delta_time,      Maximo incremento de tiempo
	     
	  ResidualBasedCentralDiferencesStrategy(
	                ModelPart& model_part, 
			const Constraint_Enforcement& CE,  
			const int        dimension,
			const double     damping_ratio, /// para calcular la matriz de amortiguamiento proporcional a la masa
                        const double     fraction_delta_time,
                        const double     max_delta_time,
		        const double     penalty_factor, 
			const bool       CalculateReactions,
			const bool       ComputeContactConditions,
			const bool       MoveMeshFlag
			)
	  : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
	      {
		
	        std::cout <<"DYNAMIC SOLVER ANALYSIS FOR COMBINED FINITE AND DISCRET ELEMENT METHODS "<< std::endl;
                std::cout <<"TIME INTEGRATION METHOD  =  CENTRAL DIFFERENCES    "<< std::endl;
                std::cout <<"IMPLEMENTED BY           =  ING. NELSON LAFONTAINE "<< std::endl;
		
		
	        mdimension                 = dimension; 
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
		mClearContactConditions    = false;
		mLocalSearchInitialize     = false;
		minitial_conditions_size   = 0;
		mcontact_conditions_size   = 0; 
		mtimestep                  = 0.00;
		mCE                        = CE;
	  	mpBCCU_Pointer             = typename BoundaryAndContactType::Pointer (new BoundaryAndContactType(model_part, mdimension, penalty_factor));
		
		if(mCE==Lagrange_Multiplier_Methods){
		 mpLagrangianMultiplier     = typename ForwardIncrementLagrangeMultiplierScheme::Pointer (new ForwardIncrementLagrangeMultiplierScheme(model_part, mdimension) ); 
		}
	
		mdamping_coeficients = false;
		mdamping_ratio       = damping_ratio;
		malpha_damp          = 0.00;
		mbeta_damp           = 0.00; 
		mpenalty_factor      = penalty_factor;
		
		mDTU.CreateJoints(model_part,mdimension);
		
	      }

	  virtual ~ResidualBasedCentralDiferencesStrategy () {}
	           
		   



double Solve() 
      {
	KRATOS_TRY
        std::cout<<std::fixed<<std::scientific<<std::setprecision(6);
	
	#ifdef _OPENMP
	double start_prod = omp_get_wtime();   
	#endif
	std::cout<<"INITIALIZE SOLVE"<<std::endl;
	
	ModelPart& r_model_part              = BaseType::GetModelPart();
	ConditionsArrayType& pConditions     = r_model_part.Conditions();

	
	
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
	   minitial_conditions_size = pConditions.size();
	   ComputeInitialConditions();
	   GetForce();
	}
        /// Predict displacements
	/// Computa el nuevo desplazamiento a partir de los pasos anteriores en n+1
	ComputeIntermedialVelocityAndNewDisplacement();
	

	
	/// Computa los valores del paso de las velocidades y aceleraciones	
       ComputeOldVelocitiesAndAccelerations();
	
	if(mCalculateReactionsFlag)
	    CalculateReaction(); 
	
	/// Computa los boundaries del modelo antes de la atualizacion de los desplazamientos
	/*
 	if(mComputeContactConditions==true )
 	  { 
	    //mpBCCU_Pointer->CreateBoundaries(minitial_conditions_size);   ///crea las superficies de contacto 
// 	    if(mCe==Lagrange_Multiplier_Methods) 
//  	       mpBCCU_Pointer->LocalSearch();
 	  }
 	  */
	    
	/// Actualizacion de los desplazamientos
 	if(BaseType::MoveMeshFlag() == true) 
 	   BaseType::MoveMesh();
   	
        ///Computing The new internal and external forces  for n+1 step
	GetForce();
	
	/// Discontinum Mechanics
	ComputeInterfaceForces();

	
	/// Fragmentation and fracture for DEM
	//Heuristic_Formula(mDTU.Begin(), mDTU.End());
        
	//WARNING = To be checked
	//ComputeDampingForces();
	
	/// Compute contact force and displacements correcions
	/// realizando el bounding box
	/// verifico si hay condiciones de contacto y corrijo el desplazamiento
	if(mComputeContactConditions==true)
	  { 
	    if(mCE==Penalty_Methods)
	    {
	       //const bool rflag = false;
	       //ResetFlagComputeBoundaryContour(rflag)
	       //mpBCCU_Pointer->Clear(minitial_conditions_size); //CreateBoundaries(minitial_conditions_size); 
	       mpBCCU_Pointer->CreateBoundaries(minitial_conditions_size); /// lista de los elementos de contorno 
	       mpBCCU_Pointer->ComputeContactForce();
	    }
	    
	    if(mCE==Lagrange_Multiplier_Methods)  
	    {
 	    ConditionsContainerIterator end_previos;   
            ConditionsContainerIterator end_actual;
            mpBCCU_Pointer->CreateBoundaries(minitial_conditions_size);                                        /// Crea las superficies de contacto
 	    mcontact_conditions_size = pConditions.size();                                                     /// Master Condition size + condiciones de inicio  
 	    mpBCCU_Pointer->CreateLinkingConditionsBasedOnLocalSearch(minitial_conditions_size);               /// Crea las nuevas condiciones de contacto
 	    int dist = CheckObjectContact(mcontact_conditions_size, end_previos, end_actual);                  /// crea los iteradores del contacto
 	    
 	    if(dist!=0)
 	       Update(end_previos,end_actual);  
	    }
	  }
	  
	FinalizeSolutionStep();
	
	if(mComputeContactConditions==true){   
	    if(mCE==Lagrange_Multiplier_Methods)  
	    {
	     ComputeNormalContactReactions(mcontact_conditions_size);
	     CheckConditionsStatus(mcontact_conditions_size);
	    }
	}
	
	
	/// Computing energies
	CalculateEnergies();

	
	#ifdef _OPENMP
	double stop_prod = omp_get_wtime();
	std::cout << "TIME SOLVING                   = "<<   stop_prod - start_prod    << "  SECONDS" << std::endl; 
	#endif
	std::cout << "FINISHED SOLVE"<<std::endl;
	return 0.00;   
	
	KRATOS_CATCH("")
      }



void ComputeDampingForces()
{
   ComputeViscousDampingForces();
}

void ComputeViscousDampingForces()
{
   ComputeDampingForcesWithMass();
   ComputeDampingForcesWithStiffness();
}


void ComputeNonViscousDampingForces()
{
  
  KRATOS_TRY
      ModelPart& r_model_part          = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();  
      //ElementsArrayType& pNodes        = r_model_part.Nodes(); ///  WARNING

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
       #endif

      vector<unsigned int> nodes_partition;
      CreatePartition(number_of_threads, pNodes.size(), nodes_partition);

      const double damp  = 0.00;
      array_1d<double, 3> DampingForces;
      
      typename NodesArrayType::iterator vec_begin = r_model_part.Nodes().ptr_begin()

      #pragma omp parallel for private(DampingForces) 
      for(int k=0; k<number_of_threads; k++)
      {
	typename NodesArrayType::iterator it_begin= vec_begin + nodes_partition[k];
	typename NodesArrayType::iterator it_end= vec_begin  + nodes_partition[k+1];
	for (NodesArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	      array_1d<double,3>& rhs = it->FastGetSolutionStepValue(RHS); 
	      array_1d<double,3>& Vel = it->FastGetSolutionStepValue(VELOCITY); 
	      noalias(DampingForces)  = -(damp * norm_2(rhs) / norm_2(Vel) ) * Vel;
              noalias(rhs) += DampingForces;
	  }  
	}
	    KRATOS_CATCH("")
    }
 



void ComputeDampingForcesWithMass()
{
  
    KRATOS_TRY
      ModelPart& r_model_part          = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();  
      ElementsArrayType& pElements     = r_model_part.Elements(); 

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
       #endif

      vector<unsigned int> element_partition;
      CreatePartition(number_of_threads, pElements.size(), element_partition);

      unsigned int index = 0; 
      unsigned int dim_2 = 0; 
      Matrix rLeftHandSideMatrix;
      Vector rRightHandSideVector;
      Vector Velocities;
      
      #pragma omp parallel for private(index, dim_2, rLeftHandSideMatrix,Velocities, rRightHandSideVector) 
      for(int k=0; k<number_of_threads; k++)
      {
	typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
	for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	    /// Las fuerzas de contacto a los elementos que estan en contacto
	    if((it)->GetValue(IS_TARGET)==true)
	    {
		it->MassMatrix(rLeftHandSideMatrix, CurrentProcessInfo);
		Element::GeometryType& geom         = it->GetGeometry();
		const unsigned int& dim             = geom.WorkingSpaceDimension();
		const unsigned int& number_of_nodes = geom.size();
		index                               = 0;
		dim_2                               =  dim * number_of_nodes;
		Velocities.resize(dim_2);
		rRightHandSideVector.resize(dim_2);
		for (unsigned int i = 0; i <geom.size(); i++)
		  {
		      array_1d<double,3>& Vel = geom[i].FastGetSolutionStepValue(VELOCITY); 
		      Velocities[index]       = Vel[0];
		      Velocities[index+1]     = Vel[1];
		      if(dim==2){
		           index = index + 2;
		           }
		      else{
			  Velocities[index+2]  = Vel[2];
			  index = index + 3;
		        }
		  }
		  
		  noalias(rRightHandSideVector) = -malpha_damp * prod(rLeftHandSideMatrix, Velocities);
		  for (unsigned int i = 0; i <geom.size(); i++){
		      array_1d<double,3>& rhs = geom[i].FastGetSolutionStepValue(RHS); 
		      geom[i].SetLock();
		      index  = i*dim;
		      rhs[0] = rhs[0] + rRightHandSideVector[index];
		      rhs[1] = rhs[1] + rRightHandSideVector[index+1];
		      rhs[2] = 0.00; 
		      if(dim==3)
			rhs[2]=  rhs[2] + rRightHandSideVector[index+2];
		      geom[i].UnSetLock();
		  }
		  
	      }
	}
      }
	    KRATOS_CATCH("")
  
}


/// for discontinum Galerking methods  
void ComputeInterfaceForces()
{
    KRATOS_TRY
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();  
    ElementsArrayType& pElements     = r_model_part.Elements(); 

    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> element_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);

    Vector rRightHandSideVector;
    
    
    #pragma omp parallel for private(rRightHandSideVector)
    for(int k=0; k<number_of_threads; k++)
    {
    typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
    typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
    for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
     {
        KRATOS_WATCH(it->Id())
        //it->Calculate(INTERFACE_FORCES, rRightHandSideVector, CurrentProcessInfo);
     }
    }
    
    KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAA")
    KRATOS_CATCH("")
}
 
  
  
/// Computa las fuerzas viscosas de amortiguamiento alfa * K_elem * Velocidad
void ComputeDampingForcesWithStiffness()
{
  
      KRATOS_TRY
      ModelPart& r_model_part          = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();  
      ElementsArrayType& pElements     = r_model_part.Elements(); 

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
       #endif

      vector<unsigned int> element_partition;
      CreatePartition(number_of_threads, pElements.size(), element_partition);

      unsigned int index = 0; 
      unsigned int dim_2 = 0; 
      Matrix rLeftHandSideMatrix;
      Vector rRightHandSideVector;
      Vector Velocities;
      
      #pragma omp parallel for private(index, dim_2, rLeftHandSideMatrix,Velocities, rRightHandSideVector) 
      for(int k=0; k<number_of_threads; k++)
      {
	typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
	typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
	for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	    /// Las fuerzas de contacto a los elementos que estan en contacto
	    if((it)->GetValue(IS_TARGET)==true)
	    {
		it->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, CurrentProcessInfo);
		Element::GeometryType& geom         = it->GetGeometry();
		const unsigned int& dim             = geom.WorkingSpaceDimension();
		const unsigned int& number_of_nodes = geom.size();
		index                               = 0;
		dim_2                               =  dim * number_of_nodes;
		Velocities.resize(dim_2);
		rRightHandSideVector.resize(dim_2);
		for (unsigned int i = 0; i <geom.size(); i++)
		  {
		      array_1d<double,3>& Vel = geom[i].FastGetSolutionStepValue(VELOCITY); 
		      Velocities[index]       = Vel[0];
		      Velocities[index+1]     = Vel[1];
		      if(dim==2){
		           index = index + 2;
		           }
		      else{
			  Velocities[index+2]  = Vel[2];
			  index = index + 3;
		        }
		  }
		  
		  noalias(rRightHandSideVector) = -mbeta_damp * prod(rLeftHandSideMatrix, Velocities);
		  for (unsigned int i = 0; i <geom.size(); i++){
		      array_1d<double,3>& rhs = geom[i].FastGetSolutionStepValue(RHS); 
		      geom[i].SetLock();
		      index  = i*dim;
		      rhs[0] = rhs[0] + rRightHandSideVector[index];
		      rhs[1] = rhs[1] + rRightHandSideVector[index+1];
		      rhs[2] = 0.00; 
		      if(dim==3)
			rhs[2]=  rhs[2] + rRightHandSideVector[index+2];
		      geom[i].UnSetLock();
		  }
		  
	      }
	}
      }
	    KRATOS_CATCH("")
  
}

// Calcula fuerza de contacto
void ComputeNormalContactReactions(const unsigned int& initial_conditions_size)
{
      KRATOS_TRY     
      ModelPart& r_model_part                    = BaseType::GetModelPart();
      ConditionsContainerType& pConditions       = r_model_part.ConditionsArray();
      ConditionsContainerIterator  end_previos   = pConditions.begin() + initial_conditions_size; 
      ConditionsContainerIterator  end_actual    = pConditions.end();
      ProcessInfo& CurrentProcessInfo            = r_model_part.GetProcessInfo();
      array_1d<double,3> Output;
      
      for(ConditionsContainerIterator rcond = end_previos; rcond!=end_actual; rcond++)
	  (*rcond)->Calculate(NORMAL, Output, CurrentProcessInfo);
      
      
      KRATOS_CATCH("") 
}


// Borra las condiciones de contacto anteriores sin borrar las otras condiciones de las estructura
void CheckConditionsStatus(const unsigned int& initial_conditions_size)
{
      KRATOS_TRY    
  
      ModelPart& r_model_part                    = BaseType::GetModelPart();
      ConditionsContainerType& pConditions       = r_model_part.ConditionsArray();
      ConditionsContainerIterator  end_previos   = pConditions.begin() + initial_conditions_size; 
      ConditionsContainerIterator  end_actual    = pConditions.end();
      pConditions.erase(end_previos, end_actual);
      
      KRATOS_CATCH("")
}


/// me permite los iteradores de los links contacts
int CheckObjectContact(const unsigned int& initial_conditions_size,
		       ConditionsContainerIterator& end_previos,   
                       ConditionsContainerIterator& end_actual)
{
  
      KRATOS_TRY
      

      ModelPart& r_model_part                     = BaseType::GetModelPart();
      ConditionsContainerType& pConditions        = r_model_part.ConditionsArray();
      end_previos  = pConditions.begin() + initial_conditions_size;  
      end_actual   = pConditions.end();      
      return (std::distance(end_previos, end_actual)); 
      
     KRATOS_CATCH("")
}

/// comuta el D(n+1) a partir de V(n+1/2) 
void ComputeIntermedialVelocityAndNewDisplacement()
{
  
    KRATOS_TRY
    
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
	  /// Nota = Se ha agregado a la ecuacion la velocidad. Antes no habia  actual_velocity;   
	  //array_1d<double,3>& actual_velocity       = i->FastGetSolutionStepValue(VELOCITY);        /// Estamos en paso  T(n+1)
	  array_1d<double,3>& actual_displacement   = i->FastGetSolutionStepValue(DISPLACEMENT);     /// Estamos en paso  T(n+1)
	  array_1d<double,3>& current_displacement  = i->FastGetSolutionStepValue(DISPLACEMENT,1);   /// U(n)
	  array_1d<double,3>& old_displacement      = i->FastGetSolutionStepValue(DISPLACEMENT,2);   /// U(n-1)
	  array_1d<double,3>& current_Rhs           = i->FastGetSolutionStepValue(RHS);              /// Fext(n) - Fint(n)  	  
	  
	  const double& nodal_mass    =  i->FastGetSolutionStepValue(NODAL_MASS); 
	  const double nodal_damping  =  malpha_damp * nodal_mass; 
	  const double factor_a       =  2.00 *  nodal_mass + nodal_damping *  current_delta_time;    
	  const double factor_b       =  2.00 *  nodal_mass - nodal_damping *  molddelta_time;  
	  const double factor_c       =  2.00 *  mid_delta_time;
	  
	  noalias(mid_neg_velocity)   = (1.00/molddelta_time) * (current_displacement - old_displacement);  
	  noalias(mid_pos_velocity)   = (factor_b/factor_a) * mid_neg_velocity + (factor_c/factor_a) * current_Rhs; //+ actual_velocity;
  
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
    
    //std::cout<< "  OLD TIMESTEP                 = "<< molddelta_time     << "  SECONDS" << std::endl; 
    //std::cout<< "  CURRENT TIMESTEP             = "<< current_delta_time << "  SECONDS" << std::endl; 
    //std::cout<< "  AVERAGE TIMESTEP             = "<< mid_delta_time     << "  SECONDS" << std::endl; 
      
    KRATOS_CATCH("")
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
	    (it)->MassMatrix(MassMatrix, CurrentProcessInfo);
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
     //r_model_part.GetCommunicator().AssembleCurrentData(NODAL_MASS); 
     mElementsAreInitialized   = true;
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
    
      
      mConditionsAreInitialized = true;
      
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

void ComputeCriticalTime()
{
    KRATOS_TRY
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();

    //double delta_time_penalty  = 1.00;
    double delta_time_computed = 1.00;
    double delta_time_used     = 1.00; 
    double time                = 0.00;
    const int step             = CurrentProcessInfo[TIME_STEPS];  

    
    delta_time_computed = ComputeTime();
    
    //Calculo los factores de alfa y beta para amortiguar la estructura
    // 3d static and dynamic analisis damping and Energy Disispation 19-7
   if(mdamping_coeficients==false){
      double wmax = 2.00 / delta_time_computed;
      mbeta_damp  = mdamping_ratio / wmax ; 
      malpha_damp = mdamping_ratio * wmax ; 
      mdamping_coeficients = true;
   }
    
    
    if(mCE==Penalty_Methods){
      double time_penalty =  ComputeTimePenalty();
      delta_time_computed = (delta_time_computed<time_penalty) ? delta_time_computed:time_penalty;
    }
    
    delta_time_computed = Truncar_Delta_Time(delta_time_computed);
    
    
    if(delta_time_computed>mmax_delta_time) {delta_time_computed = mmax_delta_time;}
    delta_time_used = mfraction_delta_time * delta_time_computed;
    CurrentProcessInfo[DELTA_TIME] = delta_time_used; 
    
    mtimestep = mtimestep + delta_time_used; 
    r_model_part.CloneTimeStep(mtimestep);
    time = CurrentProcessInfo[TIME];
    
    
    if(mCalculateOldTime==false) {
        molddelta_time    = delta_time_used;
	mCalculateOldTime = true;
    }
    
    std::cout<< "  BETA_DAMPING FOR STIFFNESS   = "<< mbeta_damp             << "         " << std::endl;
    std::cout<< "  ALPHA_DAMPING FOR MASS       = "<< malpha_damp            << "         " << std::endl;
    std::cout<< "  TIME STEPS                   = "<< step                   << "         " << std::endl;
    std::cout<< "  FACTOR DELTA CRITICAL TIME   = "<< mfraction_delta_time   << "         "  << std::endl;
    std::cout<< "  DELTA CRITICAL TIME COMPUTED = "<< delta_time_computed    << "  SECONDS" << std::endl; 
    std::cout<< "  DELTA TIME USED              = "<< delta_time_used        << "  SECONDS" << std::endl;
    std::cout<< "  CURRENT TIME                 = "<< time                   << "  SECONDS" << std::endl;  
    KRATOS_CATCH("")

}

//***************************************************************************
//***************************************************************************

double ComputeTime()
{
  
    KRATOS_TRY
  
    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    ElementsArrayType& pElements     = r_model_part.Elements();     
    //double step                    =  CurrentProcessInfo[TIME_STEPS];  

    #ifdef _OPENMP
       int number_of_threads = omp_get_max_threads();
       //thread_id         = omp_get_thread_num();
    #else
       int number_of_threads = 1;
       //thread_id = 0;
    #endif
    
    vector<unsigned int> element_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);

    Vector dts(number_of_threads);
    double delta_time_a = 0.00;
    for(int i = 0; i < number_of_threads; i++)
       dts[i] = mmax_delta_time;
    
    
    // WARNING = Los threads Id siempre empiezan por cero.
    #pragma omp parallel for private(delta_time_a)  //shared(CurrentProcessInfo)
    for(int k=0; k<number_of_threads; k++){
    typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
    typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
    for(ElementsArrayType::iterator it=it_begin; it!= it_end; it++){
	  it-> Calculate(DELTA_TIME, delta_time_a, CurrentProcessInfo);
	  if(delta_time_a>0.00) 
	     if(delta_time_a < dts[k])
	        dts[k] = delta_time_a;
       } }
       
    
   return  (*std::min_element(dts.begin(), dts.end())); 
    
   KRATOS_CATCH("")
}


///WARNING = Cuidado con los nodos que sueltos. Aquellos que no tienen masa
double ComputeTimePenalty()
{
      ModelPart& r_model_part         = BaseType::GetModelPart();
      //ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
      NodesArrayType& pNodes          = r_model_part.Nodes(); 
      //ElementsArrayType& pElements    = r_model_part.Elements(); 
      
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      Vector Min_Mass_Nodal(number_of_threads); 
      vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);
      
      Properties&  prop  = r_model_part.GetProperties(1);
      double aux         = 0.00;
      double Penalty     = (prop)[YOUNG_MODULUS];
      double Mass        = (r_model_part.Nodes()(1))->FastGetSolutionStepValue(NODAL_MASS); 
      std::size_t nprop  = r_model_part.NumberOfProperties();
            
      for(std::size_t i = 1; i<= nprop; i++)
      {
	 Properties&  prop_aux  = r_model_part.GetProperties(i);
	 aux                    = (prop_aux)[YOUNG_MODULUS];
	 Penalty                = ( Penalty > aux ) ? aux : Penalty; 
      }

      Penalty = mpenalty_factor * Penalty;
      
      for(int k=0; k<number_of_threads; k++)
        Min_Mass_Nodal[k] = Mass; 
	
      #pragma omp parallel for 
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
	  {
	    
	    double& mass = ((i)->FastGetSolutionStepValue(NODAL_MASS)); 
	    if(mass==0.00) KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAAA")
	    Min_Mass_Nodal[k] =  (Min_Mass_Nodal[k] < mass) ? Min_Mass_Nodal[k] : mass;
	  }
	}
	
	Mass = (*std::min_element(Min_Mass_Nodal.begin(), Min_Mass_Nodal.end()));
	double result = 0.50 * std::sqrt(2.00 * Mass / Penalty); 
	return result;
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
	    const array_1d<double,3>& CurrentAcceleration  =  (i->FastGetSolutionStepValue(ACCELERATION)) + (i->FastGetSolutionStepValue(GRAVITY)); 
	    array_1d<double,3>&       OldDisplacement      =  (i->FastGetSolutionStepValue(DISPLACEMENT,2));
	  
	    /// D_1 7.145 Libro de ""Estructuras Sometiadas a Acciones Sismicas"   
	    noalias(OldDisplacement) =  0.5*DeltaTime*DeltaTime*CurrentAcceleration - DeltaTime*CurrentVelocity + CurrentDisplacement;  
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
      
      /// Set to zero de RHS
      SetToZeroRHS();
      
      /// Compute the global external nodal force.
      Calculate_Conditions_RHS_and_Add();
      
      /// Compute the stress and body force of the element.
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
            //array_1d<double,3>& reaction  = (i->FastGetSolutionStepValue(REACTION));
	    array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(RHS)); 
	    
	    noalias(normal)   = ZeroVector(3);
	    //noalias(reaction) = ZeroVector(3);
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
	       {  geom(i)->SetLock();
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

      Vector rhs_elem;

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> element_partition;
      CreatePartition(number_of_threads, pElements.size(), element_partition);

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
		    { geom(i)->SetLock();      
		      node_rhs[kk] = node_rhs[kk] + rhs_elem[index+kk];
                      geom(i)->UnSetLock();
		    }

	      }
	}
    }
    
     KRATOS_CATCH("")
}


//***************************************************************************
//***************************************************************************

/*
void Calculate_Final_Force_Contribution(ModelPart::NodeIterator& i)  
{
        KRATOS_TRY
        array_1d<double,3> Final_Force; 
        const double nodal_damping         =  malpha_damp*((i)->FastGetSolutionStepValue(NODAL_MASS)); 
        array_1d<double,3>& node_rhs       =  i->FastGetSolutionStepValue(RHS);    
	const array_1d<double,3>& velocity =  i->FastGetSolutionStepValue(VELOCITY);
	noalias(Final_Force)               =  nodal_damping * velocity;    
	for(unsigned int i=0; i<3; i++)
 	        node_rhs[i] = node_rhs[i]-Final_Force[i];  
	
	KRATOS_CATCH("")
}

*/

//***************************************************************************
//***************************************************************************

/// Encuentra las velocidades y aceleraciones del paso n.
void ComputeOldVelocitiesAndAccelerations()
{ 
  
    KRATOS_TRY
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
	  array_1d<double,3>& normal                = i->FastGetSolutionStepValue(NORMAL);
	  noalias(normal)                           = ZeroVector(3);
	  array_1d<double,3>& displacement          = i->FastGetSolutionStepValue(DISPLACEMENT);   /// Estamos en paso  T(n+1)
	  array_1d<double,3>& current_displacement  = i->FastGetSolutionStepValue(DISPLACEMENT,1);   /// U(n)
	  array_1d<double,3>& olddisplacement       = i->FastGetSolutionStepValue(DISPLACEMENT,2);   /// U(n-1)
	  array_1d<double,3>& veryolddisplacement   = i->FastGetSolutionStepValue(DISPLACEMENT,3);   /// U(n-1)
	  
	  noalias(mid_neg_velocity)       = (1.00/molddelta_time)     * (current_displacement   - olddisplacement);
	  noalias(mid_pos_velocity)       = (1.00/current_delta_time) * (displacement           - current_displacement);
	  noalias(mid_neg_velocity_old)   = (1.00/molddelta_time)     * (olddisplacement        - veryolddisplacement);
	  noalias(mid_pos_velocity_old)   = (1.00/current_delta_time) * (current_displacement   - olddisplacement);

            
	    array_1d<double,3>& velocity     = (i->FastGetSolutionStepValue(VELOCITY));  
 	    array_1d<double,3>& acceleration = (i->FastGetSolutionStepValue(ACCELERATION));
	    //KRATOS_WATCH(velocity)
	    
	    /// X 
	    if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
	    { 
	       //if(velocity[0]==0.00)
	       {
	          velocity[0]      =  inv_sum_delta_time * (displacement[0] - olddisplacement[0]);  
	          acceleration[0]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[0] - mid_neg_velocity[0]);  
	        }
	   }
            else
	    {
	       //velocity[0]      =  inv_sum_delta_time * (displacement[0] - olddisplacement[0]);
	       //acceleration[0]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[0] - mid_neg_velocity[0]);  
	       velocity[0]         =  inv_sum_delta_time * (current_displacement[0] - veryolddisplacement[0]);
	       acceleration[0]     =  (1.00 / mid_delta_time ) * (mid_pos_velocity_old[0] - mid_neg_velocity_old[0]);  
	    }
            

	    
 	    if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
 	    {
	       //if(velocity[1]==0.00)
	       {
	       velocity[1]      =  inv_sum_delta_time * (displacement[1] - olddisplacement[1]);
	       acceleration[1]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[1] - mid_neg_velocity[1]);  
	       }
 	    }
             else
 	    {
	       //velocity[1]      =  inv_sum_delta_time * (displacement[1] - olddisplacement[1]);
	       //acceleration[1]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[1] - mid_neg_velocity[1]);
 	       velocity[1]         =  inv_sum_delta_time * (current_displacement[1] - veryolddisplacement[1]);
 	       acceleration[1]     =  (1.00 / mid_delta_time ) * (mid_pos_velocity_old[1] - mid_neg_velocity_old[1]);  
      
 	    }

            /// Z 
	    if( i->HasDofFor(DISPLACEMENT_Z)){
 	    if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
 	    {
	      //if(velocity[2]==0.00)
	      {
	      velocity[2]      =  inv_sum_delta_time * (displacement[2] - olddisplacement[2]);
	      acceleration[2]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[2] - mid_neg_velocity[2]);  
	      }
 	    }
 	    else
 	    {
	       //velocity[2]      =  inv_sum_delta_time * (displacement[2] - olddisplacement[2]);
	       //acceleration[2]  =  (1.00 / mid_delta_time ) * (mid_pos_velocity[2] - mid_neg_velocity[2]);  
 	       velocity[2]         =  inv_sum_delta_time * (current_displacement[2] - veryolddisplacement[2]);
 	       acceleration[2]     =  (1.00 / mid_delta_time ) * (mid_pos_velocity_old[2] - mid_neg_velocity_old[2]);  
 	      
 	    }  
	  }  
	  //KRATOS_WATCH(velocity)
	  //KRATOS_WATCH("-----------------------")
	  
// 	  if(i->Id()==374)
// 	     KRATOS_WATCH(velocity) 
      }    
  }
  
  KRATOS_CATCH("")
}
 

 
/// Computes the contact force and displacement corrections 
void Update(const ConditionsContainerIterator& end_previos,   
            const ConditionsContainerIterator& end_actual)
{
   KRATOS_TRY
   mpLagrangianMultiplier->CalculateContactForceAndDisplacementCorrections(end_previos, end_actual); 
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
	    (it) -> GetValue(IS_TARGET)=false;
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



void ChangeContactConditions(const bool& contact)
{ 
  KRATOS_TRY
    mComputeContactConditions = contact;
  KRATOS_CATCH("")
}

// Crea nuevamente las condiciones de contacto.
// WARNING = Hay que evitar que borre las existentes 
void RecalculateBoundaryContours(const bool& rflag)
{
   KRATOS_TRY 
     std::cout<<"RECOMPUTING BOUNDARY CONTOURS " << std::endl;
     mpBCCU_Pointer->ResetFlagComputeBoundaryContour(rflag);
   KRATOS_CATCH("")
}

void ChangeFractionDeltaTime(const double& new_fraction_delta_time)
{
   KRATOS_TRY
   mfraction_delta_time = new_fraction_delta_time;
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
	    array_1d<double,3>& reaction           = (i->FastGetSolutionStepValue(REACTION));
            const double& mass                           = (i->FastGetSolutionStepValue(NODAL_MASS));
            const array_1d<double,3>& rhs                = (i->FastGetSolutionStepValue(RHS));
            const array_1d<double,3>& acceleration       = (i->FastGetSolutionStepValue(ACCELERATION));
            const array_1d<double,3>& velocity           = (i->FastGetSolutionStepValue(VELOCITY));
            //array_1d<double,3> dif                 =  rhs; 

            if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == true)
 		{reaction[0] =  -rhs[0] + mass * acceleration[0] + mass * malpha_damp * velocity[0];}
                
            if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == true )
		{reaction[1] = -rhs[1] + mass * acceleration[1] + mass * malpha_damp * velocity[1];}
          
            if( i->HasDofFor(DISPLACEMENT_Z))
	        {
		if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == true )
		{reaction[2] = -rhs[2] + mass * acceleration[2] + mass * malpha_damp * velocity[2];}
                }
            }
      } 
 }


private:


unsigned int    mdimension;
unsigned int    minitial_conditions_size;  
unsigned int    mcontact_conditions_size;  
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
bool   mClearContactConditions;
bool   mLocalSearchInitialize;
bool   mdamping_coeficients;
double mdamping_ratio;
double malpha_damp;
double mbeta_damp; 
double mfraction_delta_time;
double mmax_delta_time;
double molddelta_time;
double mtimestep;  /// la suma de los delta time
double mpenalty_factor;
Constraint_Enforcement mCE; 
Disconnect_Triangle_Utilities mDTU;


/// Contact Condition Iterations
//typename ProofType::Pointer mProof;
typename BoundaryAndContactType::Pointer mpBCCU_Pointer;
typename LagrangeMultiplierType::Pointer mpLagrangianMultiplier;



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


void CalculateEnergies()
{
  CalculatePotecialEnergy();
  CalculateKineticEnergy(); 
  //CalculateDeformationEnergy();
  
}

void CalculateDeformationEnergy()
{
}


/// Asumimos que la gravedad esta en direccion Y 
void CalculatePotecialEnergy()
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
    
    double h_efe;
    #pragma omp parallel for  private(h_efe)
    for(int k=0; k<number_of_threads; k++)
      {
	typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	 { 
	    if(i->Y()>0.00){
	    double& y_coord              = i->Y0();
	    const double& displ          = i->FastGetSolutionStepValue(DISPLACEMENT_Y);  
	    const double& gravity        = i->FastGetSolutionStepValue(GRAVITY_Y); 
            const double& nodal_mass     = i->FastGetSolutionStepValue(NODAL_MASS);
	    double& potencial_energy     = i->FastGetSolutionStepValue(POTENCIAL_ENERGY); 
	    potencial_energy             = 0.00;
	    h_efe                        = y_coord + displ;
	    potencial_energy             = std::fabs(nodal_mass * gravity * h_efe );
	    } 
	  }
        }
  KRATOS_CATCH("")
}
 

void CalculateKineticEnergy()
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
    double vel = 0.00; 
    //double h   = 0.1;
    #pragma omp parallel for  private(vel) 
    for(int k=0; k<number_of_threads; k++)
      {
	typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	 { 
	    array_1d<double,3>& Velocity = i->FastGetSolutionStepValue(VELOCITY); 
	    double& mass                 = i->FastGetSolutionStepValue(NODAL_MASS);  
	    double& kinetic_energy       = i->FastGetSolutionStepValue(KINETIC_ENERGY);
	    vel                          = norm_2(Velocity); 
	    /// WARNING = solo valido para matriz de masa diagonal
	    /// Ek = 0.50 * V^t * M * V 
	    kinetic_energy               = 0.50 * mass * vel * vel; 
	 }
      }
  KRATOS_CATCH("")
}

void Heuristic_Formula(std::vector<Joint2D>::iterator Begin, std::vector<Joint2D>::iterator End)
        {
	 KRATOS_TRY
	  double dpefa = 0.63;
	  double dpefb = 1.8;
	  double dpefc = 6.0;
	  double dpefm = 0.0;
	  double small,sabs,o,s,o1,o2,s1,s2,op,sp,ot,st,z,sigma,tau;
	  double e1x,e1y,h,area;
	  int nfail;
	  //int nsoft;
          double d1nccx[4];
	  double d1nccy[4];
	  
	  ModelPart& r_model_part = BaseType::GetModelPart();
	  Properties&  prop       = r_model_part.GetProperties(1);
	  
	  
	  const double& dpeft = (prop)[TENSILE_STRENGTH];
	  const double dpepe  =  mpenalty_factor * (prop)[YOUNG_MODULUS];
	  const double& dpefs = (prop)[SHEAR_STRENGTH];
	  const double& dpegf = (prop)[FRACTURE_ENERGY];
	  
	  small=1E-9; //nsoft=0;
	  for(std::vector<Joint2D>::iterator Joint = Begin; Joint != End; ++Joint) 
	    { 
	      if((Joint)->IsFail()==false){
	      array_1d<double,3>& node_rhs_0 =  (*Joint)[0]->FastGetSolutionStepValue(RHS);
	      array_1d<double,3>& node_rhs_1 =  (*Joint)[1]->FastGetSolutionStepValue(RHS);
	      array_1d<double,3>& node_rhs_2 =  (*Joint)[2]->FastGetSolutionStepValue(RHS);
	      array_1d<double,3>& node_rhs_3 =  (*Joint)[3]->FastGetSolutionStepValue(RHS);
	      
	      d1nccx[0] = (*Joint)[0]->X();  
	      d1nccx[1] = (*Joint)[1]->X();  
	      d1nccx[2] = (*Joint)[2]->X();  
	      d1nccx[3] = (*Joint)[3]->X();  
	      
	      d1nccy[0] = (*Joint)[0]->Y();  
	      d1nccy[1] = (*Joint)[1]->Y();  
	      d1nccy[2] = (*Joint)[2]->Y();  
	      d1nccy[3] = (*Joint)[3]->Y();  
	      
	      e1x=0.50*(d1nccx[1]+d1nccx[2]-d1nccx[0]-d1nccx[3]);
	      e1y=0.50*(d1nccy[1]+d1nccy[2]-d1nccy[0]-d1nccy[3]);
	      h=std::sqrt(e1x*e1x+e1y*e1y);
	      
	      e1x=e1x/(h+small);
	      e1y=e1y/(h+small);
	      s1=(d1nccy[0]-d1nccy[3])*e1y+(d1nccx[0]-d1nccx[3])*e1x;
	      s2=(d1nccy[1]-d1nccy[2])*e1y+(d1nccx[1]-d1nccx[2])*e1x;
	      o1=(d1nccy[0]-d1nccy[3])*e1x-(d1nccx[0]-d1nccx[3])*e1y;
	      o2=(d1nccy[1]-d1nccy[2])*e1x-(d1nccx[1]-d1nccx[2])*e1y;
	      
	      
	      op=2.00*h*dpeft/dpepe;
	      sp=2.00*h*dpefs/dpepe;
	      ot=std::max((2.00*op),(3.00*dpegf/dpeft));
	      st=std::max((2.00*sp),(3.00*dpegf/dpefs));
	      //nfail=0;
	       
	      for(int integ=0;integ<3;integ++)
	      {  if(integ==0)
	         { o=o1; s=s1;
	         }
	         else if(integ==2)
	         { o=o2; s=s2;
	         } 
	         else
	         { o=0.50*(o1+o2); s=0.50*(s1+s2);
	         }
	         
	         sabs=std::fabs(s);
	         if((o>op)&&(sabs>sp))
	         { z=std::sqrt(((o-op)/ot)*((o-op)/ot)+((sabs-sp)/st)*((sabs-sp)/st));
	         }
	         else if(o>op)
	         { z=(o-op)/ot;
	         }
	         else if(sabs>sp)
	         { z=(sabs-sp)/st;
	         }
	         else
	         { z=0.00;
	         }
	         
	         if(z>=1.00)    
	         {
		   nfail=nfail+1;
		   if((nfail>1))
		      (Joint)->SetFail();
	            z=1.00;
	         }
	         
	         z=(1.00 - ((dpefa+dpefb-1.00)/(dpefa+dpefb))*exp(z*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(1.00-dpefa-dpefb))))*(dpefa*(1.00-z)+dpefb*pow((1.00-z),dpefc));        
	        

		if(o<0.00)                 /* normal stress*/ 
	        { sigma=2.00*o*dpeft/op;   /* sigma=R0; */
	        }
	        else if(o>op)
	        { sigma=dpeft*z; //nsoft=nsoft+1;
	        }
	        else
	        { sigma=(2.00*o/op-(o/op)*(o/op))*z*dpeft;
	        }
	        if((sigma>0.00)&&(sabs>sp))           /* shear stress */ 
	        { tau=z*dpefs;
	        }
	        else if(sigma>0.00)
	        { tau=(2.00*(sabs/sp)-(sabs/sp)*(sabs/sp))*z*dpefs;
	        }
	        else if(sabs>sp)
	        { tau=z*dpefs-dpefm*sigma;
	        }
	        else
	        { tau=(2.00*(sabs/sp)-(sabs/sp)*(sabs/sp))*(z*dpefs-dpefm*sigma);
	        }  
		if(s<0.00)tau=-tau;
		if(integ==0)  /* nodal forces */
		{ 
		  area=h/6.00; /* area=h/6.0; */
		  node_rhs_0[0] = node_rhs_0[0] - area*(tau*e1x-sigma*e1y); 
		  node_rhs_0[1] = node_rhs_0[1] - area*(tau*e1y+sigma*e1x);
		  node_rhs_3[0] = node_rhs_3[0] + area*(tau*e1x-sigma*e1y);
		  node_rhs_3[1] = node_rhs_3[1] + area*(tau*e1y+sigma*e1x);
		}
		else if(integ==1)
		{ 
		  area=h/3.00;  /* area=h/3.0; */
		  node_rhs_0[0] = node_rhs_0[0]-area*(tau*e1x-sigma*e1y);
		  node_rhs_0[1] = node_rhs_0[1]-area*(tau*e1y+sigma*e1x);
		  node_rhs_3[1] = node_rhs_3[1]+area*(tau*e1y+sigma*e1x);
		  node_rhs_3[0] = node_rhs_3[0]+area*(tau*e1x-sigma*e1y);
		  node_rhs_1[0] = node_rhs_1[0]-area*(tau*e1x-sigma*e1y);
		  node_rhs_1[1] = node_rhs_1[1]-area*(tau*e1y+sigma*e1x);
		  node_rhs_2[1] = node_rhs_2[1]+area*(tau*e1y+sigma*e1x);
		  node_rhs_2[0] = node_rhs_2[0]+area*(tau*e1x-sigma*e1y);
		}
		else
		{ 
		  area=h/6.00; /* area=h/6.0; */
		  node_rhs_1[0]=node_rhs_1[0]-area*(tau*e1x-sigma*e1y);
		  node_rhs_1[1]=node_rhs_1[1]-area*(tau*e1y+sigma*e1x);
		  node_rhs_2[1]=node_rhs_2[1]+area*(tau*e1y+sigma*e1x); 
		  node_rhs_2[0]=node_rhs_2[0]+area*(tau*e1x-sigma*e1y);
	      } } } } 
	  KRATOS_CATCH("")
	}

};  

} /* namespace Kratos.*/
#endif /* KRATOS_RESIDUALBASED_CENTRAL_DIFERENCES_STRATEGY */


