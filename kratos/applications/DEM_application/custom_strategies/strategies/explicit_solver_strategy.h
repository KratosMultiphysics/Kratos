//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//

#if !defined(KRATOS_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_EXPLICIT_SOLVER_STRATEGY

// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"

#include "custom_elements/Particle_Contact_Element.h"
#include "includes/variables.h"
#include "DEM_application.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <time.h> 

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/integration_scheme.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_utilities/dem_fem_utilities.h"

////Cfeng
#include "custom_utilities/dem_fem_search.h"

/* Timer defines */
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /*bool compare_ids(const Element::WeakPointer& i, const Element::WeakPointer& j) {
    
     if( (i.lock())->Id() <  (j.lock())->Id()) return true;
     if( (i.lock())->Id() >  (j.lock())->Id()) return false;
  
     std::cout<<"Two elements with the same Id!! problems can occur here!!"<<std::endl<<std::flush;           
     return false;               
     
         
  }  */
    
  /// Short class definition.
  /** Detail class definition.
  */
  template<
  class TSparseSpace,
  class TDenseSpace,
  class TLinearSolver>
  class ExplicitSolverStrategy: public  SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {
      public:
      ///@name Type Definitions
      ///@{

      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>   BaseType;
      typedef ModelPart::NodesContainerType                             NodesArrayType;
      typedef ModelPart::ElementsContainerType                          ElementsArrayType;
      typedef ElementsArrayType::iterator                               ElementsIterator;

      typedef ModelPart::ConditionsContainerType                        ConditionsArrayType;

      typedef ModelPart::NodesContainerType::ContainerType              NodesContainerType;
      typedef ModelPart::ElementsContainerType::ContainerType           ElementsContainerType;
      typedef ModelPart::ConditionsContainerType::ContainerType         ConditionsContainerType;

      typedef SpatialSearch::ResultElementsContainerType                ResultElementsContainerType;
      typedef SpatialSearch::VectorResultElementsContainerType          VectorResultElementsContainerType;

      typedef SpatialSearch::RadiusArrayType                            RadiusArrayType;
      typedef SpatialSearch::DistanceType                               DistanceType;
      typedef SpatialSearch::VectorDistanceType                         VectorDistanceType;
      
      
      //Cfeng
      typedef SpatialSearch::ResultConditionsContainerType              ResultConditionsContainerType;
      typedef SpatialSearch::VectorResultConditionsContainerType        VectorResultConditionsContainerType;
 
      typedef PointerVectorSet<Properties, IndexedObject>               PropertiesContainerType;
      typedef typename PropertiesContainerType::iterator                PropertiesIterator;
 

      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      ExplicitSolverStrategy(){}

      ExplicitSolverStrategy(ModelPart& r_model_part,
                         ModelPart& fem_model_part,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool move_mesh_flag,
                             const int delta_option,
                             const double search_tolerance,
                             const double coordination_number,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename IntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch
      ): BaseType(r_model_part, move_mesh_flag)
      {

          mElementsAreInitialized        = false;
          mInitializeWasPerformed        = false;
          mDeltaOption                   = delta_option;
          mSearchTolerance               = search_tolerance;
          mCoordinationNumber            = coordination_number;
          mpParticleCreatorDestructor    = p_creator_destructor;
          mpScheme                       = pScheme;
          mpSpSearch                     = pSpSearch;
          mMaxTimeStep                   = max_delta_time;
          mNStepSearch                   = n_step_search;
          mSafetyFactor                  = safety_factor;
          mNumberOfElementsOldRadiusList = 0;
          
          mpDem_model_part               = &r_model_part;
          mpFem_model_part               = &fem_model_part;
          
      }

      /// Destructor.
      virtual ~ExplicitSolverStrategy()
      {
          Timer::SetOuputFile("TimesPartialRelease");
          Timer::PrintTimingInformation();
      }            
      
      
      void CreatePropertiesProxies(){
          KRATOS_TRY
                  
          ModelPart& r_model_part             = BaseType::GetModelPart();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();
          
          
          int number_of_properties = r_model_part.NumberOfProperties();
          mFastProperties.resize(number_of_properties);
          
          PropertiesProxy aux_props;
          int i = 0;
          
          
          for (PropertiesIterator props_it = r_model_part.GetMesh(0).PropertiesBegin(); 
                                  props_it!= r_model_part.GetMesh(0).PropertiesEnd();   props_it++ )
          {
              aux_props.SetId( props_it->GetId() );

              double* aux_pointer = &( props_it->GetValue(YOUNG_MODULUS) );
              aux_props.SetYoungFromProperties( aux_pointer );
              
              aux_pointer = &( props_it->GetValue(POISSON_RATIO) );
              aux_props.SetPoissonFromProperties(aux_pointer);
              
              if ( r_model_part.GetProcessInfo()[ROLLING_FRICTION_OPTION] )  {
                aux_pointer = &( props_it->GetValue(ROLLING_FRICTION) );
                aux_props.SetRollingFrictionFromProperties(aux_pointer);
              }
              else {
                aux_props.SetRollingFrictionFromProperties(NULL);
              }
              
              aux_pointer = &( props_it->GetValue(PARTICLE_FRICTION) );
              aux_props.SetTgOfFrictionAngleFromProperties(aux_pointer);
              
              aux_pointer = &( props_it->GetValue(LN_OF_RESTITUTION_COEFF) );
              aux_props.SetLnOfRestitCoeffFromProperties(aux_pointer);
              
              aux_pointer = &( props_it->GetValue(PARTICLE_DENSITY) );
              aux_props.SetDensityFromProperties(aux_pointer);
              
              mFastProperties[i] = aux_props;
              i++;
              
          }
          
          RebuildPropertiesProxyPointers(pElements);
          ElementsArrayType& pGhostElements        = r_model_part.GetCommunicator().GhostMesh().Elements();
          RebuildPropertiesProxyPointers(pGhostElements);
           
          return;          
          KRATOS_CATCH("")
      }
      void RebuildPropertiesProxyPointers(ElementsArrayType& pElements){
          
          KRATOS_TRY          
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
          
          #pragma omp parallel for
          for (int k = 0; k < this->GetNumberOfThreads(); k++){
              
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (typename ElementsArrayType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it){
                
                Kratos::SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*particle_pointer_it);
                                                  
                int general_properties_id = spheric_particle.GetProperties().Id();  
                for (unsigned int i = 0; i < mFastProperties.size(); i++){
                    int fast_properties_id = mFastProperties[i].GetId(); 
                    if( fast_properties_id == general_properties_id ){                  
                        spheric_particle.SetFastProperties( &(mFastProperties[i]) );
                        break;
                    }
                }                    
            }                                
         }  
          
         return;          
         KRATOS_CATCH("")
      }
      
      virtual void Initialize()
      {
          KRATOS_TRY

          ModelPart& r_model_part            = BaseType::GetModelPart();
          
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();
          
          // Omp initializations
          GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
          mNeighbourCounter.resize(this->GetNumberOfThreads());

          CreatePropertiesProxies();
          
          int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

          GetResults().resize(number_of_elements);
          GetResultsDistances().resize(number_of_elements);          
                    
          GetRigidFaceResults().resize(number_of_elements);
          GetRigidFaceResultsDistances().resize(number_of_elements);

          GetBoundingBoxOption() = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

          InitializeSolutionStep();
          InitializeElements();
          InitializeFEMElements();

          mInitializeWasPerformed = true;
          
          ApplyPrescribedBoundaryConditions();
          
          // 0. Set search radius
          
          SetOriginalRadius(r_model_part);
          SetSearchRadius(r_model_part, 1.0);
          
          SearchNeighbours();
          
          if(mDeltaOption == 2)
          {
            
            SetCoordinationNumber(r_model_part);
            
          }
          
          ComputeNewNeighboursHistoricalData();  
          
          ///Cfeng RigidFace search
          SearchRigidFaceNeighbours();
          ComputeNewRigidFaceNeighboursHistoricalData();
          
          // 3. Finding overlapping of initial configurations

          if (rCurrentProcessInfo[CLEAN_INDENT_OPTION]){
              CalculateInitialMaxIndentations();
          }

          
           // 5. Finalize Solution Step.
          //FinalizeSolutionStep();
          KRATOS_WATCH(r_model_part.GetNodalSolutionStepVariablesList())
                    

      KRATOS_CATCH("")
      }// Initialize()

      virtual double Solve()
      {
          KRATOS_TRY

          ModelPart& r_model_part            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();                    

          int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

          this->GetResults().resize(number_of_elements);
          this->GetResultsDistances().resize(number_of_elements);
          this->GetRadius().resize(number_of_elements);          
          
          int time_step = rCurrentProcessInfo[TIME_STEPS];
                    
          //Cfeng
          this->GetRigidFaceResults().resize(number_of_elements);
          this->GetRigidFaceResultsDistances().resize(number_of_elements);

          // 1. Here we initialize member variables that depend on the rCurrentProcessInfo          
          InitializeSolutionStep();
          
          // 2. Neighbouring search. Every N times. + destruction of particles outside the bounding box                   
    
          if ((time_step + 1) % mNStepSearch == 0 && time_step > 0){
              if (this->GetBoundingBoxOption()){
                  BoundingBoxUtility();
              }
              
              SetSearchRadius(r_model_part, 1.0);              
              SearchNeighbours();
              ComputeNewNeighboursHistoricalData();  
              
              SetOriginalRadius(r_model_part);              
              SearchRigidFaceNeighbours();
              ComputeNewRigidFaceNeighboursHistoricalData();
              
          }
                
          // 3. Get and Calculate the forces
          GetForce();
          //FastGetForce();
          Calculate_Conditions_RHS_and_Add();
          
          // 4. Synchronize (should be just FORCE and TORQUE)
          SynchronizeSolidMesh(r_model_part);
          
          // 5. Motion Integration
          PerformTimeIntegrationOfMotion(rCurrentProcessInfo); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot                  
          
          FinalizeSolutionStep();         

          return 0.00;

          KRATOS_CATCH("")
      }//Solve()

      void InitialTimeStepCalculation()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          ElementsIterator it_begin = pElements.ptr_begin();
          ElementsIterator it_end   = pElements.ptr_end();

          double& process_info_delta_time = rCurrentProcessInfo[DELTA_TIME];
          double temp_time_step           = mMaxTimeStep;
          double elem_critical_time_step  = temp_time_step;

          for (ElementsIterator it = it_begin; it != it_end; it++){
              it->Calculate(DELTA_TIME, elem_critical_time_step, rCurrentProcessInfo);

              if (elem_critical_time_step < temp_time_step){
                  temp_time_step = elem_critical_time_step;
              }

          }

          temp_time_step /= mSafetyFactor;
          process_info_delta_time = temp_time_step;

          std::cout<< "****************** Calculated time step is " << temp_time_step << " ******************" << "\n" << std::endl;

          KRATOS_CATCH("")
      }

      void GetForce()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();
          double dt = rCurrentProcessInfo[DELTA_TIME];
          const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY];

          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

          Vector rhs_elem;

          #pragma omp parallel for private(rhs_elem)
          for (int k = 0; k < this->GetNumberOfThreads(); k++){

              if (rhs_elem.size() != 6){
                rhs_elem.resize(6,false);
              }

              typename ElementsArrayType::iterator it_begin   = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end     = pElements.ptr_begin() + this->GetElementPartition()[k+1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*it);
                  spheric_particle.CalculateRightHandSide(rhs_elem, rCurrentProcessInfo, dt, gravity);                  

              } //loop over particles

          } // loop threads OpenMP

          KRATOS_CATCH("")
      }
      void FastGetForce()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();
          double dt = rCurrentProcessInfo[DELTA_TIME];
          const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY];

          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

          #pragma omp parallel for
          for (int k = 0; k < this->GetNumberOfThreads(); k++){

              typename ElementsArrayType::iterator it_begin   = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end     = pElements.ptr_begin() + this->GetElementPartition()[k+1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*it);
                  spheric_particle.FirstCalculateRightHandSide(rCurrentProcessInfo, dt);                  
              } //loop over particles
              
              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*it);
                  spheric_particle.CollectCalculateRightHandSide(rCurrentProcessInfo);                  
              } //loop over particles
              
              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*it);
                  spheric_particle.FinalCalculateRightHandSide(rCurrentProcessInfo, dt, gravity);                  
              } //loop over particles

          } // loop threads OpenMP

          KRATOS_CATCH("")
      }

      virtual void PerformTimeIntegrationOfMotion(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY
          
          ModelPart& r_model_part = BaseType::GetModelPart();

          GetScheme()->Calculate(r_model_part);

          KRATOS_CATCH("")
      }

      void InitializeSolutionStep()
      {
          KRATOS_TRY

          // SPHERE MODEL PART

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

          #pragma omp parallel for
          for (int k = 0; k < this->GetNumberOfThreads(); k++){
              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                (it)->InitializeSolutionStep(rCurrentProcessInfo); 
                
              } // loop over particles

          } // loop threads OpenMP

        KRATOS_CATCH("")
      }

      virtual void BoundingBoxUtility()
      {
          KRATOS_TRY
          
          ModelPart& r_model_part = BaseType::GetModelPart();
          mpParticleCreatorDestructor->MarkDistantParticlesForErasing(r_model_part);
          mpParticleCreatorDestructor->DestroyParticles(r_model_part);

          KRATOS_CATCH("")
      }

      void MoveMesh(){}

      void FinalizeSolutionStep()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

          #pragma omp parallel for
          for (int k = 0; k < this->GetNumberOfThreads(); k++){
              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  (it)->FinalizeSolutionStep(rCurrentProcessInfo); //we use this function to call the set initial contacts and the add continuum contacts.
              } //loop over particles

          } // loop threads OpenMP

          KRATOS_CATCH("")
      }

      
    void SetCoordinationNumber(ModelPart& r_model_part) 
    {
      double in_coordination_number = mCoordinationNumber;
      double out_coordination_number = ComputeCoordinationNumber();
      int iteration = 0;
      int maxiteration = 100;
 
      if(in_coordination_number <= 0.0)
      {
        KRATOS_ERROR(std::runtime_error,"The specified Coordination Number is less or equal to zero, N.C. = ",in_coordination_number)
      }
      else
      {
          while ( fabs(out_coordination_number/in_coordination_number-1.0) > 1e-3 )
          {              
            if(iteration>=maxiteration)   break;
              
            iteration++;

            mSearchTolerance *= in_coordination_number/out_coordination_number;
            
            SetSearchRadius(r_model_part, 1.0);
            
            SearchNeighbours();
             
            out_coordination_number = ComputeCoordinationNumber();
   
          }//while

          if(iteration<maxiteration) std::cout<< "Coordination Number iteration converged after "<<iteration<< " iterations, to value " <<out_coordination_number<< " using an extension of " << mSearchTolerance <<". "<<"\n"<<std::endl;
          else
            {   
              std::cout<<"Coordination Number iteration did NOT converge after "<<iteration<<" iterations. Coordination number reached is "<<out_coordination_number<<". "<<"\n"<<std::endl;
              
              KRATOS_ERROR(std::runtime_error,"Please use a Absolute tolerance instead "," ")
              
              //NOTE: if it doesn't converge, problems occur with contact mesh and rigid face contact.
            }
            
      }
      
    } //SetCoordinationNumber
    
    double ComputeCoordinationNumber()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
         
        unsigned int total_contacts = 0;
                
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
        
        #pragma omp parallel 
        {
        mNeighbourCounter[OpenMPUtils::ThisThread()] = 0;
        
        #pragma omp for //private(index, MassMatrix)  //M. proba de compilar sense mass matrix??
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){

                //mNeighbourCounter[OpenMPUtils::ThisThread()] += (it)->GetValue(NEIGHBOUR_ELEMENTS).size();
                SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*it);
                mNeighbourCounter[OpenMPUtils::ThisThread()] += spheric_particle.mNeighbourElements.size();
                
            }

        }
        } //#pragma omp parallel 
              
        for (int i = 0; i < this->GetNumberOfThreads(); i++)
        {
          //std::cout<<mNeighbourCounter[i]<<"*********";
          total_contacts += mNeighbourCounter[i];
                  
        }
        
        double coord_number = (double(total_contacts)/double(pElements.size()));
        //std::cout<<"COORDINATION NUMBER = "<<coord_number<<std::endl;
        return coord_number;
       
        KRATOS_CATCH("")
      
    }
    


    void CalculateEnergies(){}


    void InitializeElements()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for 
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                (it)->Initialize();
            }

        }

        mElementsAreInitialized = true;

        KRATOS_CATCH("")
    }
    
    void InitializeFEMElements()
    {
        KRATOS_TRY


        ConditionsArrayType& pTContitions      = mpFem_model_part->GetCommunicator().LocalMesh().Conditions(); 
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pTContitions.size(), this->GetElementPartition());

        #pragma omp parallel for 
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ConditionsArrayType::iterator it_begin = pTContitions.ptr_begin() + this->GetElementPartition()[k];
            typename ConditionsArrayType::iterator it_end   = pTContitions.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it){

                (it)->Initialize();

            }

        }

        KRATOS_CATCH("")
    }
    
    //DEMFFEM
    
    void Calculate_Conditions_RHS_and_Add()
    {
      
      KRATOS_TRY
      
      Clear_forces_FEM();
              
      ConditionsArrayType& pConditions      = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();     

      ProcessInfo& CurrentProcessInfo  = GetFemModelPart().GetProcessInfo();

      Vector rhs_cond;

      vector<unsigned int> condition_partition;
      OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pConditions.size(), condition_partition);
      unsigned int index;
      
      #pragma omp parallel for private (index, rhs_cond)     
      
      for(int k=0; k<this->GetNumberOfThreads(); k++)
      {
          typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
          typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];

          for (typename ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
              Condition::GeometryType& geom = it->GetGeometry();
              
              it->CalculateRightHandSide(rhs_cond,CurrentProcessInfo);
            
              const unsigned int& dim = geom.WorkingSpaceDimension();
              for (unsigned int i = 0; i <geom.size(); i++)
              {
                  index = i*dim;//*2;
                  array_1d<double,3>& node_rhs = geom(i)->FastGetSolutionStepValue(ELASTIC_FORCES);//TOTAL_FORCES
                 // array_1d<double,3>& node_elastic_rhs = geom(i)->FastGetSolutionStepValue(ELASTIC_FORCES);
                  
                  for(unsigned int kk=0; kk<dim; kk++)
                  {
                      geom(i)->SetLock();

                      node_rhs[kk] = node_rhs[kk] + rhs_cond[index+kk];
                      //node_elastic_rhs[kk] = node_elastic_rhs[kk] + rhs_cond[index+kk+3];
                      geom(i)->UnSetLock();
                  }
                  
              }                   
              
          }
          
      }

      KRATOS_CATCH("")
    }
    
    
    void Clear_forces_FEM()

    {
        KRATOS_TRY

        ModelPart& fem_model_part  = GetFemModelPart();
        NodesArrayType& pNodes   = fem_model_part.Nodes();

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pNodes.size(), node_partition);

        #pragma omp parallel for
        
        for(int k=0; k<this->GetNumberOfThreads(); k++)
        {
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
                array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(ELASTIC_FORCES));
                //array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(TOTAL_FORCES));
                noalias(node_rhs)             = ZeroVector(3);
            }
        }

        KRATOS_CATCH("")
    }
    
    
     void ApplyPrescribedBoundaryConditions()
    {
      
      KRATOS_TRY

      ModelPart& r_model_part           = BaseType::GetModelPart();

      for (ModelPart::MeshesContainerType::iterator mesh_it = r_model_part.GetMeshes().begin(); mesh_it != r_model_part.GetMeshes().end(); ++mesh_it)
      {

          bool fix_x = bool((*mesh_it)[IMPOSED_VELOCITY_X]);
          bool fix_y = bool((*mesh_it)[IMPOSED_VELOCITY_Y]);
          bool fix_z = bool((*mesh_it)[IMPOSED_VELOCITY_Z]);
          
         if( fix_x || fix_y || fix_z )
         {
         
          double vel_x = (*mesh_it)[IMPOSED_VELOCITY_X_VALUE];
          double vel_y = (*mesh_it)[IMPOSED_VELOCITY_Y_VALUE];  
          double vel_z = (*mesh_it)[IMPOSED_VELOCITY_Z_VALUE];  

          NodesArrayType& pNodes = mesh_it->Nodes();
        
          vector<unsigned int> node_partition;
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pNodes.size(), node_partition);

          #pragma omp parallel for
          
          for(int k=0; k<this->GetNumberOfThreads(); k++)
          {
              typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
              typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

              for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
              {
                array_1d<double, 3>& velocity = i->FastGetSolutionStepValue(VELOCITY);
                
                if(fix_x)
                {
                  velocity[0] = vel_x;
                  i->Set(DEMFlags::FIXED_VEL_X,true);
                }
                
                if(fix_y)
                {
                  velocity[1] = vel_y;
                  i->Set(DEMFlags::FIXED_VEL_Y,true);
                }
                
                if(fix_z)
                {
                  velocity[2] = vel_z;
                  i->Set(DEMFlags::FIXED_VEL_Z,true);                
                }
     
              } //loop over particles

            }// loop threads OpenMP
            
          } //if(fix_x || fix_y || fix_z)
        
      } //for each mesh
      
      KRATOS_CATCH("")
    }
     
    
    void SetSearchRadius(ModelPart& r_model_part, double amplification)
    {
        KRATOS_TRY
     
        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        
        mNumberOfElementsOldRadiusList = number_of_elements;
        
        this->GetRadius().resize(number_of_elements);

        for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it){

            this->GetRadius()[particle_pointer_it - pElements.begin()] = amplification*(mSearchTolerance + particle_pointer_it->GetGeometry()(0)->GetSolutionStepValue(RADIUS));

        }

        KRATOS_CATCH("")
    }
    
    void SetOriginalRadius(ModelPart& r_model_part)
    {
        KRATOS_TRY
     
        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        
        this->GetOriginalRadius().resize(number_of_elements);

        for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it){

            this->GetOriginalRadius()[particle_pointer_it - pElements.begin()] = particle_pointer_it->GetGeometry()(0)->GetSolutionStepValue(RADIUS);

        }

        KRATOS_CATCH("")
    }

    void KRATOS_CHECK_SIZE(const char* msg, int id)
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it){
            
            SphericParticle* spheric_central_particle = dynamic_cast<Kratos::SphericParticle*>( &(*particle_pointer_it) );
            if (particle_pointer_it->Id() == id){
              //  std::cout << msg << " " << particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS).size() << std::endl; break;                
                std::cout << msg << " " << spheric_central_particle->mNeighbourElements.size() << std::endl; break;                
            }
        }

        KRATOS_CATCH("")
    }

    virtual void SearchNeighbours()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

        if(!number_of_elements) return;
        
        for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); particle_pointer_it++){
                
            SphericParticle* spheric_central_particle = dynamic_cast<Kratos::SphericParticle*>( &(*particle_pointer_it) );
            spheric_central_particle->mNeighbourElements.clear();            
        }        
        

        mpSpSearch->SearchElementsInRadiusExclusive(r_model_part, this->GetRadius(), this->GetResults(), this->GetResultsDistances());
//         mpSpSearch->SearchElementsInRadiusExclusive(r_model_part,this->GetRadius(),this->GetResults());
        
//         for (SpatialSearch::ElementsContainerType::iterator i = r_model_part.GetCommunicator().GhostMesh().Elements().begin(); i != r_model_part.GetCommunicator().GhostMesh().Elements().end(); i++){
//             r_model_part.Elements().push_back(*i);
//         }
        
//         SynchronizeSolidMesh(r_model_part);
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
        
        #pragma omp parallel for
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            size_t ResultCounter = this->GetElementPartition()[k];                       

            for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++ResultCounter){
                SphericParticle* spheric_central_particle = dynamic_cast<Kratos::SphericParticle*>( &(*particle_pointer_it) );
                
                for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[ResultCounter].begin(); neighbour_it != this->GetResults()[ResultCounter].end(); ++neighbour_it){
                    Element* p_neighbour_element = (*neighbour_it).get();
                    SphericParticle* p_spheric_neighbour_particle = dynamic_cast<SphericParticle*>( p_neighbour_element );
                    spheric_central_particle->mNeighbourElements.push_back( p_spheric_neighbour_particle );
                }

                this->GetResults()[ResultCounter].clear();
                this->GetResultsDistances()[ResultCounter].clear();
                
            }

        }

        KRATOS_CATCH("")
    }           
    
     void ComputeNewNeighboursHistoricalData()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for 
        
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
            std::vector<int> mTempNeighboursIds;
            std::vector<array_1d<double, 3> > mTempNeighbourElasticContactForces;
            std::vector<array_1d<double, 3> > mTempNeighbourTotalContactForces;  
            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                Kratos::SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*it);           
                spheric_particle.ComputeNewNeighboursHistoricalData(mTempNeighboursIds,mTempNeighbourElasticContactForces,mTempNeighbourTotalContactForces);
                          
            }
        }

        KRATOS_CATCH("")
    }
    
     void ComputeNewRigidFaceNeighboursHistoricalData()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements(); //ERROR? Should these be all the elements??
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for 
        
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
            
            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*it);
                spheric_particle.ComputeNewRigidFaceNeighboursHistoricalData();
            }

        }

        KRATOS_CATCH("")
    }

    
    virtual void SearchRigidFaceNeighbours()
    {
        
        KRATOS_TRY
        ElementsArrayType& pElements           = mpDem_model_part->GetCommunicator().LocalMesh().Elements();    
        ConditionsArrayType& pTContitions      = mpFem_model_part->GetCommunicator().LocalMesh().Conditions(); 
        
        if(pTContitions.size() > 0)
        {
            OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
            
            #pragma omp parallel for
            for (int k = 0; k < this->GetNumberOfThreads(); k++)
            {
                typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
                typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
                
                for (SpatialSearch::ElementsContainerType::iterator E_pointer_it = it_begin; E_pointer_it != it_end; ++E_pointer_it)
                {
                    SphericParticle* spheric_particle = dynamic_cast<Kratos::SphericParticle*>( &(*E_pointer_it) );
                    
                    spheric_particle->mNeighbourRigidFaces.resize(0);
                    
                    spheric_particle->mNeighbourRigidFacesPram.resize(0);
                    
                    spheric_particle->mNeighbourRigidFacesTotalContactForce.resize(0);
                    
                    spheric_particle->mNeighbourRigidFacesElasticContactForce.resize(0);
                }
            }
           
            moDemFemSearch.SearchRigidFaceForDEMInRadiusExclusiveImplementation(pElements, pTContitions, this->GetOriginalRadius(), this->GetRigidFaceResults(), this->GetRigidFaceResultsDistances());                        
            
            #pragma omp parallel
            {
                #pragma omp for
                for (int k = 0; k < this->GetNumberOfThreads(); k++)
                {
                    typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
                    typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

                    size_t ResultCounter = this->GetElementPartition()[k];
                    
                    for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++ResultCounter)
                    {
                        SphericParticle* spheric_particle = dynamic_cast<Kratos::SphericParticle*>( &(*particle_pointer_it) );
                        
                        std::vector<DEMWall*>& neighbour_rigid_faces = spheric_particle->mNeighbourRigidFaces;
                                        
                        for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[ResultCounter].begin(); 
                            neighbour_it != this->GetRigidFaceResults()[ResultCounter].end(); ++neighbour_it)
                        {                            
                            Condition* p_neighbour_condition = (*neighbour_it).get();
                            DEMWall* p_wall = dynamic_cast<DEMWall*>( p_neighbour_condition );
                            neighbour_rigid_faces.push_back(p_wall); 

                        }

                        this->GetRigidFaceResults()[ResultCounter].clear();
                        this->GetRigidFaceResultsDistances()[ResultCounter].clear();
                    
                    }

                }

                
                //Cfeng:resize the particle-rigidface contact forces for each particles 
                
                #pragma omp for
                for (int k = 0; k < this->GetNumberOfThreads(); k++)
                {
                    typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
                    typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
                    
                    for (SpatialSearch::ElementsContainerType::iterator E_pointer_it = it_begin; E_pointer_it != it_end; ++E_pointer_it)
                    {            
                        SphericParticle* spheric_particle = dynamic_cast<Kratos::SphericParticle*>( &(*E_pointer_it) );
                        std::size_t totalno = spheric_particle->mNeighbourRigidFaces.size() * 3;
                        spheric_particle->mNeighbourRigidFacesTotalContactForce.resize(totalno);
                        spheric_particle->mNeighbourRigidFacesElasticContactForce.resize(totalno);
                        
                    }
                }
            } //end of the parallel region
            
            typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
            

          int ConditionSize = pTContitions.size();
          int ElementSize = pElements.size();
          
            #pragma omp parallel 
            {
              
              #pragma omp for
              for (int i = 0; i<ConditionSize; i++)
              {
                
                  ConditionsArrayType::iterator ic = pTContitions.begin()+i;
                  DEMWall* wall = dynamic_cast<Kratos::DEMWall*>( &(*ic) );
                  wall->mNeighbourSphericParticles.resize(0);
              }       
              
              ////Cfeng: Find The particle neighbours for each RigidFace, used for calculating FEM force
              //This loop can not be parallelized easily, two threads could be writing on the same condition!!
              
              #pragma omp for
              for (int i = 0; i<ElementSize; i++)
              {   
                
                  SpatialSearch::ElementsContainerType::iterator E_pointer_it = pElements.begin()+i;
                  SphericParticle* spheric_particle = dynamic_cast<Kratos::SphericParticle*>( &(*E_pointer_it) );
                                   
                  for(unsigned int i=0; i<spheric_particle->mNeighbourRigidFaces.size(); i++) 
                  {
                      DEMWall* p_wall = spheric_particle->mNeighbourRigidFaces[i];
                      #pragma omp critical
                      {
                      p_wall->mNeighbourSphericParticles.push_back(spheric_particle);                                                                                      
                      }
                  }
              }
              
            }//end parallel
                                 
        }
                    
        KRATOS_CATCH("")
    }
                   
    void CalculateInitialMaxIndentations()
    {

        KRATOS_TRY

        double tol                          = 10e-18 * mpParticleCreatorDestructor->GetStrictDiameter();
        double initial_max_indentation      = 1.0 + tol;
        ModelPart& r_model_part             = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
        rCurrentProcessInfo.SetValue(DISTANCE_TOLERANCE, tol);
        ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();
        Vector reduction_distances;
        reduction_distances.resize(r_model_part.Elements().size());

        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        int elem_counter;

        #pragma omp parallel for private(elem_counter)
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
            elem_counter = 0;

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){

                if (it->GetGeometry()(0)->IsNot(DEMFlags::FIXED_VEL_X) ) {
                    reduction_distances[this->GetElementPartition()[k] + elem_counter] = initial_max_indentation;
                }

                else {
                    reduction_distances[this->GetElementPartition()[k] + elem_counter] = 0.0;
                }

                elem_counter++;
            }

        }

        #pragma omp parallel for private(elem_counter)
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
            elem_counter = 0;

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                double max_indentation = reduction_distances[this->GetElementPartition()[k] + elem_counter];
                it->Calculate(MAX_INDENTATION, max_indentation, rCurrentProcessInfo);
                reduction_distances[this->GetElementPartition()[k] + elem_counter] = 0.5 * max_indentation;
                elem_counter++;
            }

        }

        #pragma omp parallel for private(elem_counter)
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
            elem_counter = 0;

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                double reduction = reduction_distances[this->GetElementPartition()[k] + elem_counter];

                if (reduction > tol){
                    (it)->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS) -= reduction;
                }

                elem_counter++;
            }

        }

        KRATOS_CATCH("")
    
    } // CalculateInitialMaxIndentations()

    void PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part)
    {
        mcontacts_model_part.GetCommunicator().SetNumberOfColors(r_model_part.GetCommunicator().GetNumberOfColors());
        mcontacts_model_part.GetCommunicator().NeighbourIndices() = r_model_part.GetCommunicator().NeighbourIndices();
    }

    void SynchronizeSolidMesh(ModelPart& r_model_part)
    {
        r_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
//         r_model_part.GetCommunicator().SynchronizeDofs();
    }
    
    
    void DoAllOperations(                
                            ModelPart& DEM_inlet_model_part,
                            ParticleCreatorDestructor& creator_destructor,
                            DEMFEMUtilities& mesh_motion,
                            DEM_Inlet& DEM_inlet,
                            const std::string& dem_inlet_element_type,
                            double final_time,
                            double OutputTimeStep,
                            int total_steps_expected,
                            double ControlTime,
                            const std::string& main_path) {
        
            KRATOS_TRY;
            double time = 0.0;
            int step = 0;
            double time_old_print = 0.0;
            double time_to_print = 0.0;
            double& dt = mpDem_model_part->GetProcessInfo()[DELTA_TIME];
            double incremental_time = 0.0;
            double initial_real_time = std::time(0);
            double prev_time=0.0;
            //bool first_print = true;

            while (time < final_time) {

                time = time + dt;
                mpDem_model_part->GetProcessInfo()[TIME] = time;
                mpDem_model_part->GetProcessInfo()[TIME_STEPS] = step;

                //walls movement:
                mesh_motion.MoveAllMeshes(*mpFem_model_part, time);

                //# _SOLVE_###########################################
                //os.chdir(main_path);
                this->Solve();
                //# TIME CONTROL######################################

                //# adding DEM elements by the inlet:
                //if (inlet_option):
                DEM_inlet.CreateElementsFromInletMesh(*mpDem_model_part, DEM_inlet_model_part, creator_destructor, dem_inlet_element_type); //#After solving, to make sure that neighbours are already set.        

                incremental_time = (std::time(0) - initial_real_time) - prev_time;

                if (incremental_time > ControlTime){
                    double percentage = 100.0 * (float(step) / total_steps_expected);
        
                    std::cout<<"Real time calculation: "<<std::time(0)-initial_real_time<<" s"<<std::endl;
                    std::cout<<"Simulation time: "<<time<<" s"<<std::endl;
                    std::cout<<"Percentage Completed: "<<percentage<<" %"<<std::endl;   
                    std::cout<<"Time Step: "<<step<<std::endl<<std::endl<<std::flush;

                    prev_time = (std::time(0) - initial_real_time);
                }    

                time_to_print = time - time_old_print;

                if (time_to_print >= OutputTimeStep) {

                    /*print("")
                    print("*******************  PRINTING RESULTS FOR GID  ***************************")
                    print("                        (", balls_model_part.NumberOfElements(0), " elements)")
                    sys.stdout.flush()
        
                    # BENCHMARK ###
                    os.chdir(data_and_results)

                    # properties_list = proc.MonitorPhysicalProperties(balls_model_part, physics_calculator, properties_list)

                    if (index_5 == 5):
                        multifile_5.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
                        index_5 = 0

                    if (index_10 == 10):
                        multifile_10.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
                        index_10 = 0

                    if (index_50 == 50):
                        multifile_50.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
                        index_50 = 0

                    index_5 += 1
                    index_10 += 1
                    index_50 += 1

                    if (DEM_parameters.Multifile == "multiple_files"):
                        gid_io.FinalizeResults()

                    os.chdir(post_path)

                    if (DEM_parameters.Multifile == "multiple_files"):
                        mixed_model_part.Elements.clear()
                        mixed_model_part.Nodes.clear()

                        post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
                        post_utility.AddModelPartToModelPart(mixed_model_part, mpFem_model_part)

                        gid_io.InitializeMesh(time)
                        gid_io.WriteSphereMesh(balls_model_part.GetMesh())
                        gid_io.WriteMesh(mpFem_model_part.GetMesh())
                        gid_io.FinalizeMesh()

                        gid_io.InitializeResults(time, mixed_model_part.GetMesh())

                    proc.PrintingGlobalVariables(gid_io, mixed_model_part, time)
                    proc.PrintingBallsVariables(gid_io, balls_model_part, time)

                    if (DEM_parameters.Multifile == "multiple_files"):
                        gid_io.FinalizeResults()

                    time_old_print = time;*/
                }

                step += 1;
            }//while

            KRATOS_CATCH("");
        }
    
    
    
    
    
    
    
    
    // Getting member variables
    
    ModelPart&                                   GetBallsModelPart(){return(*mpDem_model_part);}
    ModelPart&                                   GetFemModelPart(){return(*mpFem_model_part);}
    
    VectorResultElementsContainerType&           GetResults(){return(mResults);}
    VectorDistanceType&                          GetResultsDistances(){return(mResultsDistances);}
    RadiusArrayType&                             GetRadius(){return(mRadius);}
    RadiusArrayType&                             GetOriginalRadius(){return(mOriginalRadius);}

    bool&                                        GetElementsAreInitialized(){return (mElementsAreInitialized);}
    bool&                                        GetInitializeWasPerformed(){return (mInitializeWasPerformed);}
    
    int&                                         GetNStepSearch(){return (mNStepSearch);}
    int&                                         GetBoundingBoxOption(){return (mBoundingBoxOption);}
    int&                                         GetNumberOfThreads(){return (mNumberOfThreads);}

    double&                                      GetMaxTimeStep(){return (mMaxTimeStep);}
    double&                                      GetSafetyFactor(){return (mSafetyFactor);}
    
    int&                                         GetDeltaOption(){return (mDeltaOption);}
    double&                                      GetSearchTolerance(){return (mSearchTolerance);}
    double&                                      GetCoordinationNumber(){return (mCoordinationNumber);}
    vector<unsigned int>&                        GetNeighbourCounter(){return(mNeighbourCounter);}
                                           
    int&                                         GetNumberOfElementsOldRadiusList(){return (mNumberOfElementsOldRadiusList);}

    vector<unsigned int>&                        GetElementPartition(){return (mElementPartition);}

    typename ParticleCreatorDestructor::Pointer& GetParticleCreatorDestructor(){return (mpParticleCreatorDestructor);}
    typename IntegrationScheme::Pointer&         GetScheme(){return (mpScheme);}
    typename SpatialSearch::Pointer&             GetSpSearch(){return (mpSpSearch);}
    
    //Cfeng
    VectorResultConditionsContainerType&         GetRigidFaceResults(){return(mRigidFaceResults);}
    VectorDistanceType&                          GetRigidFaceResultsDistances(){return(mRigidFaceResultsDistances);}
    vector<unsigned int>&                        GetConditionPartition(){return (mConditionPartition);}
    
    std::vector<PropertiesProxy>                 mFastProperties;
            
    
    protected:

    // Member variables

    VectorResultElementsContainerType            mResults;
    VectorDistanceType                           mResultsDistances;
    RadiusArrayType                              mRadius;
    RadiusArrayType                              mOriginalRadius;

    bool                                         mElementsAreInitialized;
    bool                                         mInitializeWasPerformed;

    int                                          mNStepSearch;
    int                                          mBoundingBoxOption;
    int                                          mNumberOfThreads;
    int                                          mNumberOfElementsOldRadiusList;

    double                                       mMaxTimeStep;
    double                                       mSafetyFactor;
    
    int                                          mDeltaOption;
    double                                       mSearchTolerance;
    double                                       mCoordinationNumber;
    vector<unsigned int>                         mNeighbourCounter;
          
    vector<unsigned int>                         mElementPartition;
    typename ParticleCreatorDestructor::Pointer  mpParticleCreatorDestructor;
    typename IntegrationScheme::Pointer          mpScheme;
    typename SpatialSearch::Pointer              mpSpSearch;
    
    
    ////Cfeng
    VectorResultConditionsContainerType  mRigidFaceResults;
    VectorDistanceType                   mRigidFaceResultsDistances;
    DEM_FEM_Search                       moDemFemSearch;
    vector<unsigned int>                 mConditionPartition;
    ModelPart                            *mpFem_model_part;
    ModelPart                            *mpDem_model_part;
    
  }; // Class ExplicitSolverStrategy

 /*
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                    ExplicitSolverStrategy& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                    const ExplicitSolverStrategy& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
    */
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined

