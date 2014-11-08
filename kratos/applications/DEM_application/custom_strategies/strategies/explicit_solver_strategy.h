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
#include "custom_utilities/GeometryFunctions.h"

#include "custom_elements/cluster3D.h"

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
    
  class ExplicitSolverSettings {
  public:
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverSettings);
    ExplicitSolverSettings(){}
    ~ExplicitSolverSettings(){}
        ModelPart::Pointer r_model_part;
        ModelPart::Pointer contact_model_part;
        ModelPart::Pointer fem_model_part;
        ModelPart::Pointer cluster_model_part;
        ModelPart::Pointer inlet_model_part;
      };

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
      
      ExplicitSolverStrategy(ExplicitSolverSettings& settings,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool move_mesh_flag,
                             const int delta_option,
                             const double search_tolerance,
                             const double coordination_number,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename IntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch)
      :
      BaseType( *(settings.r_model_part),move_mesh_flag)
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
          
          
          mpDem_model_part = &(*(settings.r_model_part));
          if ( mpDem_model_part == NULL )
            KRATOS_ERROR(std::runtime_error,"Undefined settings.r_model_part in ExplicitSolverStrategy constructor","")
            
          mpContact_model_part            = &(*(settings.contact_model_part));
           if ( mpContact_model_part == NULL )
            KRATOS_ERROR(std::runtime_error,"Undefined settings.contact_model_part in ExplicitSolverStrategy constructor","")
          
          mpFem_model_part                = &(*(settings.fem_model_part));
           if ( mpFem_model_part == NULL )
            KRATOS_ERROR(std::runtime_error,"Undefined settings.fem_model_part in ExplicitSolverStrategy constructor","")
         
          mpCluster_model_part              = &(*(settings.cluster_model_part));
           if ( mpCluster_model_part == NULL )
            KRATOS_ERROR(std::runtime_error,"Undefined settings.cluster_model_part in ExplicitSolverStrategy constructor","")
            
          mpInlet_model_part              = &(*(settings.inlet_model_part));
           if ( mpInlet_model_part == NULL )
            KRATOS_ERROR(std::runtime_error,"Undefined settings.inlet_model_part in ExplicitSolverStrategy constructor","")
                   
      }

      ExplicitSolverStrategy(ModelPart& r_model_part,
                             ModelPart& fem_model_part,
                             ModelPart& cluster_model_part,
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
          mpCluster_model_part           = &cluster_model_part;         
      }

      /// Destructor.
      virtual ~ExplicitSolverStrategy()
      {
          Timer::SetOuputFile("TimesPartialRelease");
          Timer::PrintTimingInformation();
      }            
      
      template <class T>
      void RebuildListOfSphericParticles(ElementsArrayType& pElements, std::vector<T*>& rCustomListOfParticles){
          
          KRATOS_TRY
              
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
          
          rCustomListOfParticles.resize( pElements.size() );
          #pragma omp parallel for
          for (int k = 0; k < this->GetNumberOfThreads(); k++){
              
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
            int elem_number = this->GetElementPartition()[k];

            for (typename ElementsArrayType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++elem_number){
                T* spheric_particle = dynamic_cast<T*>( &(*particle_pointer_it) );
                rCustomListOfParticles[elem_number] = spheric_particle;
            }                                
          }              
          
          return;          
          KRATOS_CATCH("")
          
      }                        
      
      void RebuildPropertiesProxyPointers(std::vector<SphericParticle*>& rCustomListOfSphericParticles){
          //This function is called for the local mesh and the ghost mesh, so mListOfSphericElements must not be used here.
          KRATOS_TRY         
                                  
          #pragma omp parallel for
          for(int i=0; i<(int)rCustomListOfSphericParticles.size(); i++){  
              rCustomListOfSphericParticles[i]->SetFastProperties(mFastProperties);                         
          }                     
          return;            
          KRATOS_CATCH("")
      }
      
      void SendProcessInfoToClustersModelPart(){
          
          ProcessInfo& rCurrentProcessInfo         = mpDem_model_part->GetProcessInfo();
          ProcessInfo& rClustersCurrentProcessInfo = mpCluster_model_part->GetProcessInfo();
          
          rCurrentProcessInfo[CONTAINS_CLUSTERS] = false;
          rClustersCurrentProcessInfo[CONTAINS_CLUSTERS] = true;
                  
          rClustersCurrentProcessInfo[ROTATION_OPTION]     = rCurrentProcessInfo[ROTATION_OPTION];    
          rClustersCurrentProcessInfo[DELTA_TIME]          = rCurrentProcessInfo[DELTA_TIME];
          rClustersCurrentProcessInfo[VIRTUAL_MASS_OPTION] = rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
          rClustersCurrentProcessInfo[TRIHEDRON_OPTION]    = rCurrentProcessInfo[TRIHEDRON_OPTION];
          rClustersCurrentProcessInfo[NODAL_MASS_COEFF]    = rCurrentProcessInfo[NODAL_MASS_COEFF];                        
      }
      
      void UpdateMaxIdOfCreatorDestructor(){
          ModelPart& r_model_part            = BaseType::GetModelPart();
          int max_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart( r_model_part );
          int max_FEM_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart( *mpFem_model_part );
          
          max_Id = std::max(max_Id, max_FEM_Id);                    
          mpParticleCreatorDestructor->SetMaxNodeId(max_Id);
          
      }
      
      void RepairPointersToNormalProperties(std::vector<SphericParticle*>& rCustomListOfSphericParticles){                   
          
          bool found = false;    
          //#pragma omp parallel for
          for (int i=0; i<(int)rCustomListOfSphericParticles.size(); i++){
                            
              int own_properties_id = rCustomListOfSphericParticles[i]->GetProperties().Id();  
              
              for (PropertiesIterator props_it = mpDem_model_part->GetMesh(0).PropertiesBegin(); props_it!= mpDem_model_part->GetMesh(0).PropertiesEnd();   props_it++ ) {
                  int model_part_id = props_it->GetId();
                  if (own_properties_id == model_part_id) {
                      rCustomListOfSphericParticles[i]->SetProperties( *(props_it.base()) );
                      found = true;
                      break;
                  }                  
              }
              
              if(found) continue;
              
              for (PropertiesIterator props_it = mpInlet_model_part->GetMesh(0).PropertiesBegin(); props_it!= mpInlet_model_part->GetMesh(0).PropertiesEnd();   props_it++ ) {
                  int model_part_id = props_it->GetId();
                  if (own_properties_id == model_part_id) {
                      rCustomListOfSphericParticles[i]->SetProperties( *(props_it.base()) );
                      found = true;
                      break;
                  }                  
              }
              
              if(found) continue;
              
              for (PropertiesIterator props_it = mpCluster_model_part->GetMesh(0).PropertiesBegin(); props_it!= mpCluster_model_part->GetMesh(0).PropertiesEnd();   props_it++ ) {
                  int model_part_id = props_it->GetId();
                  if (own_properties_id == model_part_id) {
                      rCustomListOfSphericParticles[i]->SetProperties( *(props_it.base()) );
                      found = true;
                      break;
                  }                  
              }
              
              if(!found) KRATOS_ERROR(std::logic_error, "This particle could not find its properties!!" , "");
          }                                                                   
      }
      
      virtual void Initialize()
      {
          KRATOS_TRY

          ModelPart& r_model_part            = BaseType::GetModelPart();
          
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();
          SendProcessInfoToClustersModelPart();          
                    
          // Omp initializations
          GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
          mNeighbourCounter.resize(this->GetNumberOfThreads());
          
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
          
          CreatePropertiesProxies(mFastProperties, *mpDem_model_part, *mpInlet_model_part, *mpCluster_model_part);
          
          RepairPointersToNormalProperties(mListOfSphericParticles);  //The particles sent to this partition have their own copy of the Kratos properties they were using in the previous partition!!
          RepairPointersToNormalProperties(mListOfGhostSphericParticles);
          
          RebuildPropertiesProxyPointers(mListOfSphericParticles);
          RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);        
                   
          GetBoundingBoxOption() = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

          InitializeSolutionStep();
          InitializeElements();                    
          InitializeFEMElements();
          UpdateMaxIdOfCreatorDestructor();
          InitializeClusters(); /// This adds elements to the balls modelpart

          mInitializeWasPerformed = true;
          
          ApplyPrescribedBoundaryConditions();
          
          // Search Neighbours and related operations                             
          SetOriginalRadius(r_model_part);
          SetSearchRadius(r_model_part, 1.0);          
          SearchNeighbours();
          
          if(mDeltaOption == 2)
          {
            
            SetCoordinationNumber(r_model_part);
            
          }
          
          ComputeNewNeighboursHistoricalData();  
          
          SearchRigidFaceNeighbours();
          ComputeNewRigidFaceNeighboursHistoricalData();
          
          // Finding overlapping of initial configurations

          if (rCurrentProcessInfo[CLEAN_INDENT_OPTION]){
              CalculateInitialMaxIndentations();
          }

          
           // 5. Finalize Solution Step.
          //FinalizeSolutionStep();
          //KRATOS_WATCH(r_model_part.GetNodalSolutionStepVariablesList())
                    

      KRATOS_CATCH("")
      }// Initialize()
    
      
      virtual void InitializeClusters() {
          
          KRATOS_TRY

          ElementsArrayType& pElements = mpCluster_model_part->GetCommunicator().LocalMesh().Elements();    
          
          typename ElementsArrayType::iterator it_begin = pElements.ptr_begin();
          typename ElementsArrayType::iterator it_end = pElements.ptr_end();
          
          for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
      
              Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&>(*it);
              
              (cluster_element).Initialize();
              (cluster_element).CreateParticles(mpParticleCreatorDestructor, *mpDem_model_part);
              
          }
          
          KRATOS_CATCH("")
        
      }

        virtual void GetClustersForce() {

            KRATOS_TRY
                    
            ProcessInfo& rCurrentProcessInfo    = BaseType::GetModelPart().GetProcessInfo(); //Getting the Process Info of the Balls ModelPart!
            const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY]; 

            ElementsArrayType& pElements = mpCluster_model_part->GetCommunicator().LocalMesh().Elements();

            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin();
            typename ElementsArrayType::iterator it_end = pElements.ptr_end();

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {

                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                
                cluster_element.GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES).clear();
                cluster_element.GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT).clear();
                
                (cluster_element).GetClustersForce( gravity );
                
            }

            KRATOS_CATCH("")
        }
        
      
      virtual double Solve()
      {
          KRATOS_TRY

          ModelPart& r_model_part            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();      
          
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

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
                  RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
                  RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
              }
              
              SetSearchRadius(r_model_part, 1.0);        
              SearchNeighbours();
              ComputeNewNeighboursHistoricalData();  
              
              SetOriginalRadius(r_model_part);              
              SearchRigidFaceNeighbours();
              ComputeNewRigidFaceNeighboursHistoricalData();
              
          }
           
          RebuildPropertiesProxyPointers(mListOfSphericParticles);
          RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);        
          
          // 3. Get and Calculate the forces
          GetForce();
          //FastGetForce();
          
          GetClustersForce();
          
          Calculate_Conditions_RHS_and_Add();
          
          Calculate_Nodal_Pressures_and_Stresses();
          
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

          ProcessInfo& rCurrentProcessInfo    = BaseType::GetModelPart().GetProcessInfo();
          double dt = rCurrentProcessInfo[DELTA_TIME];
          const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY];             

          #pragma omp parallel
          {
            Vector rhs_elem;
            rhs_elem.resize(6);
            #pragma omp for
            for(int i=0; i<(int)mListOfSphericParticles.size(); i++){              
                mListOfSphericParticles[i]->CalculateRightHandSide(rhs_elem, rCurrentProcessInfo, dt, gravity);  
            }
          }

          KRATOS_CATCH("")
      }
      
      void ClusterGetForce() {
          
          KRATOS_TRY

          ElementsArrayType& pElements = mpCluster_model_part->GetCommunicator().LocalMesh().Elements();    
          
          typename ElementsArrayType::iterator it_begin = pElements.ptr_begin();
          typename ElementsArrayType::iterator it_end = pElements.ptr_end();
                            
          for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
      
              //Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&>(*it);
                  
          }
        
          KRATOS_CATCH("")
        
      }
            
      void FastGetForce()
      {
          KRATOS_TRY
          
          ProcessInfo& rCurrentProcessInfo    = BaseType::GetModelPart().GetProcessInfo();
          double dt = rCurrentProcessInfo[DELTA_TIME];
          const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY];

          #pragma omp parallel
          {
              #pragma omp for
              for(int i=0; i<(int)mListOfSphericParticles.size(); i++){ 
                  mListOfSphericParticles[i]->FirstCalculateRightHandSide(rCurrentProcessInfo, dt); 
              }
              #pragma omp for
              for(int i=0; i<(int)mListOfSphericParticles.size(); i++){ 
                  mListOfSphericParticles[i]->CollectCalculateRightHandSide(rCurrentProcessInfo); 
              }
              #pragma omp for
              for(int i=0; i<(int)mListOfSphericParticles.size(); i++){ 
                  mListOfSphericParticles[i]->FinalCalculateRightHandSide(rCurrentProcessInfo, dt, gravity); 
              }

          }

          KRATOS_CATCH("")
      }

      virtual void PerformTimeIntegrationOfMotion(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY
          
          ModelPart& r_model_part = BaseType::GetModelPart();

          GetScheme()->Calculate(r_model_part);
          GetScheme()->Calculate(*mpCluster_model_part);

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
          mpParticleCreatorDestructor->DestroyParticlesOutsideBoundingBox(*mpCluster_model_part);        
          mpParticleCreatorDestructor->DestroyParticlesOutsideBoundingBox(r_model_part);
          

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
           
        #pragma omp parallel 
        {
            mNeighbourCounter[OpenMPUtils::ThisThread()] = 0;        
            #pragma omp for
            for(int i=0; i<(int)mListOfSphericParticles.size(); i++){              
                  mNeighbourCounter[OpenMPUtils::ThisThread()] += mListOfSphericParticles[i]->mNeighbourElements.size();
            }
        }                                      
        for (int i = 0; i < this->GetNumberOfThreads(); i++)
        {
          total_contacts += mNeighbourCounter[i];                  
        }
        
        int global_total_contacts = total_contacts; 
        r_model_part.GetCommunicator().SumAll(global_total_contacts);
        int global_number_of_elements = (int)pElements.size();
        r_model_part.GetCommunicator().SumAll(global_number_of_elements);
                
        double coord_number = double(global_total_contacts)/double(global_number_of_elements);              
        
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
                CalculateNormals(it);
            }
        }

        KRATOS_CATCH("")
    }
    
    //DEMFEM
    
    void CalculateNormals(ConditionsArrayType::iterator it) {

        KRATOS_TRY
        
        array_1d<double, 3>& NormalToCondition = it->GetValue(NORMAL);
         
        //Calculate the normal vector to an element
        
        array_1d<double, 3> v1, v2;
        
        v1[0] = it->GetGeometry()[1].X() - it->GetGeometry()[0].X();
        v1[1] = it->GetGeometry()[1].Y() - it->GetGeometry()[0].Y();
        v1[2] = it->GetGeometry()[1].Z() - it->GetGeometry()[0].Z();

        v2[0] = it->GetGeometry()[2].X() - it->GetGeometry()[0].X();
        v2[1] = it->GetGeometry()[2].Y() - it->GetGeometry()[0].Y();
        v2[2] = it->GetGeometry()[2].Z() - it->GetGeometry()[0].Z();

        MathUtils<double>::CrossProduct(NormalToCondition, v1, v2);
        
        NormalToCondition /= MathUtils<double>::Norm3(NormalToCondition);
    
        KRATOS_CATCH("")
    }
    
    
    void Calculate_Conditions_RHS_and_Add() {
      
        KRATOS_TRY
      
        Clear_forces_FEM();
        ConditionsArrayType& pConditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();     
        ProcessInfo& CurrentProcessInfo = GetFemModelPart().GetProcessInfo();
        
        Vector rhs_cond;
        vector<unsigned int> condition_partition;    
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pConditions.size(), condition_partition);    
        unsigned int index;
               
        #pragma omp parallel for private (index, rhs_cond)
        
        for (int k=0; k<this->GetNumberOfThreads(); k++) {
            
            typename ConditionsArrayType::iterator it_begin = pConditions.ptr_begin() + condition_partition[k];        
            typename ConditionsArrayType::iterator it_end = pConditions.ptr_begin() + condition_partition[k+1];
            
            for (typename ConditionsArrayType::iterator it = it_begin; it!= it_end; ++it) { //each iteration refers to a different triangle or quadrilateral
                
                Condition::GeometryType& geom = it->GetGeometry();            
                double Element_Area = geom.Area();                     
                it->CalculateRightHandSide(rhs_cond, CurrentProcessInfo);           
                array_1d<double, 3> Normal_to_Element = it->GetValue(NORMAL);                                                                       
                const unsigned int& dim = geom.WorkingSpaceDimension();
                
                for (unsigned int i = 0; i <geom.size(); i++) { //talking about each of the three nodes of the condition  
                                                                //we are studying a certain condition here
                    index = i * dim;    //*2;                 
                    geom(i)->SetLock();                    
                    
                    array_1d<double, 3>& node_rhs = geom(i)->GetSolutionStepValue(ELASTIC_FORCES); 
                    array_1d<double, 3>& node_rhs_tang = geom(i)->GetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
                    double& node_pressure = geom(i)->GetSolutionStepValue(PRESSURE);                  
                    double& node_area = geom(i)->GetSolutionStepValue(NODAL_AREA);
                    array_1d<double, 3> rhs_cond_comp;
                    
                    for (unsigned int j = 0; j < dim; j++) { //talking about each coordinate x, y and z, loop on them                   
                        
                        node_rhs[j] = node_rhs[j] + rhs_cond[index+j];
                        rhs_cond_comp[j] =  rhs_cond[index+j];
                    }
                                                           
                    node_area += 0.333333333333333 * Element_Area;
                    //node_pressure actually refers to normal force. It is really computed later in function Calculate_Nodal_Pressures_and_Stresses()
                    //as a real pressure
                    node_pressure += MathUtils<double>::Abs(GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element));
                    
                    node_rhs_tang += rhs_cond_comp - GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element) * Normal_to_Element;
                    
                    
                    geom(i)->UnSetLock();
                    
                }                                 
            }          
        }

        KRATOS_CATCH("")
    }
    
    
    void Clear_forces_FEM() {
        
        KRATOS_TRY

        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& pNodes = fem_model_part.Nodes();

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pNodes.size(), node_partition);

        #pragma omp parallel for
        
        for (int k=0; k<this->GetNumberOfThreads(); k++) {
            
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin() + node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin() + node_partition[k+1];

            for (ModelPart::NodeIterator i=i_begin; i!= i_end; ++i) {
                
                array_1d<double, 3>& node_rhs = i->FastGetSolutionStepValue(ELASTIC_FORCES);
                array_1d<double, 3>& node_rhs_tang = i->FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
                double& node_pressure = i->FastGetSolutionStepValue(PRESSURE);
                double& node_area = i->FastGetSolutionStepValue(NODAL_AREA);
                double& shear_stress = i->FastGetSolutionStepValue(SHEAR_STRESS);
                
                noalias(node_rhs) = ZeroVector(3);
                noalias(node_rhs_tang) = ZeroVector(3);
                node_pressure = 0.0;
                node_area = 0.0;
                shear_stress = 0.0;
                
            }
            
        }

        KRATOS_CATCH("")
    }
    
    
    void Calculate_Nodal_Pressures_and_Stresses() {
        
        KRATOS_TRY

        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& pNodes = fem_model_part.Nodes();

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pNodes.size(), node_partition);

        #pragma omp parallel for
        
        for (int k = 0; k < this->GetNumberOfThreads(); k++) {
            
            typename NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
            typename NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k+1];

            for (ModelPart::NodeIterator i = i_begin; i!= i_end; ++i) {
                
                double& node_pressure = i->FastGetSolutionStepValue(PRESSURE);
                double& node_area = i->FastGetSolutionStepValue(NODAL_AREA);
                double& shear_stress = i->FastGetSolutionStepValue(SHEAR_STRESS);
                array_1d<double, 3>& node_rhs_tang = i->FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);

                node_pressure = node_pressure/node_area;
                shear_stress = GeometryFunctions::module(node_rhs_tang)/node_area;
                 
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

    virtual void SearchNeighbours()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();        
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        if(!number_of_elements) return;    
        
        GetResults().resize(number_of_elements);
        GetResultsDistances().resize(number_of_elements);                              
        GetRigidFaceResults().resize(number_of_elements);
        GetRigidFaceResultsDistances().resize(number_of_elements);
        
        mpSpSearch->SearchElementsInRadiusExclusive(r_model_part, this->GetRadius(), this->GetResults(), this->GetResultsDistances());
        
        #pragma omp parallel for
        for (int i=0; i<(int)mListOfSphericParticles.size(); i++){
            mListOfSphericParticles[i]->mNeighbourElements.clear();
            for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[i].begin(); neighbour_it != this->GetResults()[i].end(); ++neighbour_it){
                Element* p_neighbour_element = (*neighbour_it).get();
                SphericParticle* p_spheric_neighbour_particle = dynamic_cast<SphericParticle*>( p_neighbour_element );
                if ( mListOfSphericParticles[i]->Is(DEMFlags::BELONGS_TO_A_CLUSTER) &&  ( mListOfSphericParticles[i]->GetClusterId() == p_spheric_neighbour_particle->GetClusterId() ) ) continue;
                mListOfSphericParticles[i]->mNeighbourElements.push_back( p_spheric_neighbour_particle );
            }
            this->GetResults()[i].clear();
            this->GetResultsDistances()[i].clear();        
        }
                
        KRATOS_CATCH("")
    }           
    
     void ComputeNewNeighboursHistoricalData()
    {
        KRATOS_TRY  

        #pragma omp parallel
        {
            std::vector<int> mTempNeighboursIds;
            std::vector<array_1d<double, 3> > mTempNeighbourElasticContactForces;
            std::vector<array_1d<double, 3> > mTempNeighbourTotalContactForces;
        
            #pragma omp for
            for(int i=0; i<(int)mListOfSphericParticles.size(); i++){
                mListOfSphericParticles[i]->ComputeNewNeighboursHistoricalData(mTempNeighboursIds,mTempNeighbourElasticContactForces,mTempNeighbourTotalContactForces);
            }
        }

        KRATOS_CATCH("")
    }
    
     void ComputeNewRigidFaceNeighboursHistoricalData()
    {
        KRATOS_TRY

        #pragma omp parallel for 
        for(int i=0; i<(int)mListOfSphericParticles.size(); i++){
            mListOfSphericParticles[i]->ComputeNewRigidFaceNeighboursHistoricalData();
        }

        KRATOS_CATCH("")
    }

    
    virtual void SearchRigidFaceNeighbours()
    {
        
        KRATOS_TRY
        ElementsArrayType& pElements           = mpDem_model_part->GetCommunicator().LocalMesh().Elements();    
        ConditionsArrayType& pTContitions      = mpFem_model_part->GetCommunicator().LocalMesh().Conditions(); 
        
        if(pTContitions.size() > 0) {
        
            OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
            #pragma omp parallel for 
            for(int i=0; i<(int)mListOfSphericParticles.size(); i++){
                mListOfSphericParticles[i]->mNeighbourRigidFaces.resize(0);
                mListOfSphericParticles[i]->mNeighbourRigidFacesPram.resize(0);
                mListOfSphericParticles[i]->mNeighbourRigidFacesTotalContactForce.resize(0);
                mListOfSphericParticles[i]->mNeighbourRigidFacesElasticContactForce.resize(0);
            }
           
            moDemFemSearch.SearchRigidFaceForDEMInRadiusExclusiveImplementation(pElements, pTContitions, this->GetOriginalRadius(), this->GetRigidFaceResults(), this->GetRigidFaceResultsDistances());                        
            
            #pragma omp parallel for
            for( int i=0; i<(int)mListOfSphericParticles.size(); i++ ){
                std::vector<DEMWall*>& neighbour_rigid_faces = mListOfSphericParticles[i]->mNeighbourRigidFaces;
                for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[i].begin(); neighbour_it != this->GetRigidFaceResults()[i].end(); ++neighbour_it){                                                  
                    Condition* p_neighbour_condition = (*neighbour_it).get();
                    DEMWall* p_wall = dynamic_cast<DEMWall*>( p_neighbour_condition );
                    neighbour_rigid_faces.push_back(p_wall); 
                }
                this->GetRigidFaceResults()[i].clear();
                this->GetRigidFaceResultsDistances()[i].clear();   

                std::size_t totalno = mListOfSphericParticles[i]->mNeighbourRigidFaces.size() * 3;
                mListOfSphericParticles[i]->mNeighbourRigidFacesTotalContactForce.resize(totalno);
                mListOfSphericParticles[i]->mNeighbourRigidFacesElasticContactForce.resize(totalno);
            }
                                     
            typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
            
            #pragma omp parallel 
            {              
                #pragma omp for
                for (int i = 0; i<(int)pTContitions.size(); i++)
                {               
                    ConditionsArrayType::iterator ic = pTContitions.begin()+i;
                    DEMWall* wall = dynamic_cast<Kratos::DEMWall*>( &(*ic) );
                    wall->mNeighbourSphericParticles.resize(0);
                }       

                #pragma omp for
                for( int i=0; i<(int)mListOfSphericParticles.size(); i++ ){                                  
                    for(unsigned int j=0; j<mListOfSphericParticles[i]->mNeighbourRigidFaces.size(); j++) 
                    {
                        DEMWall* p_wall = mListOfSphericParticles[i]->mNeighbourRigidFaces[j];
                        #pragma omp critical //TODO: What about doing this with locks?
                        {
                            p_wall->mNeighbourSphericParticles.push_back(mListOfSphericParticles[i]);                                                                                      
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
        std::vector<double> indentations_list;
        indentations_list.resize(mListOfSphericParticles.size());
                
        #pragma omp parallel for
        for( int i=0; i<(int)mListOfSphericParticles.size(); i++ ){
            double indentation;
            mListOfSphericParticles[i]->CalculateMaxBallToBallIndentation(indentation);
            double max_indentation = std::max(0.0, 0.5 * indentation); // reducing the radius by half the indentation is enough
            mListOfSphericParticles[i]->CalculateMaxBallToFaceIndentation(indentation);
            max_indentation = std::max(max_indentation, indentation);
            indentations_list[i] = max_indentation;
        }
        
        #pragma omp parallel for //THESE TWO LOOPS CANNOT BE JOINED, BECAUSE THE RADII ARE CHANGING.
        for( int i=0; i<(int)mListOfSphericParticles.size(); i++ ){
            mListOfSphericParticles[i]->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS) -= indentations_list[i];
            mListOfSphericParticles[i]->SetRadius(mListOfSphericParticles[i]->GetRadius() - indentations_list[i]); 
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
    ModelPart&                                   GetContactModelPart(){return(*mpContact_model_part);}
    ModelPart&                                   GetClusterModelPart(){return(*mpCluster_model_part);}
    ModelPart&                                   GetInletModelPart(){return(*mpInlet_model_part);}
    
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
    ModelPart                            *mpInlet_model_part;
    ModelPart                            *mpContact_model_part;
    ModelPart                            *mpCluster_model_part;
    
    std::vector<SphericParticle*>        mListOfSphericParticles;
    std::vector<SphericParticle*>        mListOfGhostSphericParticles;
    
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

