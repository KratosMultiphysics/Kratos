//
// Authors: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
// Miquel Santasusana msantasusana@cimne.upc.edu
//


#if !defined(KRATOS_EXPLICIT_SOLVER_STRATEGY)
#define KRATOS_EXPLICIT_SOLVER_STRATEGY

// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "DEM_application.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
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
#include "custom_strategies/schemes/dem_integration_scheme.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_utilities/dem_fem_utilities.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/inlet.h"

#include "custom_elements/cluster3D.h"

////Cfeng
#include "custom_utilities/dem_fem_search.h"
#include "custom_utilities/rigid_face_geometrical_object_configure.h"

#ifdef USING_CGAL
#include <CGAL/spatial_sort.h>
#endif

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

 
  template<
  class TSparseSpace,
  class TDenseSpace,
  class TLinearSolver>
  class ExplicitSolverStrategy: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
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
      typedef RigidFaceGeometricalObjectConfigure<3>                    RigidFaceGeometricalConfigureType;

      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);
      
              


      ///@}
      ///@name Life Cycle
      ///@{


      /// Default constructor.
      ExplicitSolverStrategy(){}
      
      ExplicitSolverStrategy(ExplicitSolverSettings& settings,
                             const double max_delta_time,
                             const int n_step_search,
                             const double safety_factor,
                             const int delta_option,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename DEM_FEM_Search::Pointer p_dem_fem_search,
                             typename DEMIntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch)
      :
      BaseType( *(settings.r_model_part), true)
      {

          mDeltaOption                   = delta_option;
          mpParticleCreatorDestructor    = p_creator_destructor;
          mpDemFemSearch                 = p_dem_fem_search;
          mpScheme                       = pScheme;
          mpSpSearch                     = pSpSearch;
          mMaxTimeStep                   = max_delta_time;
          mNStepSearch                   = n_step_search;
          mSafetyFactor                  = safety_factor;
          mNumberOfElementsOldRadiusList = 0;
          
          
          mpDem_model_part = &(*(settings.r_model_part));
          if ( mpDem_model_part == NULL )
            KRATOS_THROW_ERROR(std::runtime_error,"Undefined settings.r_model_part in ExplicitSolverStrategy constructor","")
            
          mpContact_model_part            = &(*(settings.contact_model_part));
           if ( mpContact_model_part == NULL )
            KRATOS_THROW_ERROR(std::runtime_error,"Undefined settings.contact_model_part in ExplicitSolverStrategy constructor","")
          
          mpFem_model_part                = &(*(settings.fem_model_part));
           if ( mpFem_model_part == NULL )
            KRATOS_THROW_ERROR(std::runtime_error,"Undefined settings.fem_model_part in ExplicitSolverStrategy constructor","")
         
          mpCluster_model_part              = &(*(settings.cluster_model_part));
           if ( mpCluster_model_part == NULL )
            KRATOS_THROW_ERROR(std::runtime_error,"Undefined settings.cluster_model_part in ExplicitSolverStrategy constructor","")
            
          mpInlet_model_part              = &(*(settings.inlet_model_part));
           if ( mpInlet_model_part == NULL )
            KRATOS_THROW_ERROR(std::runtime_error,"Undefined settings.inlet_model_part in ExplicitSolverStrategy constructor","")
                   
      }

      /// Destructor.
      virtual ~ExplicitSolverStrategy()
      {
          Timer::SetOuputFile("TimesPartialRelease");
          Timer::PrintTimingInformation();
      }     
      
      
      struct LessX {
        bool operator()(const SphericParticle* p, const SphericParticle* q) const
        {
          return p->GetGeometry()[0].Coordinates()[0] < q->GetGeometry()[0].Coordinates()[0];
        }
      };
      struct LessY {
        bool operator()(const SphericParticle* p, const SphericParticle* q) const
        {
          return p->GetGeometry()[0].Coordinates()[1] < q->GetGeometry()[0].Coordinates()[1];
        }
      };
      struct LessZ {
        bool operator()(const SphericParticle* p, const SphericParticle* q) const
        {
          return p->GetGeometry()[0].Coordinates()[2] < q->GetGeometry()[0].Coordinates()[2];
        }
      };
      struct SpatialSortingTraits {
        typedef SphericParticle* Point_2;
        typedef LessX Less_x_2;
        typedef LessY Less_y_2;
        typedef LessZ Less_z_2;

        Less_x_2 less_x_2_object() const
        {
          return Less_x_2();
        }
        Less_y_2 less_y_2_object() const
        {
          return Less_y_2();
        }
        Less_z_2 less_z_2_object() const
        {
          return Less_z_2();
        }
      };                        

#ifdef USING_CGAL
      void ReorderParticles(){
          SpatialSortingTraits sst;
          CGAL::spatial_sort(mListOfSphericParticles.begin(), mListOfSphericParticles.end(), sst);          
      }
#endif
            
      template <class T>
      void RebuildListOfSphericParticles(const ElementsArrayType& pElements, std::vector<T*>& rCustomListOfParticles){
          
          KRATOS_TRY
              
          OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), this->GetElementPartition());
          
          rCustomListOfParticles.resize(pElements.size());          
          
          #pragma omp parallel for
          for (int k = 0; k < mNumberOfThreads; k++){
              
            typename ElementsArrayType::const_iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::const_iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
            int elem_number = this->GetElementPartition()[k];

            for (typename ElementsArrayType::const_iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++elem_number){
                T* spheric_particle = dynamic_cast<T*>(&(*particle_pointer_it));
                rCustomListOfParticles[elem_number] = spheric_particle;
            }                                
          }              
          
          return;          
          KRATOS_CATCH("")
          
      }
 
      void RebuildPropertiesProxyPointers(std::vector<SphericParticle*>& rCustomListOfSphericParticles){
          //This function is called for the local mesh and the ghost mesh, so mListOfSphericElements must not be used here.
          KRATOS_TRY         
         
          const int number_of_particles = (int)rCustomListOfSphericParticles.size();
          #pragma omp parallel for 
          for (int i = 0; i < number_of_particles; i++){  
              rCustomListOfSphericParticles[i]->SetFastProperties(mFastProperties);                         
          }                     
          return;            
          KRATOS_CATCH("")
      }
      
      void SendProcessInfoToClustersModelPart(){
          
          ProcessInfo& r_process_info         = mpDem_model_part->GetProcessInfo();
          ProcessInfo& rClustersr_process_info = mpCluster_model_part->GetProcessInfo();
          
          r_process_info[CONTAINS_CLUSTERS] = false;
          rClustersr_process_info[CONTAINS_CLUSTERS] = true;
                  
          rClustersr_process_info[ROTATION_OPTION]     = r_process_info[ROTATION_OPTION];    
          rClustersr_process_info[DELTA_TIME]          = r_process_info[DELTA_TIME];
          rClustersr_process_info[VIRTUAL_MASS_OPTION] = r_process_info[VIRTUAL_MASS_OPTION];
          rClustersr_process_info[TRIHEDRON_OPTION]    = r_process_info[TRIHEDRON_OPTION];
          rClustersr_process_info[NODAL_MASS_COEFF]    = r_process_info[NODAL_MASS_COEFF];                        
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
          const int number_of_particles = (int)rCustomListOfSphericParticles.size();
          #pragma omp parallel for
          for (int i=0; i<number_of_particles; i++){
                            
              int own_properties_id = rCustomListOfSphericParticles[i]->GetProperties().Id();  
              
              for (PropertiesIterator props_it = mpDem_model_part->GetMesh(0).PropertiesBegin(); props_it!= mpDem_model_part->GetMesh(0).PropertiesEnd();   props_it++ ) {
                  int model_part_id = props_it->GetId();
                  if (own_properties_id == model_part_id) {
                      rCustomListOfSphericParticles[i]->SetProperties( *(props_it.base()) );
                      found = true;
                      break;
                  }                  
              }
              
              if (found) continue;
              
              for (PropertiesIterator props_it = mpInlet_model_part->GetMesh(0).PropertiesBegin(); props_it!= mpInlet_model_part->GetMesh(0).PropertiesEnd();   props_it++ ) {
                  int model_part_id = props_it->GetId();
                  if (own_properties_id == model_part_id) {
                      rCustomListOfSphericParticles[i]->SetProperties( *(props_it.base()) );
                      found = true;
                      break;
                  }                  
              }
              
              if (found) continue;
              
              for (PropertiesIterator props_it = mpCluster_model_part->GetMesh(0).PropertiesBegin(); props_it!= mpCluster_model_part->GetMesh(0).PropertiesEnd();   props_it++ ) {
                  int model_part_id = props_it->GetId();
                  if (own_properties_id == model_part_id) {
                      rCustomListOfSphericParticles[i]->SetProperties( *(props_it.base()) );
                      found = true;
                      break;
                  }                  
              }
              
              if (!found) KRATOS_THROW_ERROR(std::logic_error, "This particle could not find its properties!!" , "");
          }                                                                   
      }
      
      virtual void Initialize()
      {
          KRATOS_TRY

          ModelPart& r_model_part            = BaseType::GetModelPart();
          
          ProcessInfo& r_process_info   = r_model_part.GetProcessInfo();
          SendProcessInfoToClustersModelPart();          
                    
          // Omp initializations
     
          mNumberOfThreads = OpenMPUtils::GetNumThreads();
          mNeighbourCounter.resize(mNumberOfThreads);
          
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
          
          CreatePropertiesProxies(mFastProperties, *mpDem_model_part, *mpInlet_model_part, *mpCluster_model_part);
          
          RepairPointersToNormalProperties(mListOfSphericParticles);  // The particles sent to this partition have their own copy of the Kratos properties they were using in the previous partition!!
          RepairPointersToNormalProperties(mListOfGhostSphericParticles);
          
          RebuildPropertiesProxyPointers(mListOfSphericParticles);
          RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);        
                   
          GetBoundingBoxOption()    = r_process_info[BOUNDING_BOX_OPTION];
          GetSearchControl()        = r_process_info[SEARCH_CONTROL];
          GetBoundingBoxStartTime() = r_process_info[BOUNDING_BOX_START_TIME];
          GetBoundingBoxStopTime()  = r_process_info[BOUNDING_BOX_STOP_TIME];
          
          InitializeDEMElements();           
          InitializeFEMElements();
          UpdateMaxIdOfCreatorDestructor();
          InitializeClusters(); // This adds elements to the balls modelpart
          
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
          RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
          
          InitializeSolutionStep();
          
          // ApplyPrescribedBoundaryConditions();
          
          // Search Neighbours and related operations                             
          SetOriginalRadius(r_model_part);

          SearchNeighbours();
          
          if(mDeltaOption == 2) {            
            SetCoordinationNumber(r_model_part);            
          }
          
          ComputeNewNeighboursHistoricalData();  
          
          SearchRigidFaceNeighbours(1); //initial search is performed with hierarchical method in any case MSI
          ComputeNewRigidFaceNeighboursHistoricalData();
          
          //set flag to 2 (search performed this timestep)
          mSearchControl = 2;
          
          // Finding overlapping of initial configurations

          if (r_process_info[CLEAN_INDENT_OPTION]) {
              for (int i = 0; i < 10; i++) CalculateInitialMaxIndentations();
          }

          if (r_process_info[CRITICAL_TIME_OPTION]) {
              InitialTimeStepCalculation();
          }
          
          double& total_friccional_work       = r_process_info[PARTICLE_INELASTIC_FRICTIONAL_WORK];
          total_friccional_work = 0.0;
          
           // 5. Finalize Solution Step.
          //FinalizeSolutionStep();
          //KRATOS_WATCH(r_model_part.GetNodalSolutionStepVariablesList())
                    

      KRATOS_CATCH("")
      }// Initialize()

        virtual void InitializeClusters() {

            KRATOS_TRY

            ElementsArrayType& pElements = mpCluster_model_part->GetCommunicator().LocalMesh().Elements();
            const int number_of_clusters = pElements.size();
            
            //mpParticleCreatorDestructor->FindAndSaveMaxNodeIdInModelPart(*mpDem_model_part); //This has been moved to python main script and checks both dem model part and walls model part (also important!)
            
            #pragma omp for schedule(dynamic, 100) //schedule(guided)
            for (int k = 0; k < number_of_clusters; k++) {

                typename ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                
                cluster_element.Initialize();
                cluster_element.CreateParticles(mpParticleCreatorDestructor.get(), *mpDem_model_part);
            }


            KRATOS_CATCH("")

        }

        virtual void GetClustersForce() {

            KRATOS_TRY
                    
            ProcessInfo& r_process_info    = BaseType::GetModelPart().GetProcessInfo(); //Getting the Process Info of the Balls ModelPart!
            const array_1d<double,3>& gravity = r_process_info[GRAVITY]; 

            ElementsArrayType& pElements = mpCluster_model_part->GetCommunicator().LocalMesh().Elements();  
            const int number_of_clusters = pElements.size();
            
            #pragma omp for schedule(dynamic, 100) //schedule(guided)
            for (int k = 0; k < number_of_clusters; k++) {
                
                typename ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                
                cluster_element.GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES).clear();
                cluster_element.GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT).clear();

                cluster_element.GetClustersForce(gravity);

            } // loop over clusters
                                                            
            KRATOS_CATCH("")                                                                                                            
        }
        
      
      virtual double Solve()
      {
          
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        
        InitializeSolutionStep();
        SearchOperations(r_model_part);
        ForceOperations(r_model_part);
        PerformTimeIntegrationOfMotion();                   
        FinalizeSolutionStep();         

        return 0.00;

        KRATOS_CATCH("")
          
      }//Solve()

      
        void SearchOperations(ModelPart& r_model_part, bool has_mpi = true)
      {
         
         ProcessInfo& r_process_info   = r_model_part.GetProcessInfo();
         int time_step = r_process_info[TIME_STEPS];
         double time = r_process_info[TIME];
         
         if ((time_step + 1) % mNStepSearch == 0 && (time_step > 0)) { //Neighboring search. Every N times. + destruction of particles outside the bounding box                   
              

            if (this->GetBoundingBoxOption() && ((time >= this->GetBoundingBoxStartTime()) && (time <= this->GetBoundingBoxStopTime()))) {
                BoundingBoxUtility();
            }
            
            RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
            RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
              
            SearchNeighbours();

            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles); 
            this->template RebuildListOfSphericParticles <SphericParticle>          (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
            RepairPointersToNormalProperties(mListOfSphericParticles);
            RepairPointersToNormalProperties(mListOfGhostSphericParticles);
            RebuildPropertiesProxyPointers(mListOfSphericParticles);
            RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

            ComputeNewNeighboursHistoricalData();  
            
            SetOriginalRadius(r_model_part);              
            SearchRigidFaceNeighbours(r_process_info[LOCAL_RESOLUTION_METHOD]);
            ComputeNewRigidFaceNeighboursHistoricalData();
            mSearchControl = 2; // Search is active and has been performed during this time step
            //ReorderParticles();
      
            
        }

        else {
          
            mSearchControl = 1; // Search is active but no search has been done this time step;
        }
          
        //RebuildPropertiesProxyPointers(mListOfSphericParticles);
        //RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);        
        
      }//SearchOperations;
      
      void ForceOperations(ModelPart& r_model_part)
      {
        
        // 3. Get and Calculate the forces
        CleanEnergies();
        
        GetForce(); // Basically only calls CalculateRightHandSide( )
        
        //FastGetForce();
        
        GetClustersForce();
        
        if (r_model_part.GetProcessInfo()[COMPUTE_FEM_RESULTS_OPTION]) {
            CalculateConditionsRHSAndAdd();
            CalculateNodalPressuresAndStressesOnWalls();
        }
                  
        // 4. Synchronize (should be just FORCE and TORQUE)
        SynchronizeSolidMesh(r_model_part);
        
      }//ForceOperations;
      
      void InitialTimeStepCalculation()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& r_process_info    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          ElementsIterator it_begin = pElements.ptr_begin();
          ElementsIterator it_end   = pElements.ptr_end();

          double& process_info_delta_time = r_process_info[DELTA_TIME];
          process_info_delta_time         = mMaxTimeStep;
          double temp_time_step           = std::numeric_limits<double>::infinity();
          double elem_critical_time_step  = temp_time_step;

          for (ElementsIterator it = it_begin; it != it_end; it++){
              it->Calculate(DELTA_TIME, elem_critical_time_step, r_process_info);

              if (elem_critical_time_step < temp_time_step){
                  temp_time_step = elem_critical_time_step;
              }

          }

          temp_time_step /= mSafetyFactor;

          if (temp_time_step < mMaxTimeStep) process_info_delta_time = temp_time_step;

          std::cout << std::scientific;
          std::cout << std::setprecision(3) << "************* Using " << process_info_delta_time << " time step. (Critical: "
                    << temp_time_step << " with a diving factor: " << mSafetyFactor <<" ) *************" << "\n" << std::endl;

          KRATOS_CATCH("")
      }

      void GetForce()
      {
          KRATOS_TRY

          ProcessInfo& r_process_info    = BaseType::GetModelPart().GetProcessInfo();
          double dt = r_process_info[DELTA_TIME];
          const array_1d<double,3>& gravity = r_process_info[GRAVITY];             

          const int number_of_particles = (int)mListOfSphericParticles.size();
          
          #pragma omp parallel
          {
              Vector rhs_elem;
              rhs_elem.resize(6);
              #pragma omp for schedule(dynamic, 100) //schedule(guided)
            
              for (int i = 0; i < number_of_particles; i++) {
                              
                  mListOfSphericParticles[i]->CalculateRightHandSide(rhs_elem, r_process_info, dt, gravity, mSearchControl);  
              }
          }

          KRATOS_CATCH("")
      }            
            
      void FastGetForce()
      {
          KRATOS_TRY
          
          ProcessInfo& r_process_info    = BaseType::GetModelPart().GetProcessInfo();
          double dt = r_process_info[DELTA_TIME];
          const array_1d<double,3>& gravity = r_process_info[GRAVITY];
          const int number_of_particles = (int)mListOfSphericParticles.size();

          #pragma omp parallel
          {
              #pragma omp for
              for (int i = 0; i < number_of_particles; i++) { 
                  mListOfSphericParticles[i]->FirstCalculateRightHandSide(r_process_info, dt, mSearchControl); 
              }
              #pragma omp for
              for (int i = 0; i < number_of_particles; i++) { 
                  mListOfSphericParticles[i]->CollectCalculateRightHandSide(r_process_info); 
              }
              #pragma omp for
              for (int i = 0; i < number_of_particles; i++) { 
                  mListOfSphericParticles[i]->FinalCalculateRightHandSide(r_process_info, dt, gravity); 
              }

          }

          KRATOS_CATCH("")
      }

      virtual void PerformTimeIntegrationOfMotion(int StepFlag = 0)
      {
          KRATOS_TRY

          GetScheme()->Calculate(BaseType::GetModelPart(), StepFlag);
          GetScheme()->Calculate(*mpCluster_model_part, StepFlag);

          KRATOS_CATCH("")
      }

      void InitializeSolutionStep()
      {
          KRATOS_TRY

          // SPHERE MODEL PART

          ModelPart& r_model_part      = BaseType::GetModelPart();
          ProcessInfo& r_process_info  = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();                    
          
          OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), this->GetElementPartition());

          #pragma omp parallel for
          for (int k = 0; k < mNumberOfThreads; k++) {
              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
                (it)->InitializeSolutionStep(r_process_info); 
                
              } // loop over particles
          } // loop threads OpenMP

          ApplyPrescribedBoundaryConditions();


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

      void MoveMesh() {}

      void FinalizeSolutionStep()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& r_process_info    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), this->GetElementPartition());

          #pragma omp parallel for
          for (int k = 0; k < mNumberOfThreads; k++) {
              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
                  (it)->FinalizeSolutionStep(r_process_info); //we use this function to call the set initial contacts and the add continuum contacts.
              } //loop over particles

          } // loop threads OpenMP

          KRATOS_CATCH("")
      }

      
    void SetCoordinationNumber(ModelPart& r_model_part) 
    {
        ProcessInfo& r_process_info      = r_model_part.GetProcessInfo();  

        const double in_coordination_number = r_process_info[COORDINATION_NUMBER];
        double out_coordination_number = ComputeCoordinationNumber();
        int iteration = 0;
        int maxiteration = 100;
        double& added_search_distance = r_process_info[SEARCH_TOLERANCE];

        std::cout<<"Setting up Coordination Number by increasing or decreasing the search radius... "<<std::endl;

        if(in_coordination_number <= 0.0) {
            KRATOS_THROW_ERROR(std::runtime_error, "The specified Coordination Number is less or equal to zero, N.C. = ",in_coordination_number)
        }

        while (fabs(out_coordination_number/in_coordination_number - 1.0) > 1e-3) {              
            if (iteration>=maxiteration) break;
            iteration++;
            added_search_distance *= in_coordination_number/out_coordination_number;
            SearchNeighbours(); //r_process_info[SEARCH_TOLERANCE] will be used inside this function, and it's the variable we are updating in this while
            out_coordination_number = ComputeCoordinationNumber();
        }//while

        if (iteration<maxiteration) std::cout<< "Coordination Number iteration converged after "<<iteration<< " iterations, to value " <<out_coordination_number<< " using an extension of " << added_search_distance <<". "<<"\n"<<std::endl;
        else {   
            std::cout << "Coordination Number iteration did NOT converge after "<<iteration<<" iterations. Coordination number reached is "<<out_coordination_number<<". "<<"\n"<<std::endl;
            KRATOS_THROW_ERROR(std::runtime_error,"Please use a Absolute tolerance instead "," ")
            //NOTE: if it doesn't converge, problems occur with contact mesh and rigid face contact.
        }                  
      
    } //SetCoordinationNumber
    
    double ComputeCoordinationNumber()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
         
        unsigned int total_contacts = 0;  
        const int number_of_particles = (int)mListOfSphericParticles.size();
           
        #pragma omp parallel 
        {
            mNeighbourCounter[OpenMPUtils::ThisThread()] = 0;        
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {              
                  mNeighbourCounter[OpenMPUtils::ThisThread()] += mListOfSphericParticles[i]->mNeighbourElements.size();
            }
        }                                      
        for (int i = 0; i < mNumberOfThreads; i++)
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
    
    void InitializeElements()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), this->GetElementPartition());

        #pragma omp parallel for 
        for (int k = 0; k < mNumberOfThreads; k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                (it)->Initialize();
            }

        }
        
        KRATOS_CATCH("")
    }
    
    void InitializeDEMElements() {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ProcessInfo& r_process_info      = r_model_part.GetProcessInfo();
        const int number_of_particles = (int)mListOfSphericParticles.size();
        #pragma omp parallel for
        for(int i=0; i<number_of_particles; i++){
            mListOfSphericParticles[i]->Initialize(r_process_info);
        }
                
        KRATOS_CATCH("")
    }
    
    void InitializeFEMElements()
    {
        KRATOS_TRY

        ConditionsArrayType& pTConditions      = mpFem_model_part->GetCommunicator().LocalMesh().Conditions(); 
        
        OpenMPUtils::CreatePartition(mNumberOfThreads, pTConditions.size(), this->GetElementPartition());

        #pragma omp parallel for 
        for (int k = 0; k < mNumberOfThreads; k++){
            typename ConditionsArrayType::iterator it_begin = pTConditions.ptr_begin() + this->GetElementPartition()[k];
            typename ConditionsArrayType::iterator it_end   = pTConditions.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it){
                (it)->Initialize();
            }
        }

        KRATOS_CATCH("")
    }
    

    void CalculateConditionsRHSAndAdd() {
      
        KRATOS_TRY
      
        ClearFEMForces();
        ConditionsArrayType& pConditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();     
        ProcessInfo& r_process_info = GetFemModelPart().GetProcessInfo();
        
        Vector rhs_cond;
        Vector rhs_cond_elas;
        vector<unsigned int> condition_partition;    
        OpenMPUtils::CreatePartition(mNumberOfThreads, pConditions.size(), condition_partition);    
        unsigned int index;
               
        #pragma omp parallel for private (index, rhs_cond, rhs_cond_elas)
        for (int k = 0; k < mNumberOfThreads; k++) {
            
            typename ConditionsArrayType::iterator it_begin = pConditions.ptr_begin() + condition_partition[k];        
            typename ConditionsArrayType::iterator it_end = pConditions.ptr_begin() + condition_partition[k+1];
            
            for (typename ConditionsArrayType::iterator it = it_begin; it!= it_end; ++it) { //each iteration refers to a different triangle or quadrilateral
                
                Condition::GeometryType& geom = it->GetGeometry();            
                double Element_Area = geom.Area();                     
                it->CalculateRightHandSide(rhs_cond, r_process_info);
                DEMWall* p_wall = dynamic_cast<DEMWall*>(&(*it));
                p_wall->CalculateElasticForces(rhs_cond_elas, r_process_info);
                array_1d<double, 3> Normal_to_Element;
                p_wall->CalculateNormal(Normal_to_Element);                                                                                       
                const unsigned int& dim = geom.WorkingSpaceDimension();
                
                for (unsigned int i = 0; i < geom.size(); i++) { //talking about each of the three nodes of the condition  
                                                                 //we are studying a certain condition here
                    index = i * dim;    //*2;                 
                                                           
                    array_1d<double, 3>& node_rhs      = geom[i].FastGetSolutionStepValue(CONTACT_FORCES);
                    array_1d<double, 3>& node_rhs_elas = geom[i].FastGetSolutionStepValue(ELASTIC_FORCES);
                    array_1d<double, 3>& node_rhs_tang = geom[i].FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
                    double& node_pressure = geom[i].FastGetSolutionStepValue(DEM_PRESSURE);                  
                    double& node_area = geom[i].FastGetSolutionStepValue(DEM_NODAL_AREA);
                    array_1d<double, 3> rhs_cond_comp;
                    
                    geom[i].SetLock(); 
                    
                    for (unsigned int j = 0; j < dim; j++) { //talking about each coordinate x, y and z, loop on them                                           
                        node_rhs[j]      +=  rhs_cond[index+j];
                        node_rhs_elas[j] +=  rhs_cond_elas[index+j];
                        rhs_cond_comp[j] = rhs_cond[index+j];
                    }
                                                           
                    node_area += 0.333333333333333 * Element_Area;
                    
                    //node_pressure actually refers to normal force. Pressure is actually computed later in function Calculate_Nodal_Pressures_and_Stresses()
                    
                    node_pressure += MathUtils<double>::Abs(GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element));
                    
                    noalias(node_rhs_tang) += rhs_cond_comp - GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element) * Normal_to_Element;
                    
                    geom[i].UnSetLock();
                }                                 
            }          
        }

        KRATOS_CATCH("")
    }
    
    
    void ClearFEMForces() {
        
        KRATOS_TRY

        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& pNodes = fem_model_part.Nodes();

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(mNumberOfThreads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for (int k=0; k<mNumberOfThreads; k++) {
            
            typename NodesArrayType::iterator i_begin=pNodes.ptr_begin() + node_partition[k];
            typename NodesArrayType::iterator i_end=pNodes.ptr_begin() + node_partition[k+1];

            for (ModelPart::NodeIterator i=i_begin; i!= i_end; ++i) {
                
                array_1d<double, 3>& node_rhs      = i->FastGetSolutionStepValue(CONTACT_FORCES);
                array_1d<double, 3>& node_rhs_elas = i->FastGetSolutionStepValue(ELASTIC_FORCES);
                array_1d<double, 3>& node_rhs_tang = i->FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
                double& node_pressure              = i->GetSolutionStepValue(DEM_PRESSURE);
                double& node_area                  = i->GetSolutionStepValue(DEM_NODAL_AREA);
                double& shear_stress               = i->FastGetSolutionStepValue(SHEAR_STRESS);
                
                noalias(node_rhs)      = ZeroVector(3);
                noalias(node_rhs_elas) = ZeroVector(3);
                noalias(node_rhs_tang) = ZeroVector(3);
                node_pressure          = 0.0;
                node_area              = 0.0;
                shear_stress           = 0.0;
                
            }    
        }

        KRATOS_CATCH("")
    }
    
    
    void CalculateNodalPressuresAndStressesOnWalls() {
        
        KRATOS_TRY

        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& pNodes = fem_model_part.Nodes();

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(mNumberOfThreads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for (int k = 0; k < mNumberOfThreads; k++) {
            
            typename NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
            typename NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k+1];

            for (ModelPart::NodeIterator i = i_begin; i!= i_end; ++i) {
                
                double& node_pressure = i->FastGetSolutionStepValue(DEM_PRESSURE);
                double& node_area = i->FastGetSolutionStepValue(DEM_NODAL_AREA);
                double& shear_stress = i->FastGetSolutionStepValue(SHEAR_STRESS);
                array_1d<double, 3>& node_rhs_tang = i->FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);

                node_pressure = node_pressure/node_area;
                shear_stress = GeometryFunctions::module(node_rhs_tang)/node_area;
            }
        }

        KRATOS_CATCH("")
    }
      

    void ApplyPrescribedBoundaryConditions() {
      
        KRATOS_TRY
        ProcessInfo& r_process_info    = BaseType::GetModelPart().GetProcessInfo();
        double time = r_process_info[TIME];

        bool impose = false; double vel_x = 0.0;  double vel_y = 0.0;  double vel_z = 0.0;
        double ang_vel_x = 0.0; double ang_vel_y = 0.0; double ang_vel_z = 0.0;

        ModelPart& r_model_part = BaseType::GetModelPart();

        for (ModelPart::MeshesContainerType::iterator mesh_it = r_model_part.GetMeshes().begin(); mesh_it != r_model_part.GetMeshes().end(); ++mesh_it) {
          
            bool has_vel_x = (*mesh_it).Has(IMPOSED_VELOCITY_X_VALUE);
            if(has_vel_x){
                vel_x = (*mesh_it)[IMPOSED_VELOCITY_X_VALUE];
                impose = true;
            }
            
            bool has_vel_y = (*mesh_it).Has(IMPOSED_VELOCITY_Y_VALUE);
            if(has_vel_y){
                vel_y = (*mesh_it)[IMPOSED_VELOCITY_Y_VALUE];
                impose = true;
            }

            bool has_vel_z = (*mesh_it).Has(IMPOSED_VELOCITY_Z_VALUE);
            if(has_vel_z){
                vel_z = (*mesh_it)[IMPOSED_VELOCITY_Z_VALUE];
                impose = true;
            }

            bool has_ang_vel_x = (*mesh_it).Has(IMPOSED_VELOCITY_X_VALUE);
            if(has_ang_vel_x){
                vel_x = (*mesh_it)[IMPOSED_VELOCITY_X_VALUE];
                impose = true;
            }

            bool has_ang_vel_y = (*mesh_it).Has(IMPOSED_VELOCITY_Y_VALUE);
            if(has_ang_vel_y){
                vel_y = (*mesh_it)[IMPOSED_VELOCITY_Y_VALUE];
                impose = true;
            }

            bool has_ang_vel_z = (*mesh_it).Has(IMPOSED_VELOCITY_Z_VALUE);
            if(has_ang_vel_z){
                vel_z = (*mesh_it)[IMPOSED_VELOCITY_Z_VALUE];
                impose = true;
            }

            double vel_start = (*mesh_it)[VELOCITY_START_TIME];
            double vel_stop = (*mesh_it)[VELOCITY_STOP_TIME];
            if(time<vel_start || time>vel_stop){
                impose = false;
            }

            NodesArrayType& pNodes = mesh_it->Nodes();
            vector<unsigned int> node_partition;
            OpenMPUtils::CreatePartition(mNumberOfThreads, pNodes.size(), node_partition);

            #pragma omp parallel for

            for (int k=0; k<mNumberOfThreads; k++) {

                typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
                typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

                for (ModelPart::NodeIterator i=i_begin; i!= i_end; ++i) {

                      i->Set(DEMFlags::FIXED_VEL_X, false);
                      i->Set(DEMFlags::FIXED_VEL_Y, false);
                      i->Set(DEMFlags::FIXED_VEL_Z, false);
                      i->Set(DEMFlags::FIXED_ANG_VEL_X, false);
                      i->Set(DEMFlags::FIXED_ANG_VEL_Y, false);
                      i->Set(DEMFlags::FIXED_ANG_VEL_Z, false);

                } // loop over particles
            } // loop threads OpenMP

            if(impose){
              #pragma omp parallel for
            
              for (int k=0; k<mNumberOfThreads; k++) {
                      
                  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
                  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

                  for (ModelPart::NodeIterator i=i_begin; i!= i_end; ++i) {
                      
                      array_1d<double, 3>& velocity = i->FastGetSolutionStepValue(VELOCITY);
                      array_1d<double, 3>& angular_velocity = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);

                      if(has_vel_x){
                        velocity[0] = vel_x;
                        i->Set(DEMFlags::FIXED_VEL_X, true);}

                      if(has_vel_y){
                        velocity[1] = vel_y;
                        i->Set(DEMFlags::FIXED_VEL_Y, true);}
                      
                      if(has_vel_z){
                        velocity[2] = vel_z;
                        i->Set(DEMFlags::FIXED_VEL_Z, true);}

                      if(has_ang_vel_x){
                        angular_velocity[0] = ang_vel_x;
                        i->Set(DEMFlags::FIXED_ANG_VEL_X, true);}

                      if(has_ang_vel_y){
                        angular_velocity[1] = ang_vel_y;
                        i->Set(DEMFlags::FIXED_ANG_VEL_Y, true);}

                      if(has_ang_vel_z){
                        angular_velocity[2] = ang_vel_z;
                        i->Set(DEMFlags::FIXED_ANG_VEL_Z, true);}

                  }     // loop over particles
              }         // loop threads OpenMP
            }           // If there is some value to be imposed
        }               // for each mesh
      
        KRATOS_CATCH("")
    }
    
    void SetSearchRadius(ModelPart& r_model_part, const double added_search_distance=0.0, const double amplification=1.0)
    {
        KRATOS_TRY
            
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        
        mNumberOfElementsOldRadiusList = number_of_elements;
        
        this->GetRadius().resize(number_of_elements);                

	#pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++ ){

            this->GetRadius()[i] = amplification*(added_search_distance + mListOfSphericParticles[i]->GetRadius());

        }

        KRATOS_CATCH("")
    }
    
       
      
    void SetOriginalRadius(ModelPart& r_model_part)
    {
        KRATOS_TRY
     
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        
        this->GetOriginalRadius().resize(number_of_elements);

        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++ ){
        
            SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin() + i;
 
            this->GetOriginalRadius()[i] = particle_pointer_it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);

        }

        KRATOS_CATCH("")
    }

    virtual void SearchNeighbours(const double amplification = 1.0)
    {
        KRATOS_TRY

        ModelPart& r_model_part            = BaseType::GetModelPart();
        ProcessInfo& r_process_info   = r_model_part.GetProcessInfo();
        const double added_search_distance = r_process_info[SEARCH_TOLERANCE];        
        
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        if(!number_of_elements) return; 

        SetSearchRadius(r_model_part, added_search_distance, amplification);
        
        r_model_part.GetCommunicator().GhostMesh().ElementsArray().clear();
        
        GetResults().resize(number_of_elements);
        GetResultsDistances().resize(number_of_elements);                              
        
        mpSpSearch->SearchElementsInRadiusExclusive(r_model_part, this->GetRadius(), this->GetResults(), this->GetResultsDistances());
        const int number_of_particles = (int)mListOfSphericParticles.size();
        
        #pragma omp for schedule(dynamic, 100) //schedule(guided)
        for (int i=0; i<number_of_particles; i++){
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
    
    virtual void ComputeNewNeighboursHistoricalData() {
        
        KRATOS_TRY  
        const int number_of_particles = (int)mListOfSphericParticles.size();

        #pragma omp parallel
        {
            std::vector<int> mTempNeighboursIds;
            std::vector<array_1d<double, 3> > mTempNeighbourElasticContactForces;
            std::vector<array_1d<double, 3> > mTempNeighbourTotalContactForces;
        
            #pragma omp for
            for(int i=0; i<number_of_particles; i++){
                mListOfSphericParticles[i]->ComputeNewNeighboursHistoricalData(mTempNeighboursIds,
                                                                               mTempNeighbourElasticContactForces,
                                                                               mTempNeighbourTotalContactForces);
            }
        }

        KRATOS_CATCH("")
    }
    
     void ComputeNewRigidFaceNeighboursHistoricalData()
    {
        KRATOS_TRY
        const int number_of_particles = (int)mListOfSphericParticles.size();

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++){
            mListOfSphericParticles[i]->ComputeNewRigidFaceNeighboursHistoricalData();
        }

        KRATOS_CATCH("")
    }
    
    virtual void SearchRigidFaceNeighbours(int local_resolution_method)
    {
        KRATOS_TRY
        
        ElementsArrayType&   pElements    = mpDem_model_part->GetCommunicator().LocalMesh().Elements();    
        ConditionsArrayType& pTConditions = mpFem_model_part->GetCommunicator().LocalMesh().Conditions(); 

        if (pTConditions.size() > 0) {
        
            const int number_of_particles = (int)mListOfSphericParticles.size();
            
            this->GetRigidFaceResults().resize(number_of_particles);
            this->GetRigidFaceResultsDistances().resize(number_of_particles);
             
            #pragma omp parallel for 
            for (int i = 0; i < number_of_particles; i++){
                mListOfSphericParticles[i]->mNeighbourRigidFaces.resize(0);
                mListOfSphericParticles[i]->mNeighbourRigidFacesPram.resize(0);
            }
            mpDemFemSearch->SearchRigidFaceForDEMInRadiusExclusiveImplementation(pElements, pTConditions, this->GetOriginalRadius(), this->GetRigidFaceResults(), this->GetRigidFaceResultsDistances());
            
            #pragma omp parallel for
            for (int i = 0; i < number_of_particles; i++ ){
              std::vector<DEMWall*>& neighbour_rigid_faces = mListOfSphericParticles[i]->mNeighbourRigidFaces;
              for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[i].begin(); neighbour_it != this->GetRigidFaceResults()[i].end(); ++neighbour_it){                                                  
                    Condition* p_neighbour_condition = (*neighbour_it).get();
                    DEMWall* p_wall = dynamic_cast<DEMWall*>( p_neighbour_condition );
                    RigidFaceGeometricalConfigureType::DoubleHierarchyMethod(mListOfSphericParticles[i], p_wall); 
                }
              
                std::vector<double>& RF_Param = mListOfSphericParticles[i]->mNeighbourRigidFacesPram;
                size_t neigh_size = neighbour_rigid_faces.size();
                std::vector<DEMWall*> temporal_neigh (0); 
                std::vector<double>  RF_Param_temp (0);

                for (unsigned int n = 0; n < neigh_size; n++ ){

                  if(RF_Param[n * 16 + 15] != -1) //if(it is not a -1 contact neighbour, we copy it)
                  {
                    temporal_neigh.push_back(neighbour_rigid_faces[n]);
                    int current_RF_temp_size = RF_Param_temp.size();
                    RF_Param_temp.resize(current_RF_temp_size + 16);
                    for (int k=0; k<16; k++ )
                    {
                      RF_Param_temp[current_RF_temp_size + k] = RF_Param[n * 16 + k];
                    }

                  }//if(it is not a -1 contact neighbour, we copy it)

                }//loop over temporal neighbours

                //swap
                temporal_neigh.swap(neighbour_rigid_faces);
                RF_Param_temp.swap(RF_Param);

                this->GetRigidFaceResults()[i].clear();
                this->GetRigidFaceResultsDistances()[i].clear();   
            }
                                     
            //typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
            const int number_of_conditions = (int)pTConditions.size();
           
            #pragma omp parallel 
            {              
                #pragma omp for
                for (int i = 0; i < number_of_conditions; i++)
                {               
                    ConditionsArrayType::iterator ic = pTConditions.begin() + i;
                    DEMWall* wall = dynamic_cast<Kratos::DEMWall*>( &(*ic) );
                    wall->mNeighbourSphericParticles.resize(0);
                }       

                #pragma omp for
                for (int i = 0; i < number_of_particles; i++ ){                                  
                    for (unsigned int j = 0; j < mListOfSphericParticles[i]->mNeighbourRigidFaces.size(); j++) 
                    {
                        DEMWall* p_wall = mListOfSphericParticles[i]->mNeighbourRigidFaces[j];
                        #pragma omp critical
                        {
                            p_wall->mNeighbourSphericParticles.push_back(mListOfSphericParticles[i]);                                                                                      
                        }
                    }
                 }
              
              }//end parallel      
              
        }
       
                    
        KRATOS_CATCH("")
    }
    
                   
    /* This should work only with one iteration, but it with mpi does not */
    void CalculateInitialMaxIndentations()
    {
 
        KRATOS_TRY
        std::vector<double> indentations_list, indentations_list_ghost;
        indentations_list.resize(mListOfSphericParticles.size());
        indentations_list_ghost.resize(mListOfGhostSphericParticles.size());
        
        const int number_of_particles = (int)mListOfSphericParticles.size();
               
        #pragma omp parallel for
        for( int i=0; i<number_of_particles; i++ ){
            double indentation;
            mListOfSphericParticles[i]->CalculateMaxBallToBallIndentation(indentation);
            double max_indentation = std::max(0.0, 0.5 * indentation); // reducing the radius by half the indentation is enough
 
            mListOfSphericParticles[i]->CalculateMaxBallToFaceIndentation(indentation);
            max_indentation = std::max(max_indentation, indentation);
            indentations_list[i] = max_indentation;
        }

        #pragma omp parallel for //THESE TWO LOOPS CANNOT BE JOINED, BECAUSE THE RADII ARE CHANGING.
        for( int i=0; i<number_of_particles; i++ ){
            mListOfSphericParticles[i]->SetInteractionRadius(mListOfSphericParticles[i]->GetInteractionRadius() - indentations_list[i]);
        }

        SynchronizeSolidMesh(BaseType::GetModelPart());
        const int number_of_ghost_particles = (int)mListOfGhostSphericParticles.size();

        #pragma omp parallel for //THESE TWO LOOPS CANNOT BE JOINED, BECAUSE THE RADII ARE CHANGING.
        for( int i=0; i<number_of_ghost_particles; i++ ){
            mListOfGhostSphericParticles[i]->SetInteractionRadius(mListOfGhostSphericParticles[i]->GetInteractionRadius() - indentations_list_ghost[i]);
        }
       
        #pragma omp parallel for
        for( int i=0; i<number_of_particles; i++ ){
            double indentation;
            mListOfSphericParticles[i]->CalculateMaxBallToBallIndentation(indentation);
        }
 
        KRATOS_CATCH("")
   
    } // CalculateInitialMaxIndentations()

    void PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part)
    {
        mcontacts_model_part.GetCommunicator().SetNumberOfColors(r_model_part.GetCommunicator().GetNumberOfColors());
        mcontacts_model_part.GetCommunicator().NeighbourIndices() = r_model_part.GetCommunicator().NeighbourIndices();
    }
    
    void PrepareElementsForPrinting() {
       
        ProcessInfo& r_process_info    = (*mpDem_model_part).GetProcessInfo();
        ElementsArrayType& pElements        = (*mpDem_model_part).GetCommunicator().LocalMesh().Elements();
               
        vector<unsigned int> element_partition;

        OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), element_partition);

        #pragma omp parallel for

        for(int k=0; k<mNumberOfThreads; k++)
        { 
            typename ElementsArrayType::iterator it_begin=pElements.ptr_begin() + element_partition[k];
            typename ElementsArrayType::iterator it_end=pElements.ptr_begin() + element_partition[k+1];

            for (typename ElementsArrayType::iterator it= it_begin; it!=it_end; ++it) {
                Element* raw_p_element = &(*it);
                SphericParticle* p_sphere = dynamic_cast<SphericParticle*>( raw_p_element );    
                p_sphere->PrepareForPrinting(r_process_info);
            } //loop over ELEMENTS
        }// loop threads OpenMP        
    } //PrepareElementsForPrinting

    void SynchronizeSolidMesh(ModelPart& r_model_part)
    {
        r_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
//         r_model_part.GetCommunicator().SynchronizeDofs();
    }
            
    void CleanEnergies()
     {
          KRATOS_TRY

          ProcessInfo& r_process_info    = BaseType::GetModelPart().GetProcessInfo();
                  
          double& total_elastic_energy        = r_process_info[PARTICLE_ELASTIC_ENERGY];
          double& total_damping_energy        = r_process_info[PARTICLE_INELASTIC_VISCODAMPING_ENERGY];

          total_elastic_energy        = 0.0;
          total_damping_energy        = 0.0;
          //total_friccional_work       = 0.0;

          KRATOS_CATCH("")
     }            
    

     void GlobalDamping() {   // flagged for deletion

         KRATOS_TRY

         ModelPart& r_model_part = BaseType::GetModelPart();
         ElementsArrayType& pElements = GetElements(r_model_part);

         OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), this->GetElementPartition());

         #pragma omp parallel for
         for (int k = 0; k < mNumberOfThreads; k++) {
             typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
             typename ElementsArrayType::iterator it_end = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

             for (typename ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
                 ModelPart::NodeType pNode = it->GetGeometry()[0];

                 array_1d<double, 3>& total_force = pNode.FastGetSolutionStepValue(TOTAL_FORCES); //Includes all elastic, damping, but not external (gravity)
                 array_1d<double, 3>& velocity = pNode.FastGetSolutionStepValue(VELOCITY);

                 const double global_damping         = 0.0;

                 if (pNode.IsNot(DEMFlags::FIXED_VEL_X)) {
                     total_force[0] = total_force[0] - global_damping * fabs(total_force[0]) * GeometryFunctions::sign(velocity[0]);
                 }

                 if (pNode.IsNot(DEMFlags::FIXED_VEL_Y)) {
                     total_force[1] = total_force[1] - global_damping * fabs(total_force[1]) * GeometryFunctions::sign(velocity[1]);
                 }

                 if (pNode.IsNot(DEMFlags::FIXED_VEL_Z)) {
                     total_force[2] = total_force[2] - global_damping * fabs(total_force[2]) * GeometryFunctions::sign(velocity[2]);
                 }

             } // loop over particles

         }   // loop threads OpenMP

         KRATOS_CATCH("")
     }       // GlobalDamping


    
    //*******************************************************************************************************
    
    
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
    
    int&                                         GetNStepSearch(){return (mNStepSearch);}
    int&                                         GetSearchControl(){return mSearchControl;}
    int&                                         GetBoundingBoxOption(){return (mBoundingBoxOption);}
    int&                                         GetNumberOfThreads(){return (mNumberOfThreads);}
    double&                                      GetBoundingBoxStartTime(){return (mBoundingBoxStartTime);}
    double&                                      GetBoundingBoxStopTime(){return (mBoundingBoxStopTime);}
    
    double&                                      GetMaxTimeStep(){return (mMaxTimeStep);}
    double&                                      GetSafetyFactor(){return (mSafetyFactor);}
    
    int&                                         GetDeltaOption(){return (mDeltaOption);}
    vector<unsigned int>&                        GetNeighbourCounter(){return(mNeighbourCounter);}
                                           
    int&                                         GetNumberOfElementsOldRadiusList(){return (mNumberOfElementsOldRadiusList);}

    vector<unsigned int>&                        GetElementPartition(){return (mElementPartition);}

    typename ParticleCreatorDestructor::Pointer& GetParticleCreatorDestructor(){return (mpParticleCreatorDestructor);}
    typename DEMIntegrationScheme::Pointer&         GetScheme(){return (mpScheme);}
    typename SpatialSearch::Pointer&             GetSpSearch(){return (mpSpSearch);}
    
    //Cfeng
    VectorResultConditionsContainerType&         GetRigidFaceResults(){return(mRigidFaceResults);}
    VectorDistanceType&                          GetRigidFaceResultsDistances(){return(mRigidFaceResultsDistances);}
    vector<unsigned int>&                        GetConditionPartition(){return (mConditionPartition);}
    
    std::vector<PropertiesProxy>                 mFastProperties;

    virtual ElementsArrayType& GetElements(ModelPart& r_model_part)
    {
        return r_model_part.GetCommunicator().LocalMesh().Elements();
    }
    
    protected:

    // Member variables

    VectorResultElementsContainerType            mResults;
    VectorDistanceType                           mResultsDistances;
    RadiusArrayType                              mRadius;
    RadiusArrayType                              mOriginalRadius;

    int                                          mNStepSearch;
    int                                          mSearchControl;
    int                                          mBoundingBoxOption;
    double                                       mBoundingBoxStartTime;
    double                                       mBoundingBoxStopTime;
    int                                          mNumberOfThreads;
    int                                          mNumberOfElementsOldRadiusList;

    double                                       mMaxTimeStep;
    double                                       mSafetyFactor;
    
    int                                          mDeltaOption;
    vector<unsigned int>                         mNeighbourCounter;
          
    vector<unsigned int>                         mElementPartition;
    typename ParticleCreatorDestructor::Pointer  mpParticleCreatorDestructor;
    typename DEM_FEM_Search::Pointer             mpDemFemSearch;
    typename DEMIntegrationScheme::Pointer          mpScheme;
    typename SpatialSearch::Pointer              mpSpSearch;
    
    
    ////Cfeng
    VectorResultConditionsContainerType  mRigidFaceResults;
    VectorDistanceType                   mRigidFaceResultsDistances;
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

