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

#include "custom_elements/spheric_swimming_particle.h"
#include "custom_elements/Particle_Contact_Element.h"

#include "includes/variables.h"
#include "DEM_application.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#define CUSTOMTIMER 1  // ACTIVATES AND DISABLES ::TIMER:::::

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/integration_scheme.h"
#include "custom_utilities/create_and_destroy.h"

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
 

      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      ExplicitSolverStrategy(){}

      ExplicitSolverStrategy(ModelPart& r_model_part,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool move_mesh_flag,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename IntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch
      ): BaseType(r_model_part, move_mesh_flag)
      {

          mElementsAreInitialized      = false;
          mInitializeWasPerformed      = false;
          mpParticleCreatorDestructor  = p_creator_destructor;
          mpScheme                     = pScheme;
          mpSpSearch                   = pSpSearch;
          mMaxTimeStep                 = max_delta_time;
          mNStepSearch                 = n_step_search;
          mSafetyFactor                = safety_factor;
      }

      /// Destructor.
      virtual ~ExplicitSolverStrategy()
      {
          Timer::SetOuputFile("TimesPartialRelease");
          Timer::PrintTimingInformation();
      }

      virtual void Initialize()
      {
          KRATOS_TRY
          
          KRATOS_TIMER_START("INITIALIZE")
          
          ModelPart& r_model_part            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();

          int NumberOfElements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

          this->GetResults().resize(NumberOfElements);
          this->GetResultsDistances().resize(NumberOfElements);
          this->GetRadius().resize(NumberOfElements);

          // Omp initializations
          this->GetNumberOfThreads() = OpenMPUtils::GetNumThreads();

          // 0. Set search radius
          
          
          SetSearchRadius(r_model_part, rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION]);

          // 1. Search Neighbours with tolerance (Not in mpi.)
          this->GetBoundingBoxOption() = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

          // 2. Search Neighbours with tolerance (after first repartition process)
          SearchNeighbours(r_model_part);

          // 3. Finding overlapping of initial configurations

          if (rCurrentProcessInfo[CLEAN_INDENT_OPTION]){
              CalculateInitialMaxIndentations();
          }

          // 4. Initializing elements and perform the repartition
          if (!mElementsAreInitialized){
              InitializeElements();
          }

          mInitializeWasPerformed = true;
          
           // 5. Finalize Solution Step.
          FinalizeSolutionStep();
          

          KRATOS_TIMER_STOP("INITIALIZE")

      KRATOS_CATCH("")
      }// Initialize()

      virtual double Solve()
      {
          KRATOS_TRY

          KRATOS_TIMER_START("BEGIN")
          ModelPart& r_model_part            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();

          int time_step = rCurrentProcessInfo[TIME_STEPS];
          KRATOS_TIMER_STOP("BEGIN")
          
          // 1. Here we initialize member variables that depend on the rCurrentProcessInfo
          KRATOS_TIMER_START("InitializeSolutionStep")
          InitializeSolutionStep();
                
          // 2. Get and Calculate the forces
          KRATOS_TIMER_START("GetForce")
          GetForce();
                    
          // 3. Motion Integration
          KRATOS_TIMER_START("PerformTimeIntegrationOfMotion")
          PerformTimeIntegrationOfMotion(rCurrentProcessInfo); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot
          
          // 4. Synchronize
          KRATOS_TIMER_START("SynchronizeSolidMesh")
          SynchronizeSolidMesh(r_model_part);
                  
          // 5. Neighbouring search. Every N times. + destruction of particles outside the bounding box                   
          
          if ((time_step + 1) % mNStepSearch == 0 && time_step > 0){
              if (this->GetBoundingBoxOption()){
                  BoundingBoxUtility();
              }

              SearchNeighbours(r_model_part);
          }
          
          KRATOS_TIMER_STOP("SearchNeighbours")

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
          KRATOS_WATCH(mMaxTimeStep)

          std::cout<< "****************** Calculated time step is " << temp_time_step << " ******************" << std::endl;

          KRATOS_CATCH("")
      }

      void GetForce()
      {
          KRATOS_TRY

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

          Vector rhs_elem;

          #pragma omp parallel for private(rhs_elem)
          for (int k = 0; k < this->GetNumberOfThreads(); k++){

              if (rhs_elem.size() != 6){
                rhs_elem.resize(6);
              }

              typename ElementsArrayType::iterator it_begin   = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end     = pElements.ptr_begin() + this->GetElementPartition()[k+1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  Element::GeometryType& geom = it->GetGeometry();

                  (it)->CalculateRightHandSide(rhs_elem, rCurrentProcessInfo);

                  array_1d<double,3>& applied_force  = geom(0)->FastGetSolutionStepValue(TOTAL_FORCES);
                  array_1d<double,3>& applied_moment = geom(0)->FastGetSolutionStepValue(PARTICLE_MOMENT);

                  for (int i = 0; i < 3; i++){
                      applied_force[i]  = rhs_elem[i];
                      applied_moment[i] = rhs_elem[3 + i];
                  }

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
                (it)->InitializeSolutionStep(rCurrentProcessInfo); // we use this function to call the set initial contacts and the add continuum contacts.        
              } // loop over particles

          } // loop threads OpenMP

        KRATOS_CATCH("")
      }

      void BoundingBoxUtility()
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

      void FixVelocities()
      {
          KRATOS_TRY

          KRATOS_WATCH("")
          KRATOS_WATCH("FIXING VELOCITIES!")

          ModelPart& r_model_part             = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo    = r_model_part.GetProcessInfo();
          ElementsArrayType& pElements        = r_model_part.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

          #pragma omp parallel for

          for (int k = 0; k < this->GetNumberOfThreads(); k++){
              typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){

                  if (it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 1){
                      (it)->GetGeometry()(0)->Fix(VELOCITY_Y);
                      (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y) = rCurrentProcessInfo[FIXED_VEL_TOP];
                  }

                  if (it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 2){
                      (it)->GetGeometry()(0)->Fix(VELOCITY_Y);
                      (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y) = rCurrentProcessInfo[FIXED_VEL_BOT];
                  }

              } // for over particles

          } // for threads OpenMP

          KRATOS_CATCH("")
      }

      void FreeVelocities()
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

                  if (it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 1){
                    (it)->GetGeometry()(0)->Free(VELOCITY_Y);
                    rCurrentProcessInfo[FIXED_VEL_TOP] = (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y); //cutre way yeah!
                    //I only store one value for every ball in the group ID
                    (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
                  }

                  if (it->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) == 2){
                    (it)->GetGeometry()(0)->Free(VELOCITY_Y);
                    rCurrentProcessInfo[FIXED_VEL_BOT] = (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y); //cutre way yeah!
                    //I only store one value for every ball in the group ID
                    (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
                  }

              } // loop over particles

          } // loop threads OpenMP

          KRATOS_CATCH("")
      }


    void CalculateEnergies(){}


    void InitializeElements()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo      = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        int trihedron_OPTION = rCurrentProcessInfo[TRIHEDRON_OPTION];

        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for //private(index, MassMatrix)  //M. proba de compilar sense mass matrix??
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
            //  Element::GeometryType& geom = it->GetGeometry(); ///WARNING: COMMENTED AVOIDING WARNING COMPILATION
                (it)->Initialize();
              // 4. Set the Local Initial Axes for the trihedron Option

                if (trihedron_OPTION == 1){
                  double dummy = 0.0;
                  (it)->Calculate(DUMMY_LOCAL_AXES, dummy, rCurrentProcessInfo);
                }

            }

        }

        mElementsAreInitialized = true;

        KRATOS_CATCH("")
    }

    void SetSearchRadius(ModelPart& r_model_part, double radiusExtend)
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        //ProcessInfo& rCurrentProcessInfo      = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it){
            this->GetRadius()[particle_pointer_it - pElements.begin()] = (1.0 + radiusExtend) * particle_pointer_it->GetGeometry()(0)->GetSolutionStepValue(RADIUS); //if this is changed, then compobation before adding neighbours must change also.

        }

        KRATOS_CATCH("")
    }

    void KRATOS_CHECK_SIZE(const char* msg, int id)
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it){

            if (particle_pointer_it->Id() == id){
                std::cout << msg << " " << particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS).size() << std::endl; break;
            }
        }

        KRATOS_CATCH("")
    }

    virtual void SearchNeighbours(ModelPart& r_model_part)
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        
        for (SpatialSearch::ElementsContainerType::iterator i = pElements.begin(); i != pElements.end(); i++){
            i->GetValue(NEIGHBOUR_ELEMENTS).clear();
        }
        
        

        mpSpSearch->SearchElementsInRadiusExclusive(r_model_part, this->GetRadius(), this->GetResults(), this->GetResultsDistances());
//         mpSpSearch->SearchElementsInRadiusExclusive(r_model_part,this->GetRadius(),this->GetResults());
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
        
        #pragma omp parallel for
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            size_t ResultCounter = this->GetElementPartition()[k];

            for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++ResultCounter){

                for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[ResultCounter].begin(); neighbour_it != this->GetResults()[ResultCounter].end(); ++neighbour_it){
                    particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);         
                }

                this->GetResults()[ResultCounter].clear();
                this->GetResultsDistances()[ResultCounter].clear();
            }

        }

        KRATOS_CATCH("")
    }

    void CalculateInitialMaxIndentations(){

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

                if (!(it->GetGeometry()(0)->pGetDof(VELOCITY_X)->IsFixed())){
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
        r_model_part.GetCommunicator().SynchronizeDofs();
    }
    
    // Getting member variables

    VectorResultElementsContainerType&           GetResults(){return(mResults);}
    VectorDistanceType&                          GetResultsDistances(){return(mResultsDistances);}
    RadiusArrayType&                             GetRadius(){return(mRadius);}

    bool&                                        GetElementsAreInitialized(){return (mElementsAreInitialized);}
    bool&                                        GetInitializeWasPerformed(){return (mInitializeWasPerformed);}
    
    int&                                         GetNStepSearch(){return (mNStepSearch);}
    int&                                         GetBoundingBoxOption(){return (mBoundingBoxOption);}
    int&                                         GetNumberOfThreads(){return (mNumberOfThreads);}

    double&                                      GetMaxTimeStep(){return (mMaxTimeStep);}
    double&                                      GetSafetyFactor(){return (mSafetyFactor);}

    vector<unsigned int>&                        GetElementPartition(){return (mElementPartition);}

    typename ParticleCreatorDestructor::Pointer& GetParticleCreatorDestructor(){return (mpParticleCreatorDestructor);}
    typename IntegrationScheme::Pointer&         GetScheme(){return (mpScheme);}
    typename SpatialSearch::Pointer&             GetSpSearch(){return (mpSpSearch);}
    
    protected:

    // Member variables

    VectorResultElementsContainerType            mResults;
    VectorDistanceType                           mResultsDistances;
    RadiusArrayType                              mRadius;

    bool                                         mElementsAreInitialized;
    bool                                         mInitializeWasPerformed;

    int                                          mNStepSearch;
    int                                          mBoundingBoxOption;
    int                                          mNumberOfThreads;

    double                                       mMaxTimeStep;
    double                                       mSafetyFactor;

    vector<unsigned int>                         mElementPartition;
    typename ParticleCreatorDestructor::Pointer  mpParticleCreatorDestructor;
    typename IntegrationScheme::Pointer          mpScheme;
    typename SpatialSearch::Pointer              mpSpSearch;
    
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

