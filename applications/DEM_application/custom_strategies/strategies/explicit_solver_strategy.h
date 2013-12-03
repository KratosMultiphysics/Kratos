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

  bool compare_ids(const boost::weak_ptr<Element>& i, const boost::weak_ptr<Element>& j) {
    
     if( (i.lock())->Id() <  (j.lock())->Id()) return true;
     if( (i.lock())->Id() >  (j.lock())->Id()) return false;
  
     std::cout<<"Two elements with the same Id!! problems can occur here!!"<<std::endl<<std::flush;           
     return false;               
     
         
  }  
    
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
          mDeltaOption                   = delta_option,
          mSearchTolerance               = search_tolerance,
          mCoordinationNumber            = coordination_number;
          mpParticleCreatorDestructor    = p_creator_destructor;
          mpScheme                       = pScheme;
          mpSpSearch                     = pSpSearch;
          mMaxTimeStep                   = max_delta_time;
          mNStepSearch                   = n_step_search;
          mSafetyFactor                  = safety_factor;
          mNumberOfElementsOldRadiusList = 0;
		  
	  //Cfengs
	  mpDem_model_part             = &r_model_part;
	  mpFem_model_part             = &fem_model_part;
		  
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

          ModelPart& r_model_part            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();

          int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

          this->GetResults().resize(number_of_elements);
          this->GetResultsDistances().resize(number_of_elements);          
		  
		  //Cfeng
		  this->GetRigidFaceResults().resize(number_of_elements);
		  this->GetRigidFaceResultsDistances().resize(number_of_elements);

          // Omp initializations
          this->GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
          mNeighbourCounter.resize(this->GetNumberOfThreads());

          // 0. Set search radius
                    
          SetSearchRadius(r_model_part, 1.0);

          this->GetBoundingBoxOption() = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

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

          // 4. Initializing elements and perform the repartition
          //if (!mElementsAreInitialized){
            InitializeSolutionStep();
            InitializeElements();
              
          //}

          mInitializeWasPerformed = true;
          
           // 5. Finalize Solution Step.
          FinalizeSolutionStep();

      KRATOS_CATCH("")
      }// Initialize()

      virtual double Solve()
      {
          KRATOS_TRY

          ModelPart& r_model_part            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo   = r_model_part.GetProcessInfo();                    

          int NumberOfElements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

          this->GetResults().resize(NumberOfElements);
          this->GetResultsDistances().resize(NumberOfElements);
          this->GetRadius().resize(NumberOfElements);          
          
          int time_step = rCurrentProcessInfo[TIME_STEPS];
          
          // 1. Here we initialize member variables that depend on the rCurrentProcessInfo
          //InitializeSolutionStep();
		  
		 // 2. Neighbouring search. Every N times. + destruction of particles outside the bounding box                   
	
          if ((time_step + 1) % mNStepSearch == 0 && time_step > 0){
              if (this->GetBoundingBoxOption()){
                  BoundingBoxUtility();
              }

              SearchNeighbours();
              ComputeNewNeighboursHistoricalData();  
              
              ///Cfeng RigidFace search
              SearchRigidFaceNeighbours();
              ComputeNewRigidFaceNeighboursHistoricalData();
              
          }

                
          // 3. Get and Calculate the forces
          GetForce();
               
          // 4. Motion Integration
          PerformTimeIntegrationOfMotion(rCurrentProcessInfo); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot

           ////Cfeng, compute rigid face movement
		  Compute_RigidFace_Movement();
		  

          // 5. Synchronize
          SynchronizeSolidMesh(r_model_part);
	          
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

          std::cout<< "****************** Calculated time step is " << temp_time_step << " ******************" << "\n" << std::endl;

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
                rhs_elem.resize(6,false);
              }

              typename ElementsArrayType::iterator it_begin   = pElements.ptr_begin() + this->GetElementPartition()[k];
              typename ElementsArrayType::iterator it_end     = pElements.ptr_begin() + this->GetElementPartition()[k+1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  Element::GeometryType& geom = it->GetGeometry();

                  (it)->CalculateRightHandSide(rhs_elem, rCurrentProcessInfo);
				  

                  array_1d<double,3>& total_forces  = geom(0)->FastGetSolutionStepValue(TOTAL_FORCES);
                  array_1d<double,3>& total_moment = geom(0)->FastGetSolutionStepValue(PARTICLE_MOMENT);

                  for (int i = 0; i < 3; i++){
                      total_forces[i] = rhs_elem[i];
                      total_moment[i] = rhs_elem[3 + i];
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

                  if (it->GetGeometry()(0)->FastGetSolutionStepValue(GROUP_ID) == 1){
                      (it)->GetGeometry()(0)->Fix(VELOCITY_Y);
                      (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y) = rCurrentProcessInfo[FIXED_VEL_TOP];
                  }

                  if (it->GetGeometry()(0)->FastGetSolutionStepValue(GROUP_ID) == 2){
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

                  if (it->GetGeometry()(0)->FastGetSolutionStepValue(GROUP_ID) == 1){
                    (it)->GetGeometry()(0)->Free(VELOCITY_Y);
                    rCurrentProcessInfo[FIXED_VEL_TOP] = (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y); //cutre way yeah!
                    //I only store one value for every ball in the group ID
                    (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
                  }

                  if (it->GetGeometry()(0)->FastGetSolutionStepValue(GROUP_ID) == 2){
                    (it)->GetGeometry()(0)->Free(VELOCITY_Y);
                    rCurrentProcessInfo[FIXED_VEL_BOT] = (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y); //cutre way yeah!
                    //I only store one value for every ball in the group ID
                    (it)->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
                  }

              } // loop over particles

          } // loop threads OpenMP

          KRATOS_CATCH("")
      }
	  
	  
    void SetCoordinationNumber(ModelPart& r_model_part) 
    {
      double in_coordination_number = mCoordinationNumber;
      double out_coordination_number = ComputeCoordinationNumber();
      int iteration = 0;
 
      if(in_coordination_number <= 0.0)
      {
        KRATOS_ERROR(std::runtime_error,"The specified Coordination Number is less or equal to zero, N.C. = ",in_coordination_number)
      }
      else
      {
          while ( fabs(out_coordination_number/in_coordination_number-1.0) > 1e-3)
          {
    
            
            iteration++;

            mSearchTolerance *= in_coordination_number/out_coordination_number;
            
            SetSearchRadius(r_model_part, 1.0);
            
            SearchNeighbours();
             
            out_coordination_number = ComputeCoordinationNumber();
   
          }//while
          
          std::cout<< "Coordination Number iteration converged after "<<iteration<< " iterations, to value " <<out_coordination_number<<". "<<"\n"<<std::endl;
          
            
      }
      
    } //SetCoordinationNumber
    
    double ComputeCoordinationNumber()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        unsigned int size = 0;
        unsigned int total_contacts = 0;
        mNeighbourCounter[OpenMPUtils::ThisThread()] = 0.0;
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for //private(index, MassMatrix)  //M. proba de compilar sense mass matrix??
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){

                size = (it)->GetValue(NEIGHBOUR_ELEMENTS).size();
                mNeighbourCounter[OpenMPUtils::ThisThread()]+=size;
                
            }

        }
              
        for (int i = 0; i < this->GetNumberOfThreads(); i++)
        {
          
          total_contacts += mNeighbourCounter[OpenMPUtils::ThisThread()];
                  
        }
        
        return (double(total_contacts)/double(pElements.size()));
       
        KRATOS_CATCH("")
      
    }
    


    void CalculateEnergies(){}


    void InitializeElements()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for //private(index, MassMatrix)  //M. proba de compilar sense mass matrix??
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
            //  Element::GeometryType& geom = it->GetGeometry(); ///WARNING: COMMENTED AVOIDING WARNING COMPILATION
                (it)->Initialize();

            }

        }

        mElementsAreInitialized = true;

        KRATOS_CATCH("")
    }
    
    
    void SetSearchRadius(ModelPart& r_model_part, double amplification)
    {
        KRATOS_TRY
     
        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        
       // if(mNumberOfElementsOldRadiusList == number_of_elements) return;  //TODO:: this can fail when one particle is created and one is destroyed for example.
       // else 
        mNumberOfElementsOldRadiusList = number_of_elements;
        
        this->GetRadius().resize(number_of_elements);

        for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it){

            this->GetRadius()[particle_pointer_it - pElements.begin()] = amplification*(mSearchTolerance + particle_pointer_it->GetGeometry()(0)->GetSolutionStepValue(RADIUS));

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

    virtual void SearchNeighbours()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();

        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

        if(!number_of_elements) return;
        
        for (SpatialSearch::ElementsContainerType::iterator i = pElements.begin(); i != pElements.end(); i++){
            i->GetValue(NEIGHBOUR_ELEMENTS).clear();
        }        
        

        mpSpSearch->SearchElementsInRadiusExclusive(r_model_part, this->GetRadius(), this->GetResults(), this->GetResultsDistances());
//         mpSpSearch->SearchElementsInRadiusExclusive(r_model_part,this->GetRadius(),this->GetResults());
        
        for (SpatialSearch::ElementsContainerType::iterator i = r_model_part.GetCommunicator().GhostMesh().Elements().begin(); i != r_model_part.GetCommunicator().GhostMesh().Elements().end(); i++){
            r_model_part.Elements().push_back(*i);
        }
        
        SynchronizeSolidMesh(r_model_part);
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
        
        #pragma omp parallel for
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            size_t ResultCounter = this->GetElementPartition()[k];

            for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++ResultCounter){
                WeakPointerVector<Element>& neighbour_elements = particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS);
                for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[ResultCounter].begin(); neighbour_it != this->GetResults()[ResultCounter].end(); ++neighbour_it){
                    neighbour_elements.push_back(*neighbour_it);         
                }

                this->GetResults()[ResultCounter].clear();
                this->GetResultsDistances()[ResultCounter].clear();
                
                //SORTING NEIGHBOURS BY ID: (if you activate this, do not forget to activate it also in searchinitialneighbours
                //std::sort( (particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS)).ptr_begin(), (particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS)).ptr_end(), compare_ids );
                /*KRATOS_WATCH("new ball")
                for (WeakPointerVector<Element >::iterator nei_i = (particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS)).begin(); nei_i != (particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS)).end(); nei_i++){
                    KRATOS_WATCH(nei_i->Id())  
                }*/
            }

        }

        KRATOS_CATCH("")
    }        	
	
	 void ComputeNewNeighboursHistoricalData()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        ProcessInfo& rCurrentProcessInfo      = r_model_part.GetProcessInfo();
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for 
        
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            double dummy;
            
            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
           
                (it)->Calculate(CALCULATE_COMPUTE_NEW_NEIGHBOURS_HISTORICAL_DATA, dummy, rCurrentProcessInfo);
            
              
            }

        }

        KRATOS_CATCH("")
    }
    
     void ComputeNewRigidFaceNeighboursHistoricalData()
    {
        KRATOS_TRY

        ModelPart& r_model_part               = BaseType::GetModelPart();
        ElementsArrayType& pElements          = r_model_part.GetCommunicator().LocalMesh().Elements();
        ProcessInfo& rCurrentProcessInfo      = r_model_part.GetProcessInfo();
        
        OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());

        #pragma omp parallel for 
        
        for (int k = 0; k < this->GetNumberOfThreads(); k++){
            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            double dummy;
            
            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
           
                (it)->Calculate(CALCULATE_COMPUTE_NEW_RIGID_FACE_NEIGHBOURS_HISTORICAL_DATA, dummy, rCurrentProcessInfo);

            }

        }

        KRATOS_CATCH("")
    }

	
	////Cfeng
    virtual void SearchRigidFaceNeighbours()
    {
		
        KRATOS_TRY
		ElementsArrayType& pElements           = mpDem_model_part->GetCommunicator().LocalMesh().Elements();	
		ConditionsArrayType& pTContitions      = mpFem_model_part->GetCommunicator().LocalMesh().Conditions();     
		
		if(pTContitions.size() > 0)
		{
			OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
			
			////Cfeng: clear and swap the vector,initializaiton 
			#pragma omp parallel for
			for (int k = 0; k < this->GetNumberOfThreads(); k++)
			{
				typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
				typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
				
				for (SpatialSearch::ElementsContainerType::iterator E_pointer_it = it_begin; E_pointer_it != it_end; ++E_pointer_it)
				{
					WeakPointerVector<Condition > tempP;
					E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).swap(tempP);
					
					Vector tempV;
					E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_PRAM).swap(tempV);
					
					Vector tempV1;
					E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE).swap(tempV1);
				}
			}
		   

			moDemFemSearch.SearchRigidFaceForDEMInRadiusExclusiveImplementation(pElements, pTContitions, this->GetRadius(), this->GetRigidFaceResults(), this->GetRigidFaceResultsDistances());
			
			#pragma omp parallel for
			for (int k = 0; k < this->GetNumberOfThreads(); k++)
			{
				typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
				typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

				size_t ResultCounter = this->GetElementPartition()[k];
				
				for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++ResultCounter)
				{
                                    WeakPointerVector<Condition>& neighbour_rigid_faces = particle_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES);
					for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[ResultCounter].begin(); 
						 neighbour_it != this->GetRigidFaceResults()[ResultCounter].end(); ++neighbour_it)
					{
						neighbour_rigid_faces.push_back(*neighbour_it);   
					}

					this->GetRigidFaceResults()[ResultCounter].clear();
					this->GetRigidFaceResultsDistances()[ResultCounter].clear();
				
				}

			}
			
			
			
			
			//****************************************************************************
			/////************Cfeng: Below for DEM FEM coupling********************
			//*****************************************************************************
			
			SpatialSearch::ElementsContainerType::iterator E_pointer_it;
			
			//Cfeng:resize the particle-rigidface contact forces for each particles 
			#pragma omp parallel for
			for (int k = 0; k < this->GetNumberOfThreads(); k++)
			{
				typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
				typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
				
				for (SpatialSearch::ElementsContainerType::iterator E_pointer_it = it_begin; E_pointer_it != it_end; ++E_pointer_it)
				{				
					std::size_t totalno = E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).size() * 3;
					E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE).resize(totalno);
				}
			}
			
			typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
			
			///Cfeng:clear  neighbours for rigidfaces
			for (ConditionsArrayType::iterator ic = pTContitions.begin(); ic != pTContitions.end(); ic++)
			{
				WeakPointerVector<Element > tempP;
				ic->GetValue(NEIGHBOUR_PARTICLE_OF_RIGID_FACE).swap(tempP);
			}		
			
			////Cfeng: Find The particle neighbours for each RigidFace, used for calculating FEM force
			for (E_pointer_it = pElements.begin(); E_pointer_it != pElements.end(); ++E_pointer_it)
			{					
				for(ConditionWeakIteratorType ineighbour = E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).begin(); 
					 ineighbour != E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).end(); ineighbour++)
				{
					ineighbour->GetValue(NEIGHBOUR_PARTICLE_OF_RIGID_FACE).push_back(*(E_pointer_it.base()));				
				}
			}	
		}
	
				
        KRATOS_CATCH("")
    }
	
	
	
	
	////Not used now
	virtual void SearchRigidFaceNeighbours_Old_Method()
    {
        KRATOS_TRY

        ModelPart& r_model_part                = *mpDem_model_part;
		ElementsArrayType& pElements           = r_model_part.GetCommunicator().LocalMesh().Elements();
		ConditionsArrayType& pMContitions      = r_model_part.GetCommunicator().LocalMesh().Conditions();		
		ConditionsArrayType& pTContitions      = mpFem_model_part->GetCommunicator().LocalMesh().Conditions();   
		
		if(pTContitions.size() > 0)
		{
			for (SpatialSearch::ElementsContainerType::iterator E_pointer_it = pElements.begin(); E_pointer_it != pElements.end(); ++E_pointer_it)
			{
				WeakPointerVector<Condition > tempP;
				E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).swap(tempP);
				
				Vector tempV;
				E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_PRAM).swap(tempV);
				
				Vector tempV1;
				E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE).swap(tempV1);
			}
		   
		   
			for (ConditionsArrayType::iterator i = pMContitions.begin(); i != pMContitions.end(); i++)
			{
				WeakPointerVector<Condition > tempP;
				i->GetValue(NEIGHBOUR_RIGID_FACES).swap(tempP);
				
				Vector tempV;
				i->GetValue(NEIGHBOUR_RIGID_FACES_PRAM).swap(tempV);
			}
			
			


			moDemFemSearch.SearchConditionsInRadiusExclusive(r_model_part, pTContitions, this->GetRadius(), this->GetRigidFaceResults(), this->GetRigidFaceResultsDistances());
			
			
			OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pMContitions.size(), this->GetConditionPartition());
			
			
			#pragma omp parallel for
			for (int k = 0; k < this->GetNumberOfThreads(); k++)
			{
				typename ConditionsArrayType::iterator it_begin = pMContitions.ptr_begin() + this->GetConditionPartition()[k];
				typename ConditionsArrayType::iterator it_end   = pMContitions.ptr_begin() + this->GetConditionPartition()[k + 1];

				size_t ResultCounter = this->GetConditionPartition()[k];
				
				for (SpatialSearch::ConditionsContainerType::iterator particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it,++ResultCounter)
				{
                                    WeakPointerVector<Condition>& neighbour_rigid_faces = particle_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES);
					for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[ResultCounter].begin(); 
						 neighbour_it != this->GetRigidFaceResults()[ResultCounter].end(); ++neighbour_it)
					{
						neighbour_rigid_faces.push_back(*neighbour_it);         
					}

					this->GetRigidFaceResults()[ResultCounter].clear();
					this->GetRigidFaceResultsDistances()[ResultCounter].clear();
				}

			}
			
			
			typename ConditionsArrayType::iterator C_pointer_it = pMContitions.ptr_begin();
			
			SpatialSearch::ElementsContainerType::iterator E_pointer_it;
			
			/////#pragma omp parallel for 
			for (E_pointer_it = pElements.begin(); E_pointer_it != pElements.end(); ++E_pointer_it, ++C_pointer_it)
			{
				E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).swap(C_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES));
				E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_PRAM).swap(C_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_PRAM));
				
				std::size_t totalno = E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).size() * 3;
				E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE).resize(totalno);
			}
			
			
			
			
			for (ConditionsArrayType::iterator ic = pTContitions.begin(); ic != pTContitions.end(); ic++)
			{
				WeakPointerVector<Element > tempP;
				ic->GetValue(NEIGHBOUR_PARTICLE_OF_RIGID_FACE).swap(tempP);
			}
			///push_back the particle point to rigid face conditon
			typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
			
			for (E_pointer_it = pElements.begin(); E_pointer_it != pElements.end(); ++E_pointer_it)
			{
                                WeakPointerVector<Condition>& neighbour_rigid_faces = E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES);
				for(ConditionWeakIteratorType ineighbour = neighbour_rigid_faces.begin(); 
					 ineighbour != neighbour_rigid_faces.end(); ineighbour++)
				{
					ineighbour->GetValue(NEIGHBOUR_PARTICLE_OF_RIGID_FACE).push_back(*(E_pointer_it.base()));				
				}
			}
			
		}

        KRATOS_CATCH("")
    }
	 
	
	
	void Compute_RigidFace_Movement()
	{
		KRATOS_TRY
		
        ProcessInfo& rCurrentProcessInfo  = mpDem_model_part->GetProcessInfo();
		
		if(rCurrentProcessInfo[RIGID_FACE_FLAG])
		{
			double delta_t                    = rCurrentProcessInfo[DELTA_TIME];
			int time_step                     = rCurrentProcessInfo[TIME_STEPS];
			
			double begin_time  = rCurrentProcessInfo[RIGID_FACE_BEGIN_TIME];
			double end_time    = rCurrentProcessInfo[RIGID_FACE_END_TIME];
			double time_now    = delta_t * time_step;
			
			if(time_now >= begin_time && time_now <= end_time)
			{
				
				int PropID         = rCurrentProcessInfo[RIGID_FACE_PROP_ID];
				
				ConditionsArrayType& pContitions = mpFem_model_part->GetCommunicator().LocalMesh().Conditions();
			
				OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pContitions.size(), this->GetConditionPartition());
				
				/////Calculate node velocity
				#pragma omp parallel for
				for (int k = 0; k < this->GetNumberOfThreads(); k++)
				{
					typename ConditionsArrayType::iterator it_begin = pContitions.ptr_begin() + this->GetConditionPartition()[k];
					typename ConditionsArrayType::iterator it_end   = pContitions.ptr_begin() + this->GetConditionPartition()[k + 1];
					
					for (ConditionsArrayType::iterator ip = it_begin; ip != it_end; ++ip)
					{
						if(static_cast<int>(ip->GetProperties().Id()) == PropID)
						{
							Vector VelocityArray;
							ip->Calculate(RIGID_FACE_COMPUTE_MOVEMENT, VelocityArray, rCurrentProcessInfo);
							
							Condition::GeometryType& geom = ip->GetGeometry();
							
							const unsigned int& dim       = geom.WorkingSpaceDimension();
							
							for (unsigned int i = 0; i <geom.size(); i++)
							{
								unsigned int index = i * dim;
								
								array_1d<double,3>& node_vel = geom(i)->FastGetSolutionStepValue(VELOCITY);
								
								for(unsigned int kk=0; kk < dim; kk++)
								{
									geom(i)->SetLock();
									node_vel[kk] = VelocityArray[index+kk];
									geom(i)->UnSetLock();
								}
							}
							
						}
					}
				}
				
				
				/////Feng Chun::Calculate Node Movement
				
				vector<unsigned int> node_partition;
						
				NodesArrayType& pNodes = mpFem_model_part->GetCommunicator().LocalMesh().Nodes();
				OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pNodes.size(), node_partition);
				
				#pragma omp parallel for
				for(int k = 0; k < this->GetNumberOfThreads(); k++)
				{
					  NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
					  NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
					 
					  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
					  {      
						  array_1d<double, 3 > & vel             = i->FastGetSolutionStepValue(VELOCITY);
						  array_1d<double, 3 > & displ           = i->FastGetSolutionStepValue(DISPLACEMENT);
						  array_1d<double, 3 > & coor            = i->Coordinates();
						  array_1d<double, 3 > & initial_coor    = i->GetInitialPosition();
						  						
						  displ += vel * delta_t;
						  coor   = initial_coor + displ;
					  }
				}
				
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
//         r_model_part.GetCommunicator().SynchronizeDofs();
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

