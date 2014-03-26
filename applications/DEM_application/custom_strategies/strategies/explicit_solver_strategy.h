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
          OpenMPUtils::CreatePartition(this->GetNumberOfThreads(), pElements.size(), this->GetElementPartition());
          
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
              
              aux_pointer = &( props_it->GetValue(ROLLING_FRICTION) );
              aux_props.SetRollingFrictionFromProperties(aux_pointer);
              
              aux_pointer = &( props_it->GetValue(PARTICLE_FRICTION) );
              aux_props.SetTgOfFrictionAngleFromProperties(aux_pointer);
              
              aux_pointer = &( props_it->GetValue(LN_OF_RESTITUTION_COEFF) );
              aux_props.SetLnOfRestitCoeffFromProperties(aux_pointer);
              
              //ROLLING FRICTION???                            
              mFastProperties[i] = aux_props;
              i++;
              
      }

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

          //CreatePropertiesProxies();
          
          int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

          this->GetResults().resize(number_of_elements);
          this->GetResultsDistances().resize(number_of_elements);          
          
          //Cfeng
          this->GetRigidFaceResults().resize(number_of_elements);
      this->GetRigidFaceResultsDistances().resize(number_of_elements);

          // Omp initializations
          this->GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
          mNeighbourCounter.resize(this->GetNumberOfThreads());

          this->GetBoundingBoxOption() = rCurrentProcessInfo[BOUNDING_BOX_OPTION];

          InitializeSolutionStep();
          InitializeElements();

          mInitializeWasPerformed = true;
          
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
          FinalizeSolutionStep();
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

            //WeakPointerVector<Element>& neighbour_elements = particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS);
            //neighbour_elements.clear();
                
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
                //WeakPointerVector<Element>& neighbour_elements = particle_pointer_it->GetValue(NEIGHBOUR_ELEMENTS);
                //neighbour_elements.clear();                
                SphericParticle* spheric_central_particle = dynamic_cast<Kratos::SphericParticle*>( &(*particle_pointer_it) );
                //spheric_central_particle->mNeighbourElements.clear();
                
                for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[ResultCounter].begin(); neighbour_it != this->GetResults()[ResultCounter].end(); ++neighbour_it){
                    //neighbour_elements.push_back( *neighbour_it );   
                    Element* p_neighbour_element = (*neighbour_it).get();
                    SphericParticle* p_spheric_neighbour_particle = static_cast<SphericParticle*>( p_neighbour_element );
                    spheric_central_particle->mNeighbourElements.push_back( p_spheric_neighbour_particle );

                    if( (*neighbour_it)->Id() != (*p_spheric_neighbour_particle).Id() ) KRATOS_ERROR(std::logic_error,"ELEEEEEEEEE",(*neighbour_it)->Id() );                                                
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
        //ProcessInfo& rCurrentProcessInfo      = mpDem_model_part->GetProcessInfo();
        
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
                    E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_TOTAL_CONTACT_FORCE).swap(tempV1);
                    
                    Vector tempV2;
                    E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_ELASTIC_CONTACT_FORCE).swap(tempV2);
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
            
            //SpatialSearch::ElementsContainerType::iterator E_pointer_it;
            
            //Cfeng:resize the particle-rigidface contact forces for each particles 
            
            #pragma omp for
            for (int k = 0; k < this->GetNumberOfThreads(); k++)
            {
                typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
                typename ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];
                
                for (SpatialSearch::ElementsContainerType::iterator E_pointer_it = it_begin; E_pointer_it != it_end; ++E_pointer_it)
                {               
                    std::size_t totalno = E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES).size() * 3;
                    
                    E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_TOTAL_CONTACT_FORCE).resize(totalno);
                    E_pointer_it->GetValue(NEIGHBOUR_RIGID_FACES_ELASTIC_CONTACT_FORCE).resize(totalno);
                    
                }
            }
            } //end of the parallel region
            
            typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
            
            //MSIMSI 3: parallelize these two loops
            
            ///Cfeng:clear  neighbours for rigidfaces
            for (ConditionsArrayType::iterator ic = pTContitions.begin(); ic != pTContitions.end(); ic++)
            {
                WeakPointerVector<Element > tempP;
                ic->GetValue(NEIGHBOUR_PARTICLE_OF_RIGID_FACE).swap(tempP);
            }       
            
            ////Cfeng: Find The particle neighbours for each RigidFace, used for calculating FEM force
            //This loop can not be parallelized easily, two threads could be writing on the same condition!!
            for (SpatialSearch::ElementsContainerType::iterator E_pointer_it = pElements.begin(); E_pointer_it != pElements.end(); ++E_pointer_it)
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

                //if (!(it->GetGeometry()(0)->pGetDof(VELOCITY_X)->IsFixed())){
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

