//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//

#if !defined(KRATOS_CONTINUUM_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_CONTINUUM_EXPLICIT_SOLVER_STRATEGY

#include "custom_strategies/strategies/explicit_solver_strategy.h"

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::


namespace Kratos
{
  
  template<
  class TSparseSpace,
  class TDenseSpace,
  class TLinearSolver>
  class ContinuumExplicitSolverStrategy: public ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
  {
      public:
      ///@name Type Definitions
      ///@{
        
      typedef ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>   BaseType;
   
      typedef typename BaseType::NodesArrayType                             NodesArrayType;
      typedef typename BaseType::ElementsArrayType                          ElementsArrayType;
      typedef typename BaseType::ElementsIterator                           ElementsIterator;
      typedef typename BaseType::ConditionsArrayType                        ConditionsArrayType;
     
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ContinuumExplicitSolverStrategy);

      /// Default constructor.
      ContinuumExplicitSolverStrategy(){}

      ContinuumExplicitSolverStrategy(
                             ModelPart& model_part,
                             ModelPart& contacts_model_part, 
                             const double enlargement_factor,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool MoveMeshFlag,
                             const bool delta_option,
                             const bool continuum_simulating_option,
                             typename IntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch
      ): ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, enlargement_factor, max_delta_time, n_step_search,safety_factor,MoveMeshFlag,pScheme,pSpSearch), mcontacts_model_part(contacts_model_part)
      {
          mdelta_option                 = delta_option;
          mcontinuum_simulating_option  = continuum_simulating_option;
          
      }

      /// Destructor.
      virtual ~ContinuumExplicitSolverStrategy()
      {
         Timer::SetOuputFile("TimesPartialRelease");
         Timer::PrintTimingInformation();
      }

      virtual void Initialize()
      {
        
         KRATOS_TIMER_START("INITIALIZE")
          KRATOS_WATCH("---------------------CONTINUUM EXPLICIT SOLVER STRATEGY-------------------------------")
               
          KRATOS_TRY

          ModelPart& rModelPart            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
                
          int NumberOfElements = rModelPart.GetCommunicator().LocalMesh().ElementsArray().end() - rModelPart.GetCommunicator().LocalMesh().ElementsArray().begin();
          
          this->GetResults().resize(NumberOfElements);
          this->GetResultsDistances().resize(NumberOfElements);
          this->GetRadius().resize(NumberOfElements);
          
          // Omp initializations
          this->GetNumberOfThreads() = OpenMPUtils::GetNumThreads();
     
          // 0. Set search radius
          BaseType::SetSearchRadius(rModelPart,rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION],0);
          
          // 1. Search Neighbours with tolerance (Not in mpi.)
          bool extension_option            = false;
          this->GetBoundingBoxOption()     = rCurrentProcessInfo[BOUNDING_BOX_OPTION];
          extension_option                 = rCurrentProcessInfo[CASE_OPTION];

          // 2. Initializing elements and perform the repartition
          if (this->GetElementsAreInitialized() == false){
              BaseType::InitializeElements();
          }

          
          this->GetInitializeWasPerformed() = true;

          // 3. Search Neighbours with tolerance (after first repartition process)
          this->SearchNeighbours(rModelPart, extension_option);
          
          //the search radius is modified for the next steps.
        
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();
          for (SpatialSearch::ElementsContainerType::iterator particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
          {
           
              this->GetRadius()[particle_pointer_it - pElements.begin()] *= (rCurrentProcessInfo[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION]); 
             
          }
          
          
          // 4. Set Initial Contacts
           BaseType::InitializeSolutionStep();
          
           this->SetInitialContacts(); //MSIMSI 8 Yet to define

          // 5. Calculate bounding box
          this->GetParticleCreatorDestructor().CalculateSurroundingBoundingBox(rModelPart, this->mEnlargementFactor);

          KRATOS_CATCH("")
          
           KRATOS_TIMER_STOP("INITIALIZE")

      }// Initialize()

      virtual double Solve()
      {
            
          KRATOS_TRY
          
             KRATOS_TIMER_START("SOLVE")
          
          
          KRATOS_TIMER_START("BEGIN")
          ModelPart& rModelPart            = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

          bool extension_option            = false;

          if (rCurrentProcessInfo[CASE_OPTION]){
              extension_option             = true;
          }

          int time_step = rCurrentProcessInfo[TIME_STEPS];
          KRATOS_TIMER_STOP("BEGIN")

          //STRATEGY:
          // 1. Here we initialize member variables that depend on the rCurrentProcessInfo
          KRATOS_TIMER_START("InitializeSolutionStep")
          BaseType::InitializeSolutionStep();
          KRATOS_TIMER_STOP("InitializeSolutionStep")
          
          // 2. Get and Calculate the forces
          KRATOS_TIMER_START("GetForce")
          BaseType::GetForce();
          KRATOS_TIMER_STOP("GetForce")
           
          // 3. Motion Integration
          KRATOS_TIMER_START("PerformTimeIntegrationOfMotion")
          BaseType::PerformTimeIntegrationOfMotion(); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot
          KRATOS_TIMER_STOP("PerformTimeIntegrationOfMotion")
          
          // 4. Synchronize  
          KRATOS_TIMER_START("SynchronizeSolidMesh")
          BaseType::SynchronizeSolidMesh(rModelPart);
          KRATOS_TIMER_STOP("SynchronizeSolidMesh")
          

          // 5. Neighbouring search. Every N times. + destruction of particles outside the bounding box
          KRATOS_TIMER_START("SearchNeighbours")
          if (rCurrentProcessInfo[ACTIVATE_SEARCH] == 1){

              if ((time_step + 1)%this->GetNStepSearch() == 0 && time_step > 0){

                  if (this->GetBoundingBoxOption() == 1){
                      BaseType::BoundingBoxUtility();
                  }

                   this->SearchNeighbours(rModelPart, extension_option); //the amplification factor has benn modified after the first search. 
              }

          }
          KRATOS_TIMER_STOP("SearchNeighbours")
          
          KRATOS_TIMER_START("FinalizeSolutionStep")
          BaseType::FinalizeSolutionStep();
          KRATOS_TIMER_STOP("FinalizeSolutionStep")
          
          
      KRATOS_TIMER_STOP("SOLVE")
          
          
          return 0.00;
                
          
          KRATOS_CATCH("") 
          
          
      }//Solve()
      
      
      void SetInitialContacts(){}

   
    
    protected:
    
    ModelPart& mcontacts_model_part;
    bool   mdelta_option;
    bool   mcontinuum_simulating_option;
    

  }; // Class ContinuumExplicitSolverStrategy


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined




