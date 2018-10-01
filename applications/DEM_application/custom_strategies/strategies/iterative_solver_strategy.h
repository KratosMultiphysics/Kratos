//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_ITERATIVE_SOLVER_STRATEGY)
#define  KRATOS_ITERATIVE_SOLVER_STRATEGY
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

namespace Kratos
{
  class IterativeSolverStrategy: public ExplicitSolverStrategy
  {
      public:

      typedef ExplicitSolverStrategy  BaseType;

      typedef BaseType::NodesArrayType                             NodesArrayType;
      typedef BaseType::ElementsArrayType                          ElementsArrayType;
      typedef BaseType::ElementsIterator                           ElementsIterator;
      typedef BaseType::ConditionsArrayType                        ConditionsArrayType;

      typedef WeakPointerVector<Element> ParticleWeakVectorType;
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(IterativeSolverStrategy);

      /// Default constructor.
      IterativeSolverStrategy(){}

      IterativeSolverStrategy(
                             ExplicitSolverSettings& settings,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const int delta_option,
                             ParticleCreatorDestructor::Pointer p_creator_destructor,
                             DEM_FEM_Search::Pointer p_dem_fem_search,
                             SpatialSearch::Pointer pSpSearch,
                             const bool do_search_balls = true)
      :ExplicitSolverStrategy(settings, max_delta_time, n_step_search, safety_factor, delta_option, p_creator_destructor, p_dem_fem_search, pSpSearch, do_search_balls)
      {
          BaseType::GetParticleCreatorDestructor()   = p_creator_destructor;
      }

      /// Destructor.
      virtual ~IterativeSolverStrategy()
      {
         //Timer::SetOuputFile("TimesPartialRelease");
         //Timer::PrintTimingInformation();
      }


      virtual void Initialize() override
      {

        KRATOS_TRY
        ModelPart& r_model_part = BaseType::GetModelPart();
        BaseType::Initialize();
        BaseType::InitializeSolutionStep();
        BaseType::ForceOperations(r_model_part);
        BaseType::FinalizeSolutionStep();
        //SchemeInitialize();

        KRATOS_CATCH("")
      }// Initialize()


       void SchemeInitialize()
      {
        PerformTimeIntegrationOfMotion(0);
          //BaseType::GetScheme()->Calculate(BaseType::GetModelPart(),0);
          //BaseType::GetScheme()->Calculate(BaseType::GetClusterModelPart(),0);
      }


      void SchemePredict()
      {
          PerformTimeIntegrationOfMotion(1);
          //BaseType::GetScheme()->Calculate(BaseType::GetModelPart(),1); //TODO: better call Predict function (would be empty in general)
          //BaseType::GetScheme()->Calculate(BaseType::GetClusterModelPart(),1);
      }

      void SchemeCorrect()
      {
          PerformTimeIntegrationOfMotion(2);
        //BaseType::GetScheme()->Calculate(BaseType::GetModelPart(),2); //TODO: better call Correct function (normal operations would be in that method)
        //BaseType::GetScheme()->Calculate(BaseType::GetClusterModelPart(),2);
      }

      virtual double Solve() override
      {

        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();

        BaseType::InitializeSolutionStep();
        SchemePredict();
        BaseType::SearchDEMOperations(r_model_part);
        BaseType::SearchFEMOperations(r_model_part);
        BaseType::ForceOperations(r_model_part);
        SchemeCorrect();
        BaseType::FinalizeSolutionStep();

        return 0.00;

        KRATOS_CATCH("")

      }//Solve()

  };//ClassIterativeSolverStrategy




}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
