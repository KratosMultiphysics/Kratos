#if !defined(KRATOS_RUNGE_KUTTA_SOLVER_STRATEGY)
#define  KRATOS_RUNGE_KUTTA_SOLVER_STRATEGY
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#define CUSTOMTIMER 0

namespace Kratos
{
  class RK4SolverStrategy: public ExplicitSolverStrategy
  {
      public:

      typedef ExplicitSolverStrategy  BaseType;

      typedef BaseType::NodesArrayType                             NodesArrayType;
      typedef BaseType::ElementsArrayType                          ElementsArrayType;
      typedef BaseType::ElementsIterator                           ElementsIterator;
      typedef BaseType::ConditionsArrayType                        ConditionsArrayType;

      typedef GlobalPointersVector<Element> ParticleWeakVectorType;
      typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(RK4SolverStrategy);

      /// Default constructor.
      RK4SolverStrategy(){}

      RK4SolverStrategy(ExplicitSolverSettings& settings,
                        const double max_delta_time,
                        const double n_step_search,
                        const double safety_factor,
                        const int delta_option,
                        ParticleCreatorDestructor::Pointer p_creator_destructor,
                        DEM_FEM_Search::Pointer p_dem_fem_search,
                        SpatialSearch::Pointer pSpSearch,
                        Parameters strategy_parameters)
      :ExplicitSolverStrategy(settings, max_delta_time, n_step_search, safety_factor, delta_option, p_creator_destructor, p_dem_fem_search, pSpSearch, strategy_parameters)
      {
          BaseType::GetParticleCreatorDestructor()   = p_creator_destructor;
      }

      /// Destructor.
      virtual ~RK4SolverStrategy(){}

      virtual void Initialize() override {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        BaseType::Initialize();
        BaseType::InitializeSolutionStep();
        BaseType::ForceOperations(r_model_part);
        BaseType::FinalizeSolutionStep();

        KRATOS_CATCH("")
      }

      void SchemeRKStepInit(){
        KRATOS_WATCH("SchemeRKStepInit")
        PerformTimeIntegrationOfMotion(0);
      }

      void SchemeRKStep1(){
        KRATOS_WATCH("SchemeRKStep1")
          PerformTimeIntegrationOfMotion(1);
      }

      void SchemeRKStep2(){
        KRATOS_WATCH("SchemeRKStep2")
          PerformTimeIntegrationOfMotion(2);
      }

      void SchemeRKStep3(){
        KRATOS_WATCH("SchemeRKStep3")
          PerformTimeIntegrationOfMotion(3);
      }

      virtual double Solve() override {

        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();

        BaseType::InitializeSolutionStep();

        BaseType::SearchDEMOperations(r_model_part);
        BaseType::SearchFEMOperations(r_model_part);
        BaseType::ForceOperations(r_model_part);

        SchemeRKStepInit();
        BaseType::SearchDEMOperations(r_model_part);
        BaseType::SearchFEMOperations(r_model_part);
        BaseType::ForceOperations(r_model_part);

        SchemeRKStep1();
        BaseType::SearchDEMOperations(r_model_part);
        BaseType::SearchFEMOperations(r_model_part);
        BaseType::ForceOperations(r_model_part);

        SchemeRKStep2();
        BaseType::SearchDEMOperations(r_model_part);
        BaseType::SearchFEMOperations(r_model_part);
        BaseType::ForceOperations(r_model_part);

        SchemeRKStep3();

        BaseType::FinalizeSolutionStep();

        return 0.00;

        KRATOS_CATCH("")
      }
  };
}

#endif
