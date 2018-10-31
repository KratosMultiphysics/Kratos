//
// Authors: Ferran Arrufat farrufat@cimne.upc.edu,
//          Guillermo Casas gcasas@cimne.upc.edu
//

#if !defined(KRATOS_VERLET_SOLVER_STRATEGY)
#define  KRATOS_VERLET_SOLVER_STRATEGY
#include "custom_strategies/strategies/explicit_solver_continuum.h"
#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

namespace Kratos
{
  template<class TBaseStrategy>
  class VelocityVerletSolverStrategy: public TBaseStrategy
  {
      public:

      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(VelocityVerletSolverStrategy);

      /// Default constructor.
      VelocityVerletSolverStrategy(){}

      VelocityVerletSolverStrategy(
                             ExplicitSolverSettings& settings,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const int delta_option,
                             ParticleCreatorDestructor::Pointer p_creator_destructor,
                             DEM_FEM_Search::Pointer p_dem_fem_search,
                             SpatialSearch::Pointer pSpSearch,
                             Parameters strategy_parameters)

      :TBaseStrategy(settings, max_delta_time, n_step_search, safety_factor, delta_option, p_creator_destructor, p_dem_fem_search, pSpSearch, strategy_parameters)
      {}

      /// Destructor.
      virtual ~VelocityVerletSolverStrategy()
      {
         Timer::SetOuputFile("TimesPartialRelease");
         Timer::PrintTimingInformation();
      }

      void Initialize() override
      {

        KRATOS_TRY
        ModelPart& r_model_part = this->GetModelPart();
        TBaseStrategy::Initialize();
        this->InitializeSolutionStep();
        this->ForceOperations(r_model_part);
        this->FinalizeSolutionStep();

        KRATOS_CATCH("")
      }

      double Solve() override
      {
          KRATOS_TRY
          ModelPart& r_model_part = this->GetModelPart();

          bool has_mpi = false;
          VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
          if(r_modelpart_nodal_variables_list.Has(PARTITION_INDEX) )  has_mpi = true;

          this->InitializeSolutionStep();

          this->PerformTimeIntegrationOfMotion(1);

          this->SearchDEMOperations(r_model_part, has_mpi);
          this->SearchFEMOperations(r_model_part, has_mpi);
          this->ForceOperations(r_model_part);

          this->PerformTimeIntegrationOfMotion(2);

          this->FinalizeSolutionStep();

          KRATOS_CATCH("")

          return 0.0;

      }//Solve()

  }; // Class VelocityVerletSolverStrategy<TBaseStrategy>


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
