
#if !defined(KRATOS_VERLET_SOLVER_STRATEGY)
#define  KRATOS_VERLET_SOLVER_STRATEGY

#include "custom_strategies/strategies/explicit_solver_continuum.h"
#include "DEM_definitions.h"

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::


namespace Kratos
{

  template<
  class TSparseSpace,
  class TDenseSpace,
  class TLinearSolver>
  class VerletVelocitySolverStrategy: public ContinuumExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
  {
      public:
          
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(VerletVelocitySolverStrategy);

      /// Default constructor.
      VerletVelocitySolverStrategy(){}
     
      VerletVelocitySolverStrategy(
                             ExplicitSolverSettings& settings,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const int delta_option,
                             const double search_tolerance,
                             const double coordination_number,
                             typename ParticleCreatorDestructor::Pointer p_creator_destructor,
                             typename DEM_FEM_Search::Pointer p_dem_fem_search,
                             typename DEMIntegrationScheme::Pointer pScheme,
                             typename SpatialSearch::Pointer pSpSearch)
      :ContinuumExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(settings, max_delta_time, n_step_search, safety_factor, delta_option, search_tolerance, coordination_number, p_creator_destructor, p_dem_fem_search, pScheme, pSpSearch)
      {}

      /// Destructor.
      virtual ~VerletVelocitySolverStrategy()
      {
         Timer::SetOuputFile("TimesPartialRelease");
         Timer::PrintTimingInformation();
      }      
      

      virtual double Solve()
      {
          KRATOS_TRY
          ModelPart& r_model_part            = this->GetModelPart();

          bool has_mpi = false;
          VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
          if(r_modelpart_nodal_variables_list.Has(PARTITION_INDEX) )  has_mpi = true;

          this->InitializeSolutionStep();
          
          this->PerformTimeIntegrationOfMotion(1);
          
          this->SearchOperations(r_model_part, has_mpi);
          this->ForceOperations(r_model_part);
          
          this->PerformTimeIntegrationOfMotion(2);
          
          this->FinalizeSolutionStep();

          KRATOS_CATCH("")

          return 0.0;

      }//Solve()

    protected:
    
    bool   mcontinuum_simulating_option;
    int    mFixSwitch;
    //bool   mDempackOption;
    std::vector<SphericContinuumParticle*>  mListOfSphericContinuumParticles;
    std::vector<SphericContinuumParticle*>  mListOfGhostSphericContinuumParticles;

  }; // Class VerletVelocitySolverStrategy


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined




