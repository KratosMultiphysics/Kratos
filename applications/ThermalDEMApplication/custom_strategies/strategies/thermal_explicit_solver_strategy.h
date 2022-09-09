//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(KRATOS_THERMAL_EXPLICIT_SOLVER_STRATEGY_H_INCLUDED)
#define KRATOS_THERMAL_EXPLICIT_SOLVER_STRATEGY_H_INCLUDED

// System includes

// External includes
#include "custom_strategies/strategies/explicit_solver_strategy.h"

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) ThermalExplicitSolverStrategy : public ExplicitSolverStrategy
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(ThermalExplicitSolverStrategy);

      using ExplicitSolverStrategy::mListOfSphericParticles;

      // Constructor
      ThermalExplicitSolverStrategy();

      ThermalExplicitSolverStrategy(ExplicitSolverSettings&            settings,
                                    const double                       max_delta_time,
                                    const int                          n_step_search,
                                    const double                       safety_factor,
                                    const int                          delta_option,
                                    ParticleCreatorDestructor::Pointer p_creator_destructor,
                                    DEM_FEM_Search::Pointer            p_dem_fem_search,
                                    SpatialSearch::Pointer             pSpSearch,
                                    Parameters                         strategy_parameters);

      // Destructor
      virtual ~ThermalExplicitSolverStrategy();

      // Public derived methods
      void Initialize                          (void) override;
      void SetSearchRadiiOnAllParticles        (ModelPart& r_model_part, double added_search_distance = 0.0, double amplification = 1.0) override;
      void SetSearchRadiiWithFemOnAllParticles (ModelPart& r_model_part, double added_search_distance = 0.0, double amplification = 1.0) override;

      // Public particular methods
      double SolveSolutionStepStatic(void);

    protected:

      // Protected particular methods
      void SetSolveFrequency             (void);
      void PerformThermalTimeIntegration (void);
      void SetSearchRadii                (ModelPart & r_model_part, double added_search_distance, double amplification);

  }; // Class ThermalExplicitSolverStrategy
} // namespace Kratos

#endif // KRATOS_THERMAL_EXPLICIT_SOLVER_STRATEGY_H_INCLUDED defined
