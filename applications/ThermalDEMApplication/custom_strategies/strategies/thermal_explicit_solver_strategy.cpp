//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "thermal_explicit_solver_strategy.h"

namespace Kratos {

  //=====================================================================================================================================================================================
  // Constructor/Destructor methods

  ThermalExplicitSolverStrategy::ThermalExplicitSolverStrategy():ExplicitSolverStrategy() {}

  ThermalExplicitSolverStrategy::ThermalExplicitSolverStrategy(ExplicitSolverSettings&            settings,
                                                               const double                       max_delta_time,
                                                               const int                          n_step_search,
                                                               const double                       safety_factor,
                                                               const int                          delta_option,
                                                               ParticleCreatorDestructor::Pointer p_creator_destructor,
                                                               DEM_FEM_Search::Pointer            p_dem_fem_search,
                                                               SpatialSearch::Pointer             pSpSearch,
                                                               Parameters                         strategy_parameters):
                                                               ExplicitSolverStrategy(settings,
                                                                                      max_delta_time,
                                                                                      n_step_search,
                                                                                      safety_factor,
                                                                                      delta_option,
                                                                                      p_creator_destructor,
                                                                                      p_dem_fem_search,
                                                                                      pSpSearch,
                                                                                      strategy_parameters) {}

  ThermalExplicitSolverStrategy::~ThermalExplicitSolverStrategy() {}

  //=====================================================================================================================================================================================
  // Derived methods

  //------------------------------------------------------------------------------------------------------------
  void ThermalExplicitSolverStrategy::SetSearchRadiiOnAllParticles(ModelPart& r_model_part, double added_search_distance, double amplification) {
    SetSearchRadii(r_model_part, added_search_distance, amplification);
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalExplicitSolverStrategy::SetSearchRadiiWithFemOnAllParticles(ModelPart& r_model_part, double added_search_distance, double amplification) {
    SetSearchRadii(r_model_part, added_search_distance, amplification);
  }

  //=====================================================================================================================================================================================
  // Particular methods

  //------------------------------------------------------------------------------------------------------------
  // Solve solution step ignoring particles kinetics (forces and motion). Should be called when computing only heat transfer.
  double ThermalExplicitSolverStrategy::SolveSolutionStepStatic(void) {
    KRATOS_TRY

    ExplicitSolverStrategy::GetForce();
    PerformThermalTimeIntegration();
    return 0.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalExplicitSolverStrategy::PerformThermalTimeIntegration(void) {
    KRATOS_TRY
    
    ProcessInfo& r_process_info   = GetModelPart().GetProcessInfo();
    const int number_of_particles = (int)mListOfSphericParticles.size();

    #pragma omp parallel for
    for (int i = 0; i < number_of_particles; i++) {
      ThermalSphericParticle* particle = dynamic_cast<ThermalSphericParticle*>(mListOfSphericParticles[i]);
      particle->Move(r_process_info[DELTA_TIME], false, 0.0, 0);
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  // Set search distance for particle-particle and particle-wall by adding the additional distance due to heat transfer models.
  void ThermalExplicitSolverStrategy::SetSearchRadii(ModelPart& r_model_part, double added_search_distance, double amplification) {
    KRATOS_TRY
    
    int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

    IndexPartition<unsigned int>(number_of_elements).for_each([&](unsigned int i) {
      ThermalSphericParticle* particle = dynamic_cast<ThermalSphericParticle*>(mListOfSphericParticles[i]);

      particle->ComputeAddedSearchDistance(r_model_part.GetProcessInfo(), added_search_distance);
      particle->SetSearchRadius(amplification * (added_search_distance + mListOfSphericParticles[i]->GetRadius()));
    });

    KRATOS_CATCH("")
  }

} // namespace Kratos
