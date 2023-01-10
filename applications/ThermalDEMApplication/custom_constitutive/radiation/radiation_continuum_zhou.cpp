//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "radiation_continuum_zhou.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RadiationContinuumZhou::RadiationContinuumZhou() {}
  RadiationContinuumZhou::~RadiationContinuumZhou() {}

  //------------------------------------------------------------------------------------------------------------
  double RadiationContinuumZhou::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check if particles are close enough
    // TODO: radiation with walls not yet implemented
    if ((!particle->CheckSurfaceDistance(r_process_info[MAX_RADIATION_DISTANCE])) || (particle->mNeighborType & WALL_NEIGHBOR))
      return 0.0;

    // In continuous models, environment temperature is computed by accumulating each neighbor contribution
    // (heat flux is computed after the loop over neighbors)
    AccumulateEnvironmentTemperature(r_process_info, particle);
    return 0.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double RadiationContinuumZhou::FinalizeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check if radiation neighbors exist
    const int neighbors = particle->mRadiativeNeighbors;
    if (neighbors == 0)
      return 0.0;

    const double accum_temperature    = particle->mEnvironmentTemperature;
    const double particle_emissivity  = particle->GetParticleEmissivity();
    const double particle_surface     = particle->GetParticleSurfaceArea();
    const double particle_temperature = particle->GetParticleTemperature();
    const double porosity             = r_process_info[AVERAGE_POROSITY];
    const double fluid_temperature    = r_process_info[FLUID_TEMPERATURE];    

    // Compute environment temperature
    const double env_temperature = porosity * fluid_temperature + (1.0 - porosity) * accum_temperature / neighbors;

    // Compute heat flux
    return STEFAN_BOLTZMANN * particle_emissivity * particle_surface * (pow(env_temperature,4.0) - pow(particle_temperature,4.0));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void RadiationContinuumZhou::AccumulateEnvironmentTemperature(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Update number of radiative neighbors
    particle->mRadiativeNeighbors++;

    // Accumulate temperature
    particle->mEnvironmentTemperature += particle->GetNeighborTemperature();

    KRATOS_CATCH("")
  }

} // namespace Kratos
