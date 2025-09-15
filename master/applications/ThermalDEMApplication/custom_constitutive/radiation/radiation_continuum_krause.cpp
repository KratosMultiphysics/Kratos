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
#include "radiation_continuum_krause.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  RadiationContinuumKrause::RadiationContinuumKrause() {}
  RadiationContinuumKrause::~RadiationContinuumKrause() {}

  //------------------------------------------------------------------------------------------------------------
  double RadiationContinuumKrause::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
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
  double RadiationContinuumKrause::FinalizeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check if radiation neighbors exist
    const int neighbors = particle->mRadiativeNeighbors;
    if (neighbors == 0)
      return 0.0;

    const double accum_temperature1   = particle->mEnvironmentTemperature;
    const double accum_temperature2   = particle->mEnvironmentTempAux;
    const double particle_emissivity  = particle->GetParticleEmissivity();
    const double particle_surface     = particle->GetParticleSurfaceArea();
    const double particle_temperature = particle->GetParticleTemperature();

    // Compute environment temperature
    const double env_temperature = pow(accum_temperature1 / accum_temperature2, 0.25);

    // Compute heat flux
    return STEFAN_BOLTZMANN * particle_emissivity * particle_surface * (pow(env_temperature,4.0) - pow(particle_temperature,4.0));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void RadiationContinuumKrause::AccumulateEnvironmentTemperature(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Update number of radiative neighbors
    particle->mRadiativeNeighbors++;

    // Accumulate temperature
    const double neighbor_emissivity   = particle->GetNeighborEmissivity();
    const double neighbor_temperature  = particle->GetNeighborTemperature();
    const double neighbor_surface      = particle->GetNeighborSurfaceArea();
    particle->mEnvironmentTemperature += 0.5 * STEFAN_BOLTZMANN * neighbor_emissivity * neighbor_surface * pow(neighbor_temperature, 4.0);
    particle->mEnvironmentTempAux     += 0.5 * STEFAN_BOLTZMANN * neighbor_emissivity * neighbor_surface;

    KRATOS_CATCH("")
  }

} // namespace Kratos
