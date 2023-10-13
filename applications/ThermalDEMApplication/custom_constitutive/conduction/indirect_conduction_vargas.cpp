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
#include "indirect_conduction_vargas.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  IndirectConductionVargas::IndirectConductionVargas() {}
  IndirectConductionVargas::~IndirectConductionVargas() {}

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVargas::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check for contact
    if (!particle->mNeighborInContact)
      return 0.0;

    // Assumption 1: Formulation for a liquid (not gas) as the interstitial fluid is being used
    const double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double temp_grad          = particle->GetNeighborTemperature() - particle->GetParticleTemperature();

    // Assumption 2 : Model developed for mono-sized particles, but the average radius is being used (if neighbor is a wall, it is assumed as a particle with the same radius)
    const double particle_radius = particle->GetParticleRadius();
    const double neighbor_radius = (particle->mNeighborType & PARTICLE_NEIGHBOR) ? particle->GetNeighborRadius() : particle_radius;
    const double avg_radius      = (particle_radius + neighbor_radius) / 2.0;
    const double contact_radius  = particle->mContactRadiusAdjusted;

    // Compute heat flux
    return 2.0 * Globals::Pi * fluid_conductivity * temp_grad * (1.0 - 0.5 * pow(contact_radius / avg_radius, 2.0)) * (avg_radius - contact_radius) / (1.0 - Globals::Pi / 4.0);

    KRATOS_CATCH("")
  }

} // namespace Kratos
