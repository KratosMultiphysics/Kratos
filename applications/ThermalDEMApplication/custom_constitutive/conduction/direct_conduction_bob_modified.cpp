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
#include "direct_conduction_bob_modified.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  DirectConductionBOBModified::DirectConductionBOBModified() {}
  DirectConductionBOBModified::~DirectConductionBOBModified() {}

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBModified::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Compute interaction thermal resistance
    const double r_bob = 1.0 / DirectConductionBOBComplete::ComputeHeatTransferCoeff(r_process_info, particle);

    // Compute inner thermal resistance
    const double r_inn = ComputeInnerResistance(r_process_info, particle);

    // Compute effective thermal resistance
    const double r_eff = r_bob + r_inn;

    // Compute heat flux
    return (particle->GetNeighborTemperature() - particle->GetParticleTemperature()) / r_eff;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBModified::ComputeInnerResistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Heat transfer resistance for direct and/or indirect conduction
    if (particle->mNeighborType & PARTICLE_NEIGHBOR)
      return InnerResistanceWithParticle(r_process_info, particle);
    else if (particle->mNeighborType & WALL_NEIGHBOR)
      return InnerResistanceWithWall(r_process_info, particle);
    else
      return 0.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBModified::InnerResistanceWithParticle(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double kp    = particle->GetParticleConductivity();
    const double kn    = particle->mNeighbor_p->GetParticleConductivity();
    const double Rp    = particle->GetParticleRadius();
    const double Rn    = particle->mNeighbor_p->GetParticleRadius();
    const double Reff  = particle->ComputeEffectiveRadius();
    const double Rcond = r_process_info[CONDUCTION_RADIUS];

    // Conductive cylinder radius
    const double Rcyl = Rcond * Reff;

    // Inner thermal resistances
    const double rp = Rp / kp * Globals::Pi * Rcyl * Rcyl;
    const double rn = Rn / kn * Globals::Pi * Rcyl * Rcyl;

    // Effective inner resistance
    return rp + rn;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBModified::InnerResistanceWithWall(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double kp    = particle->GetParticleConductivity();
    const double Rp    = particle->GetParticleRadius();
    const double Reff  = particle->ComputeEffectiveRadius();
    const double Rcond = r_process_info[CONDUCTION_RADIUS];

    // Conductive cylinder radius
    const double Rcyl = Rcond * Reff;

    // Inner thermal resistances (particle only)
    return Rp / kp * Globals::Pi * Rcyl * Rcyl;

    KRATOS_CATCH("")
  }

} // namespace Kratos
