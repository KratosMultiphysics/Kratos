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
#include "direct_conduction_bob_complete.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  DirectConductionBOBComplete::DirectConductionBOBComplete() {}
  DirectConductionBOBComplete::~DirectConductionBOBComplete() {}

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBComplete::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return particle->GetParticleRadius() * r_process_info[MAX_CONDUCTION_DISTANCE];
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBComplete::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Compute heat transfer coefficient
    const double h = ComputeHeatTransferCoeff(r_process_info, particle);

    // Compute heat flux
    return h * (particle->GetNeighborTemperature() - particle->GetParticleTemperature());

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBComplete::ComputeHeatTransferCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Heat transfer coefficient for direct and/or indirect conduction
    if (particle->mNeighborInContact)
      return ContactCoeff(r_process_info, particle);
    else if (particle->mNeighborSeparation < r_process_info[MAX_CONDUCTION_DISTANCE] * particle->ComputeEffectiveRadius())
      return SeparatedCoeff(r_process_info, particle);
    else
      return 0.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBComplete::ContactCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double Rc   = particle->mContactRadiusAdjusted;
    const double Reff = particle->ComputeEffectiveRadius();
    const double kp   = particle->ComputeMeanConductivity();  // Assumption: average conductivity for different properties
    const double kf   = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double kpf  = kp / kf;

    // Dimentionless parameter
    const double eta = kpf * Rc / (2.0 * Reff);

    // Compute heat transfer coefficient
    if (eta <= ETA_MIN) {
      return 2.0 * Globals::Pi * Reff * kf * (0.17 * eta * eta + log(kpf * kpf));
    }
    else if (eta >= ETA_MAX) {
      return 4.0 * Globals::Pi * Reff * kf * (eta / Globals::Pi + log(kpf / eta));
    }
    else { // Linear interpolation
      const double a = 2.0 * Globals::Pi * Reff * kf * (0.17 * eta * eta + log(kpf * kpf));
      const double b = 4.0 * Globals::Pi * Reff * kf * (eta / Globals::Pi + log(kpf / eta));
      return a + (eta - ETA_MIN) * (b - a) / (ETA_MAX - ETA_MIN);
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionBOBComplete::SeparatedCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double Ds    = particle->mNeighborSeparation;
    const double Reff  = particle->ComputeEffectiveRadius();
    const double Rcond = r_process_info[CONDUCTION_RADIUS];
    const double kp    = particle->ComputeMeanConductivity();  // Assumption: average conductivity for different properties
    const double kf    = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double kpf2  = (kp / kf) * (kp / kf);

    // Conductive cylinder radius
    const double Rcyl = Rcond * Reff;

    // Dimentionless parameter
    const double lamb = kpf2 * Ds / (2.0 * Reff);

    // Compute heat transfer coefficient
    if (lamb <= LAMB_MIN) {
      return 2.0 * Globals::Pi * Reff * kf * log(kpf2);
    }
    else if (lamb >= LAMB_MAX) {
      return 2.0 * Globals::Pi * Reff * kf * log(1.0 + Rcyl * Rcyl / (2.0 * Reff * Ds));
    }
    else { // Minimum value
      const double a = 2.0 * Globals::Pi * Reff * kf * log(kpf2);
      const double b = 2.0 * Globals::Pi * Reff * kf * log(1.0 + Rcyl * Rcyl / (2.0 * Reff * Ds));
      return std::min(a,b);
    }

    KRATOS_CATCH("")
  }

} // namespace Kratos
