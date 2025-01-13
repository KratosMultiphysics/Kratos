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
#include "direct_conduction_collision.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  DirectConductionCollision::DirectConductionCollision() {}
  DirectConductionCollision::~DirectConductionCollision() {}

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionCollision::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check for contact
    if (!particle->mNeighborInContact)
      return 0.0;

    // Get collision time and impact normal velocity
    typename ThermalSphericParticle::ContactParams contact_params = particle->GetContactParameters();
    const double col_time = r_process_info[TIME] - contact_params.impact_time;
    const double impact_normal_velocity = std::abs(contact_params.impact_velocity[0]);

    // Compute max collision time
    double col_time_max = 0.0;
    if (impact_normal_velocity != 0.0)
      col_time_max = particle->ComputeMaxCollisionTime();
    
    // Check if collision time is smaller than max value, otherwise use static model (simple BOB)
    if (col_time < col_time_max) {
      const double temp_grad = particle->GetNeighborTemperature() - particle->GetParticleTemperature();
      const double Rc_max    = particle->ComputeMaxContactRadius(); // TODO: This should be multiplied by the correction coefficient (and not computed with real Young modulus)
      const double Fo        = particle->ComputeFourierNumber();

      const double a1 = particle->GetParticleDensity() * particle->GetParticleHeatCapacity();
      const double a2 = particle->GetNeighborDensity() * particle->GetNeighborHeatCapacity();
      const double b1 = a1 * particle->GetParticleConductivity();
      const double b2 = a2 * particle->GetNeighborConductivity();
      const double c  = a1 / a2;

      const double C1 = -2.300 * c * c +  8.909 * c - 4.235;
      const double C2 =  8.169 * c * c - 33.770 * c + 24.885;
      const double C3 = -5.758 * c * c + 24.464 * c - 20.511;

      const double C_coeff = 0.435 * (sqrt(C2 * C2 - 4.0 * C1 * (C3 - Fo)) - C2) / C1;

      return C_coeff * Globals::Pi * Rc_max * Rc_max * pow(col_time_max,-0.5) * temp_grad / (pow(b1,-0.5) + pow(b2,-0.5));
    }
    else {
      // Assumption: Use simple BOB model for static contact
      const double Rc        = particle->mContactRadiusAdjusted;
      const double keff      = particle->ComputeEffectiveConductivity();
      const double temp_grad = particle->GetNeighborTemperature() - particle->GetParticleTemperature();
      return 4.0 * keff * Rc * temp_grad;
    }

    KRATOS_CATCH("")
  }

} // namespace Kratos
