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
#include "friction_coulomb.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  FrictionCoulomb::FrictionCoulomb() {}
  FrictionCoulomb::~FrictionCoulomb() {}

  //------------------------------------------------------------------------------------------------------------
  double FrictionCoulomb::ComputeHeatGeneration(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check for contact
    if (!particle->mNeighborInContact)
      return 0.0;

    typename ThermalSphericParticle:: ContactParams contact_params = particle->GetContactParameters();
    const double velocity_tangent = contact_params.local_velocity[1];
    const double force_normal     = contact_params.local_force[0];

    if (velocity_tangent == 0 || force_normal == 0)
      return 0.0;

    const double friction_conversion = r_process_info[FRICTION_HEAT_CONVERSION];
    const double friction_coeff      = particle->GetContactDynamicFrictionCoefficient();

    // Partition coefficient
    const double partition = ComputePartitionCoeff(particle);
    
    // Compute frictional heat transfer
    return partition * friction_conversion * friction_coeff * fabs(velocity_tangent * force_normal);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double FrictionCoulomb::ComputePartitionCoeff(ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double k1 = particle->GetParticleConductivity();
    const double k2 = particle->GetNeighborConductivity();
    return k1 / (k1 + k2);

    KRATOS_CATCH("")
  }

} // namespace Kratos
