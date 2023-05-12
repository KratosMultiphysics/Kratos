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
#include "real_contact_lu.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RealContactLu::RealContactLu() {}
  RealContactLu::~RealContactLu() {}

  //------------------------------------------------------------------------------------------------------------
  void RealContactLu::AdjustContact(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Effective radius
    const double eff_radius = particle->ComputeEffectiveRadius();

    // Simulation and real values of effective Young modulus
    const double eff_young      = particle->ComputeEffectiveYoung();
    const double eff_young_real = particle->ComputeEffectiveYoungReal();

    // Simulation and real values of stiffness
    const double stiff      = 4.0 / 3.0 * sqrt(eff_radius) * eff_young;
    const double stiff_real = 4.0 / 3.0 * sqrt(eff_radius) * eff_young_real;

    // Adjusted value of contact radius
    particle->mContactRadiusAdjusted = pow(particle->mContactRadius * stiff / stiff_real, 2.0/3.0);

    // Compute adjusted distance/separation from adjusted contact radius
    particle->mNeighborDistanceAdjusted   = particle->ComputeDistanceToNeighborAdjusted();
    particle->mNeighborSeparationAdjusted = particle->ComputeSeparationToNeighborAdjusted();

    KRATOS_CATCH("")
  }

} // namespace Kratos
