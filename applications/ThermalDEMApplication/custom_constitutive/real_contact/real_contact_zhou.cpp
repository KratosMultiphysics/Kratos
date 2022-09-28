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
#include "real_contact_zhou.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RealContactZhou::RealContactZhou() {}
  RealContactZhou::~RealContactZhou() {}

  //------------------------------------------------------------------------------------------------------------
  void RealContactZhou::AdjustContact(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Simulation and real values of effective Young modulus
    const double eff_young      = particle->ComputeEffectiveYoung();
    const double eff_young_real = particle->ComputeEffectiveYoungReal();

    // Adjusted value of contact radius
    particle->mContactRadiusAdjusted = particle->mContactRadius * pow(eff_young / eff_young_real, 0.2);

    // Compute adjusted distance/separation from adjusted contact radius
    particle->mNeighborDistanceAdjusted   = particle->ComputeDistanceToNeighborAdjusted();
    particle->mNeighborSeparationAdjusted = particle->ComputeSeparationToNeighborAdjusted();

    KRATOS_CATCH("")
  }

} // namespace Kratos
