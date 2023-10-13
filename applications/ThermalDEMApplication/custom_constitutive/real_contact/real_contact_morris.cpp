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
#include "real_contact_morris.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RealContactMorris::RealContactMorris() {}
  RealContactMorris::~RealContactMorris() {}

  //------------------------------------------------------------------------------------------------------------
  void RealContactMorris::AdjustContact(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Parameters
    const double eff_young      = particle->ComputeEffectiveYoung();
    const double eff_young_real = particle->ComputeEffectiveYoungReal();
    const double eff_radius     = particle->ComputeEffectiveRadius();
    const double identation     = std::max(-1.0 * particle->mNeighborSeparation, 0.0);

    // Contact force with simulation parameters (using Hertz theory)
    const double hertz_force = 4.0 * eff_young * sqrt(eff_radius) * pow(identation, 3.0 / 2.0) / 3.0;

    // Area correction
    const double correction_area = pow(hertz_force * eff_radius / eff_young_real, 1.0 / 3.0);

    // Time correction
    double correction_time = 1.0; // TODO: Compute time correction from collision time

    // Adjusted value of contact radius
    particle->mContactRadiusAdjusted = correction_area * correction_time;

    // Compute adjusted distance/separation from adjusted contact radius
    particle->mNeighborDistanceAdjusted   = particle->ComputeDistanceToNeighborAdjusted();
    particle->mNeighborSeparationAdjusted = particle->ComputeSeparationToNeighborAdjusted();

    KRATOS_CATCH("")
  }

} // namespace Kratos
