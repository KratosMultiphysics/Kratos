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
#include "real_contact_rangel_area.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RealContactRangelArea::RealContactRangelArea() {}
  RealContactRangelArea::~RealContactRangelArea() {}

  //------------------------------------------------------------------------------------------------------------
  void RealContactRangelArea::AdjustContact(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Simulation and real values of effective Young modulus
    const double eff_young      = particle->ComputeEffectiveYoung();
    const double eff_young_real = particle->ComputeEffectiveYoungReal();

    // Area correction
    if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
      const double r1           = particle->GetParticleRadius();
      const double r2           = particle->GetNeighborRadius();
      const double d_sim        = particle->mNeighborDistance;
      const double overlap_sim  = r1 + r2 - d_sim;
      const double overlap_real = pow(eff_young/eff_young_real, 2.0/3.0) * overlap_sim;
      const double d_real       = r1 + r2 - overlap_real;
      particle->mContactRadiusAdjusted = sqrt(fabs(r1 * r1 - pow(((r1 * r1 - r2 * r2 + d_real * d_real) / (2.0 * d_real)), 2.0)));
    }
    else if (particle->mNeighborType & WALL_NEIGHBOR) {
      const double r1           = particle->GetParticleRadius();
      const double d_sim        = particle->mNeighborDistance;
      const double overlap_sim  = r1 - d_sim;
      const double overlap_real = pow(eff_young/eff_young_real, 2.0/3.0) * overlap_sim;
      const double d_real       = r1 - overlap_real;
      particle->mContactRadiusAdjusted = sqrt(r1 * r1 - d_real * d_real);
    }
    
    // Compute adjusted distance/separation from adjusted contact radius
    particle->mNeighborDistanceAdjusted   = particle->ComputeDistanceToNeighborAdjusted();
    particle->mNeighborSeparationAdjusted = particle->ComputeSeparationToNeighborAdjusted();

    KRATOS_CATCH("")
  }

} // namespace Kratos
