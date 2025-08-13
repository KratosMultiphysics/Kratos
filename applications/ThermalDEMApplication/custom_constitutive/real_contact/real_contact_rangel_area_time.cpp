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
#include "real_contact_rangel_area_time.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  RealContactRangelAreaTime::RealContactRangelAreaTime() {}
  RealContactRangelAreaTime::~RealContactRangelAreaTime() {}

  //------------------------------------------------------------------------------------------------------------
  void RealContactRangelAreaTime::AdjustContact(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Simulation and real values of effective Young modulus
    const double eff_young      = particle->ComputeEffectiveYoung();
    const double eff_young_real = particle->ComputeEffectiveYoungReal();
    
    // Get collision time
    typename ThermalSphericParticle::ContactParams contact_params = particle->GetContactParameters();
    const double col_time          = r_process_info[TIME] - contact_params.impact_time;
    const double col_time_max      = particle->ComputeMaxCollisionTime();
    const double col_time_max_real = particle->ComputeMaxCollisionTimeReal();

    // Area correction
    if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
      const double r1           = particle->GetParticleRadius();
      const double r2           = particle->GetNeighborRadius();
      const double d_sim        = particle->mNeighborDistance;
      const double overlap_sim  = r1 + r2 - d_sim;
      const double overlap_real = pow(eff_young/eff_young_real, 2.0/3.0) * overlap_sim;
      const double d_real       = r1 + r2 - overlap_real;
      particle->mContactRadiusAdjusted = sqrt(std::abs(r1 * r1 - pow(((r1 * r1 - r2 * r2 + d_real * d_real) / (2.0 * d_real)), 2.0)));
    }
    else if (particle->mNeighborType & WALL_NEIGHBOR) {
      const double r1           = particle->GetParticleRadius();
      const double d_sim        = particle->mNeighborDistance;
      const double overlap_sim  = r1 - d_sim;
      const double overlap_real = pow(eff_young/eff_young_real, 2.0/3.0) * overlap_sim;
      const double d_real       = r1 - overlap_real;
      particle->mContactRadiusAdjusted = sqrt(r1 * r1 - d_real * d_real);
    }

    // Time correction
    double correction_time = 1.0;
    if (col_time < col_time_max) {
      correction_time = pow(col_time_max_real / col_time_max, 2.0 / 3.0);
      particle->mContactRadiusAdjusted *= correction_time;
    }
    
    // Compute adjusted distance/separation from adjusted contact radius
    particle->mNeighborDistanceAdjusted   = particle->ComputeDistanceToNeighborAdjusted();
    particle->mNeighborSeparationAdjusted = particle->ComputeSeparationToNeighborAdjusted();

    KRATOS_CATCH("")
  }

} // namespace Kratos
