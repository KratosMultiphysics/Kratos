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
#include "generation_dissipation.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  GenerationDissipation::GenerationDissipation() {}
  GenerationDissipation::~GenerationDissipation() {}

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputeHeatGeneration(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check for contact
    if (!particle->mNeighborInContact)
      return 0.0;

    // Conversion and partition coefficients
    const double conversion = r_process_info[HEAT_GENERATION_RATIO];
    const double partition  = ComputePartitionCoeff(particle);

    // Add contribution from different sources of energy dissipation
    double heat_gen_sliding = 0.0;
    double heat_gen_rolling = 0.0;
    double heat_gen_damping = 0.0;

    if (r_process_info[GENERATION_SLIDING_OPTION]) {
      heat_gen_sliding = partition * conversion * ComputeHeatGenerationSlidingFriction(particle);

      if (particle->mNeighborType & PARTICLE_NEIGHBOR)
        particle->mGenerationHeatFlux_slid_particle += heat_gen_sliding;
      else if (particle->mNeighborType & WALL_NEIGHBOR)
        particle->mGenerationHeatFlux_slid_wall += heat_gen_sliding;
    }
    if (r_process_info[GENERATION_ROLLING_OPTION] && particle->Is(DEMFlags::HAS_ROTATION) && particle->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
      heat_gen_rolling = partition * conversion * ComputeHeatGenerationRollingFriction(particle);

      if (particle->mNeighborType & PARTICLE_NEIGHBOR)
        particle->mGenerationHeatFlux_roll_particle += heat_gen_rolling;
      else if (particle->mNeighborType & WALL_NEIGHBOR)
        particle->mGenerationHeatFlux_roll_wall += heat_gen_rolling;
    }
    if (r_process_info[GENERATION_DAMPING_OPTION]) {
      heat_gen_damping = partition * conversion * ComputeHeatGenerationDampingContact(particle);

      if (particle->mNeighborType & PARTICLE_NEIGHBOR)
        particle->mGenerationHeatFlux_damp_particle += heat_gen_damping;
      else if (particle->mNeighborType & WALL_NEIGHBOR)
        particle->mGenerationHeatFlux_damp_wall += heat_gen_damping;
    }

    return heat_gen_sliding + heat_gen_rolling + heat_gen_damping;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputeHeatGenerationSlidingFriction(ThermalSphericParticle* particle) {
    KRATOS_TRY

    typename ThermalSphericParticle:: ContactParams contact_params = particle->GetContactParameters();
    const double force_normal     = contact_params.local_force_total[0];
    const double velocity_tangent = contact_params.local_velocity[1];

    if (fabs(force_normal)     < std::numeric_limits<double>::epsilon() ||
        fabs(velocity_tangent) < std::numeric_limits<double>::epsilon())
      return 0.0;

    const double friction_coeff = particle->GetContactDynamicFrictionCoefficient();
    return friction_coeff * fabs(force_normal * velocity_tangent);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputeHeatGenerationRollingFriction(ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Total normal force
    typename ThermalSphericParticle:: ContactParams contact_params = particle->GetContactParameters();
    const double force_normal = fabs(contact_params.local_force_total[0]);

    // Relative angular velocity
    // ASSUMPTION: Angular velocity of walls is not considered!
    double rel_vel;

    if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
      array_1d<double, 3> rel_angular_velocity;
      noalias(rel_angular_velocity) = particle->GetParticleAngularVelocity() - particle->mNeighbor_p->GetParticleAngularVelocity();
      rel_vel = DEM_MODULUS_3(rel_angular_velocity);
    }
    else {
      rel_vel = DEM_MODULUS_3(particle->GetParticleAngularVelocity());
    }

    if (force_normal < std::numeric_limits<double>::epsilon() ||
        rel_vel      < std::numeric_limits<double>::epsilon())
      return 0.0;
    
    const double eff_radius    = particle->ComputeEffectiveRadius();
    const double rolling_coeff = particle->GetContactRollingFrictionCoefficient();

    return rel_vel * rolling_coeff * eff_radius * force_normal;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputeHeatGenerationDampingContact(ThermalSphericParticle* particle) {
    KRATOS_TRY

    typename ThermalSphericParticle:: ContactParams contact_params = particle->GetContactParameters();
    const double force_damp_normal  = fabs(contact_params.local_force_damping[0]);
    const double force_damp_tangent = fabs(contact_params.local_force_damping[1]);
    const double velocity_normal    = fabs(contact_params.local_velocity[0]);
    const double velocity_tangent   = fabs(contact_params.local_velocity[1]);

    return force_damp_normal * velocity_normal + force_damp_tangent * velocity_tangent;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationDissipation::ComputePartitionCoeff(ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double k1 = particle->GetParticleConductivity();
    const double k2 = particle->GetNeighborConductivity();
    return k1 / (k1 + k2);

    KRATOS_CATCH("")
  }

} // namespace Kratos
