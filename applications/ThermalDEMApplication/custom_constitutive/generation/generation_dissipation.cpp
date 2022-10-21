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

    // Get contact info with current neighbor
    typename ThermalSphericParticle::ContactParams contact_params = particle->GetContactParameters();

    // Time passed since last thermal solution
    double time = particle->mNumStepsEval * r_process_info[DELTA_TIME];

    // Conversion and partition coefficients
    const double conversion = r_process_info[HEAT_GENERATION_RATIO];
    const double partition  = ComputePartitionCoeff(particle);
    const double coeff      = conversion * partition;

    // Initialize contribution from different sources of heat generation
    double heat_gen;
    double heat_gen_damping_pp = 0.0;
    double heat_gen_damping_pw = 0.0;
    double heat_gen_sliding_pp = 0.0;
    double heat_gen_sliding_pw = 0.0;
    double heat_gen_rolling_pp = 0.0;
    double heat_gen_rolling_pw = 0.0;

    // Damping thermal power
    if (r_process_info[GENERATION_DAMPING_OPTION]) {
      heat_gen = coeff * contact_params.viscodamping_energy / time;

      if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
        heat_gen_damping_pp = heat_gen;
        particle->mGenerationHeatFlux_damp_particle += heat_gen;
      }
      else if (particle->mNeighborType & WALL_NEIGHBOR) {
        heat_gen_damping_pw = heat_gen;
        particle->mGenerationHeatFlux_damp_wall += heat_gen;
      }  
    }

    // Sliding friction thermal power
    if (r_process_info[GENERATION_SLIDING_OPTION]) {
      heat_gen = coeff * contact_params.frictional_energy / time;

      if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
        heat_gen_sliding_pp = heat_gen;
        particle->mGenerationHeatFlux_slid_particle += heat_gen;
      }
      else if (particle->mNeighborType & WALL_NEIGHBOR) {
        heat_gen_sliding_pw = heat_gen;
        particle->mGenerationHeatFlux_slid_wall += heat_gen;
      }
    }

    // Rolling friction thermal power
    if (r_process_info[GENERATION_ROLLING_OPTION] && particle->Is(DEMFlags::HAS_ROTATION) && particle->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
      heat_gen = coeff * contact_params.rollresist_energy / time;

      if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
        heat_gen_rolling_pp = heat_gen;
        particle->mGenerationHeatFlux_roll_particle += heat_gen;
      }
      else if (particle->mNeighborType & WALL_NEIGHBOR) {
        heat_gen_rolling_pw = heat_gen;
        particle->mGenerationHeatFlux_roll_wall += heat_gen;
      }
    }

    // Fill heat map
    if (r_process_info[HEAT_MAP_GENERATION_OPTION])
      FillDensityMap(r_process_info, particle, time, heat_gen_damping_pp, heat_gen_damping_pw, heat_gen_sliding_pp, heat_gen_sliding_pw, heat_gen_rolling_pp, heat_gen_rolling_pw);

    return heat_gen_damping_pp + heat_gen_damping_pw + heat_gen_sliding_pp + heat_gen_sliding_pw + heat_gen_rolling_pp + heat_gen_rolling_pw;

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

  //------------------------------------------------------------------------------------------------------------
  void GenerationDissipation::FillDensityMap(const ProcessInfo& r_process_info,
                                             ThermalSphericParticle* particle,
                                             const double time,
                                             const double heat_gen_damping_pp,
                                             const double heat_gen_damping_pw,
                                             const double heat_gen_sliding_pp,
                                             const double heat_gen_sliding_pw,
                                             const double heat_gen_rolling_pp,
                                             const double heat_gen_rolling_pw) {
    KRATOS_TRY

    // Get heat generation coordinates (contact point)
    const array_1d<double, 3> coord = particle->GetParticleCoordinates();
    array_1d<double, 3> dir;
    noalias(dir) = particle->GetNeighborCoordinates() - coord;

    const double dist        = DEM_MODULUS_3(dir);
    const double r1          = particle->GetRadius();
    const double r2          = particle->GetRadius();
    const double indentation = r1 + r2 - dist;
    const double length      = r1 - indentation / 2.0;
    const double ratio       = length / dist;

    const double x = coord[0] + dir[0] * ratio;
    const double y = coord[1] + dir[1] * ratio;
    const double z = coord[2] + dir[2] * ratio;

    // Get map geometry (corner points must have already been adjusted so that all coords_1 < coords_2)
    const array_1d<double, 3> coords_1     = r_process_info[HEAT_MAP_COORDINATES_1];
    const array_1d<double, 3> coords_2     = r_process_info[HEAT_MAP_COORDINATES_2];
    const array_1d<int, 3>    subdivisions = r_process_info[HEAT_MAP_SUBDIVISIONS];

    const double dx = (coords_2[0] - coords_1[0]) / subdivisions[0];
    const double dy = (coords_2[1] - coords_1[1]) / subdivisions[1];
    const double dz = (coords_2[2] - coords_1[2]) / subdivisions[2];

    // Check if generation point is inside map boundary
    if ((x < coords_1[0]) || (x > coords_2[0]) ||
        (y < coords_1[1]) || (y > coords_2[1]) ||
        (z < coords_1[2]) || (z > coords_2[2]))
      return;

    // Determine indexes
    const int i_x = (x - coords_1[0]) / dx;
    const int i_y = (y - coords_1[1]) / dy;
    const int i_z = (z - coords_1[2]) / dz;
    
    // Fill maps with energy (heat power * time)
    particle->mHeatMapGenerationDampingPP[i_x][i_y][i_z] += heat_gen_damping_pp * time;
    particle->mHeatMapGenerationDampingPW[i_x][i_y][i_z] += heat_gen_damping_pw * time;
    particle->mHeatMapGenerationSlidingPP[i_x][i_y][i_z] += heat_gen_sliding_pp * time;
    particle->mHeatMapGenerationSlidingPW[i_x][i_y][i_z] += heat_gen_sliding_pw * time;
    particle->mHeatMapGenerationRollingPP[i_x][i_y][i_z] += heat_gen_rolling_pp * time;
    particle->mHeatMapGenerationRollingPW[i_x][i_y][i_z] += heat_gen_rolling_pw * time;

    KRATOS_CATCH("")
  }

} // namespace Kratos
