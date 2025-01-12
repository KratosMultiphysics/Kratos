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
    double thermal_energy;
    double heat_gen;
    double heat_gen_damping_pp = 0.0;
    double heat_gen_damping_pw = 0.0;
    double heat_gen_sliding_pp = 0.0;
    double heat_gen_sliding_pw = 0.0;
    double heat_gen_rolling_pp = 0.0;
    double heat_gen_rolling_pw = 0.0;

    // Check type of wall neighbor
    bool is_wall_side = false;
    if (particle->mNeighborType & WALL_NEIGHBOR) {
      DEMWall* neighbor = particle->mNeighbor_w;
      array_1d<double, 3> coord_z = ZeroVector(3);
      int number_of_nodes = neighbor->GetGeometry().size();
      for (unsigned int i = 0; i < number_of_nodes; i++) {
        const array_1d<double, 3>& Nodecoord = neighbor->GetGeometry()[i].Coordinates();
        coord_z[i] = Nodecoord[2];
      }
      if (fabs(coord_z[0] - coord_z[1]) < 0.00001 && fabs(coord_z[0] - coord_z[2]) < 0.00001) {
        is_wall_side = true;
      }
    }

    // Damping thermal power
    if (r_process_info[GENERATION_DAMPING_OPTION]) {
      thermal_energy = coeff * contact_params.viscodamping_energy;
      heat_gen = thermal_energy / time;

      if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
        particle->mGenerationThermalEnergy_damp_particle += thermal_energy;
        particle->mGenerationHeatFlux_damp_particle      += heat_gen;
        heat_gen_damping_pp = heat_gen;
      }
      else if (particle->mNeighborType & WALL_NEIGHBOR) {
        particle->mGenerationThermalEnergy_damp_wall += thermal_energy;
        particle->mGenerationHeatFlux_damp_wall      += heat_gen;
        heat_gen_damping_pw = heat_gen;

        if (is_wall_side)
          particle->mGenerationThermalEnergy_wall_sides += thermal_energy;
        else
          particle->mGenerationThermalEnergy_wall_annul += thermal_energy;
      }  
    }

    // Sliding friction thermal power
    if (r_process_info[GENERATION_SLIDING_OPTION]) {
      thermal_energy = coeff * contact_params.frictional_energy;
      heat_gen = thermal_energy / time;

      if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
        particle->mGenerationThermalEnergy_slid_particle += thermal_energy;
        particle->mGenerationHeatFlux_slid_particle      += heat_gen;
        heat_gen_sliding_pp = heat_gen;
      }
      else if (particle->mNeighborType & WALL_NEIGHBOR) {
        particle->mGenerationThermalEnergy_slid_wall += thermal_energy;
        particle->mGenerationHeatFlux_slid_wall      += heat_gen;
        heat_gen_sliding_pw = heat_gen;

        if (is_wall_side)
          particle->mGenerationThermalEnergy_wall_sides += thermal_energy;
        else
          particle->mGenerationThermalEnergy_wall_annul += thermal_energy;
      }
    }

    // Rolling friction thermal power
    if (r_process_info[GENERATION_ROLLING_OPTION] && particle->Is(DEMFlags::HAS_ROTATION) && particle->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
      thermal_energy = coeff * contact_params.rollresist_energy;
      heat_gen = thermal_energy / time;

      if (particle->mNeighborType & PARTICLE_NEIGHBOR) {
        particle->mGenerationThermalEnergy_roll_particle += thermal_energy;
        particle->mGenerationHeatFlux_roll_particle      += heat_gen;
        heat_gen_rolling_pp = heat_gen;
      }
      else if (particle->mNeighborType & WALL_NEIGHBOR) {
        particle->mGenerationThermalEnergy_roll_wall += thermal_energy;
        particle->mGenerationHeatFlux_roll_wall      += heat_gen;
        heat_gen_rolling_pw = heat_gen;

        if (is_wall_side)
          particle->mGenerationThermalEnergy_wall_sides += thermal_energy;
        else
          particle->mGenerationThermalEnergy_wall_annul += thermal_energy;
      }
    }

    // Fill heat map
    if (r_process_info[HEAT_MAP_GENERATION_OPTION])
      FillHeatMap(r_process_info, particle, time, heat_gen_damping_pp, heat_gen_damping_pw, heat_gen_sliding_pp, heat_gen_sliding_pw, heat_gen_rolling_pp, heat_gen_rolling_pw);

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
  void GenerationDissipation::FillHeatMap(const ProcessInfo& r_process_info,
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
