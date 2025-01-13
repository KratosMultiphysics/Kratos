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
#include "graph_utilities.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  GraphUtilities::GraphUtilities() {
    mGraph_ParticleTempMin       = false;
    mGraph_ParticleTempMax       = false;
    mGraph_ParticleTempAvg       = false;
    mGraph_ParticleTempAvgVol    = false;
    mGraph_ParticleTempDev       = false;
    mGraph_EnergyMechanical      = false;
    mGraph_EnergyDissipated      = false;
    mGraph_EnergyThermal         = false;
    mGraph_HeatFluxContributions = false;
    mGraph_HeatGenValues         = false;
    mGraph_HeatGenContributions  = false;
  }

  GraphUtilities::~GraphUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteInitialize(bool ParticleTempAll,
                                         bool ParticleTempMin,
                                         bool ParticleTempMax,
                                         bool ParticleTempAvg,
                                         bool ParticleTempAvgVol,
                                         bool ParticleTempDev,
                                         bool EnergyMechanical,
                                         bool EnergyDissipated,
                                         bool EnergyThermal,
                                         bool HeatFluxContributions,
                                         bool HeatGenValues,
                                         bool HeatGenContributions)
  {
    KRATOS_TRY

    // Set member flags
    mGraph_ParticleTempAll               = ParticleTempAll;
    mGraph_ParticleTempMin       = ParticleTempMin;
    mGraph_ParticleTempMax       = ParticleTempMax;
    mGraph_ParticleTempAvg       = ParticleTempAvg;
    mGraph_ParticleTempAvgVol    = ParticleTempAvgVol;
    mGraph_ParticleTempDev       = ParticleTempDev;
    mGraph_EnergyMechanical      = EnergyMechanical;
    mGraph_EnergyDissipated      = EnergyDissipated;
    mGraph_EnergyThermal         = EnergyThermal;
    mGraph_HeatFluxContributions = HeatFluxContributions;
    mGraph_HeatGenValues         = HeatGenValues;
    mGraph_HeatGenContributions  = HeatGenContributions;

    // Open files
    OpenFiles();

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteFinalizeSolutionStep(ModelPart& rModelPart) {
    KRATOS_TRY
    
    if (rModelPart.GetProcessInfo()[IS_TIME_TO_PRINT])
      WriteGraphs(rModelPart);

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteFinalize(void) {
    KRATOS_TRY

    CloseFiles();

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::OpenFiles(void) {
    if (mGraph_ParticleTempAll) {
      mFile_ParticleTempAll.open("graph_particle_temp_all.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempAll) << "Could not open graph file graph_particle_temp_all!" << std::endl;
      mFile_ParticleTempAll << "1 - TIME STEP / TIME | ";
      mFile_ParticleTempAll << "2 - ID / RADIUS / COORDS XY / TEMPERATURE";
      mFile_ParticleTempAll << std::endl;
    }
    if (mGraph_ParticleTempMin) {
      mFile_ParticleTempMin.open("graph_particle_temp_min.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempMin) << "Could not open graph_particle_temp_min.txt!" << std::endl;
      mFile_ParticleTempMin << "1 - TIME STEP | ";
      mFile_ParticleTempMin << "2 - TIME | ";
      mFile_ParticleTempMin << "3 - MIN PARTICLE TEMPERATURE";
      mFile_ParticleTempMin << std::endl;
    }
    if (mGraph_ParticleTempMax) {
      mFile_ParticleTempMax.open("graph_particle_temp_max.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempMax) << "Could not open graph_particle_temp_max.txt!" << std::endl;
      mFile_ParticleTempMax << "1 - TIME STEP | ";
      mFile_ParticleTempMax << "2 - TIME | ";
      mFile_ParticleTempMax << "3 - MAX PARTICLE TEMPERATURE";
      mFile_ParticleTempMax << std::endl;
    }
    if (mGraph_ParticleTempAvg) {
      mFile_ParticleTempAvg.open("graph_particle_temp_avg.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempAvg) << "Could not open graph_particle_temp_avg.txt!" << std::endl;
      mFile_ParticleTempAvg << "1 - TIME STEP | ";
      mFile_ParticleTempAvg << "2 - TIME | ";
      mFile_ParticleTempAvg << "3 - AVERAGE PARTICLE TEMPERATURE";
      mFile_ParticleTempAvg << std::endl;
    }
    if (mGraph_ParticleTempAvgVol) {
      mFile_ParticleTempAvgVol.open("graph_particle_temp_avg_vol.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempAvgVol) << "Could not open graph_particle_temp_avg_vol.txt!" << std::endl;
      mFile_ParticleTempAvgVol << "1 - TIME STEP | ";
      mFile_ParticleTempAvgVol << "2 - TIME | ";
      mFile_ParticleTempAvgVol << "3 - AVERAGE PARTICLE TEMPERATURE (VOLUME WEIGHTED)";
      mFile_ParticleTempAvgVol << std::endl;
    }
    if (mGraph_ParticleTempDev) {
      mFile_ParticleTempDev.open("graph_particle_temp_dev.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempDev) << "Could not open graph_particle_temp_dev.txt!" << std::endl;
      mFile_ParticleTempDev << "1 - TIME STEP | ";
      mFile_ParticleTempDev << "2 - TIME | ";
      mFile_ParticleTempDev << "3 - PARTICLE TEMPERATURE STANDARD DEVIATION";
      mFile_ParticleTempDev << std::endl;
    }
    if (mGraph_EnergyMechanical) {
      mFile_EnergyMechanical.open("graph_energy_mechanical.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_EnergyMechanical) << "Could not open graph_energy_mechanical.txt!" << std::endl;
      mFile_EnergyMechanical << "1 - TIME STEP | ";
      mFile_EnergyMechanical << "2 - TIME | ";
      mFile_EnergyMechanical << "3 - GRAVITY | ";
      mFile_EnergyMechanical << "4 - ELASTIC | ";
      mFile_EnergyMechanical << "5 - KINETIC TRANSLATION | ";
      mFile_EnergyMechanical << "6 - KINETIC ROTATION | ";
      mFile_EnergyMechanical << "7 - TOTAL";
      mFile_EnergyMechanical << std::endl;
    }
    if (mGraph_EnergyDissipated) {
      mFile_EnergyDissipated.open("graph_energy_dissipated.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_EnergyDissipated) << "Could not open graph_energy_dissipated.txt!" << std::endl;
      mFile_EnergyDissipated << "1 - TIME STEP | ";
      mFile_EnergyDissipated << "2 - TIME | ";
      mFile_EnergyDissipated << "3 - SLIDING | ";
      mFile_EnergyDissipated << "4 - ROLLING | ";
      mFile_EnergyDissipated << "5 - DAMPING | ";
      mFile_EnergyDissipated << "6 - TOTAL (ACCUMULATED)";
      mFile_EnergyDissipated << std::endl;
    }
    if (mGraph_EnergyThermal) {
      mFile_EnergyThermal.open("graph_energy_thermal.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_EnergyThermal) << "Could not open graph file graph_energy_thermal.txt!" << std::endl;
      mFile_EnergyThermal << "1 - TIME STEP | ";
      mFile_EnergyThermal << "2 - TIME | ";
      mFile_EnergyThermal << "3 - SLIDING PARTICLE-PARTICLE | ";
      mFile_EnergyThermal << "4 - SLIDING PARTICLE-WALL | ";
      mFile_EnergyThermal << "5 - ROLLING PARTICLE-PARTICLE | ";
      mFile_EnergyThermal << "6 - ROLLING PARTICLE-WALL | ";
      mFile_EnergyThermal << "7 - DAMPING PARTICLE-PARTICLE | ";
      mFile_EnergyThermal << "8 - DAMPING PARTICLE-WALL | ";
      mFile_EnergyThermal << "9 - TOTAL (ACCUMULATED)";
      mFile_EnergyThermal << std::endl;
    }
    if (mGraph_HeatFluxContributions) {
      mFile_HeatFluxContributions.open("graph_heat_flux_contributions.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_HeatFluxContributions) << "Could not open graph_heat_flux_contributions.txt!" << std::endl;
      mFile_HeatFluxContributions << "1 - TIME STEP | ";
      mFile_HeatFluxContributions << "2 - TIME | ";
      mFile_HeatFluxContributions << "3 - CONDUCTION DIRECT | ";
      mFile_HeatFluxContributions << "4 - CONDUCTION INDIRECT | ";
      mFile_HeatFluxContributions << "5 - RADIATION | ";
      mFile_HeatFluxContributions << "6 - GENERATION | ";
      mFile_HeatFluxContributions << "7 - CONVECTION | ";
      mFile_HeatFluxContributions << "8 - SURFACE PRESCRIBED | ";
      mFile_HeatFluxContributions << "9 - VOLUME PRESCRIBED";
      mFile_HeatFluxContributions << std::endl;
    }
    if (mGraph_HeatGenValues) {
      mFile_HeatGenValues.open("graph_heat_generation_values.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_HeatGenValues) << "Could not open graph_heat_generation_values.txt!" << std::endl;
      mFile_HeatGenValues << "1 - TIME STEP | ";
      mFile_HeatGenValues << "2 - TIME | ";
      mFile_HeatGenValues << "3 - SLIDING PARTICLE-PARTICLE | ";
      mFile_HeatGenValues << "4 - SLIDING PARTICLE-WALL | ";
      mFile_HeatGenValues << "5 - ROLLING PARTICLE-PARTICLE | ";
      mFile_HeatGenValues << "6 - ROLLING PARTICLE-WALL | ";
      mFile_HeatGenValues << "7 - DAMPING PARTICLE-PARTICLE | ";
      mFile_HeatGenValues << "8 - DAMPING PARTICLE-WALL";
      mFile_HeatGenValues << "9 - TOTAL";
      mFile_HeatGenValues << std::endl;
    }
    if (mGraph_HeatGenContributions) {
      mFile_HeatGenContributions.open("graph_heat_generation_contributions.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_HeatGenContributions) << "Could not open graph_heat_generation_contributions.txt!" << std::endl;
      mFile_HeatGenContributions << "1 - TIME STEP | ";
      mFile_HeatGenContributions << "2 - TIME | ";
      mFile_HeatGenContributions << "3 - SLIDING PARTICLE-PARTICLE | ";
      mFile_HeatGenContributions << "4 - SLIDING PARTICLE-WALL | ";
      mFile_HeatGenContributions << "5 - ROLLING PARTICLE-PARTICLE | ";
      mFile_HeatGenContributions << "6 - ROLLING PARTICLE-WALL | ";
      mFile_HeatGenContributions << "7 - DAMPING PARTICLE-PARTICLE | ";
      mFile_HeatGenContributions << "8 - DAMPING PARTICLE-WALL";
      mFile_HeatGenContributions << std::endl;
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::WriteGraphs(ModelPart& rModelPart) {
    KRATOS_TRY

    // Initialize results
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    // Print all temperatures
    if (mFile_ParticleTempAll.is_open() && time_step % write_all_temp_freq == 0.0) {
      mFile_ParticleTempAll << "#TIME/STEP: ";
      mFile_ParticleTempAll << std::defaultfloat << time_step << " " << time << std::endl;

      ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();
      for (int i = 0; i < num_of_particles; i++) {
        ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*(it + i));
        const int    id = i + 1;
        const double r  = particle.GetRadius();
        const double x  = particle.GetGeometry()[0][0];
        const double y  = particle.GetGeometry()[0][1];
        const double t  = particle.GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
        mFile_ParticleTempAll << std::defaultfloat << id << " ";
        mFile_ParticleTempAll << std::fixed << std::setprecision(15) << r << " " << x << " " << y << " " << t << std::endl;
      }
      mFile_ParticleTempAll << std::endl;
    }

    const int num_particles              = rModelPart.NumberOfElements();
    int       num_particles_flux_contrib = 0;
    int       num_particles_gen_contrib  = 0;
    double    vol_total                  = 0.0;

    double    particle_temp_min     =  DBL_MAX;
    double    particle_temp_max     = -DBL_MAX;
    double    particle_temp_avg     =  0.0;
    double    particle_temp_avg_vol =  0.0;
    double    particle_temp_dev     =  0.0;

    double    energy_mech_gravity_all     = 0.0;
    double    energy_mech_elastic_all     = 0.0;
    double    energy_mech_translation_all = 0.0;
    double    energy_mech_rotation_all    = 0.0;
    double    energy_mech_all             = 0.0;

    double    energy_dissip_slid_all = 0.0;
    double    energy_dissip_roll_all = 0.0;
    double    energy_dissip_damp_all = 0.0;
    double    energy_dissip_all      = 0.0;

    double    energy_thermal_slid_pp_all = 0.0;
    double    energy_thermal_slid_pw_all = 0.0;
    double    energy_thermal_roll_pp_all = 0.0;
    double    energy_thermal_roll_pw_all = 0.0;
    double    energy_thermal_damp_pp_all = 0.0;
    double    energy_thermal_damp_pw_all = 0.0;
    double    energy_thermal_all         = 0.0;

    double    heat_flux_values_conducdir_all   = 0.0;
    double    heat_flux_values_conducindir_all = 0.0;
    double    heat_flux_values_rad_all         = 0.0;
    double    heat_flux_values_gen_all         = 0.0;
    double    heat_flux_values_conv_all        = 0.0;
    double    heat_flux_values_prescsurf_all   = 0.0;
    double    heat_flux_values_prescvol_all    = 0.0;
    double    heat_flux_values_all             = 0.0;

    double    heat_flux_contrib_conducdir   = 0.0;
    double    heat_flux_contrib_conducindir = 0.0;
    double    heat_flux_contrib_rad         = 0.0;
    double    heat_flux_contrib_gen         = 0.0;
    double    heat_flux_contrib_conv        = 0.0;
    double    heat_flux_contrib_prescsurf   = 0.0;
    double    heat_flux_contrib_prescvol    = 0.0;

    double    heat_gen_values_slid_pp_all = 0.0;
    double    heat_gen_values_slid_pw_all = 0.0;
    double    heat_gen_values_roll_pp_all = 0.0;
    double    heat_gen_values_roll_pw_all = 0.0;
    double    heat_gen_values_damp_pp_all = 0.0;
    double    heat_gen_values_damp_pw_all = 0.0;
    double    heat_gen_values_all         = 0.0;

    double    heat_gen_contrib_slid_pp = 0.0;
    double    heat_gen_contrib_slid_pw = 0.0;
    double    heat_gen_contrib_roll_pp = 0.0;
    double    heat_gen_contrib_roll_pw = 0.0;
    double    heat_gen_contrib_damp_pp = 0.0;
    double    heat_gen_contrib_damp_pw = 0.0;

    ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();
    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < num_particles; i++) {
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*(it+i));

      // Particle results
      const double vol  = particle.CalculateVolume();
      const double temp = particle.GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);

      double energy_mech_gravity     = 0.0;
      double energy_mech_elastic     = 0.0;
      double energy_mech_translation = 0.0;
      double energy_mech_rotation    = 0.0;

      double energy_dissip_slid = 0.0;
      double energy_dissip_roll = 0.0;
      double energy_dissip_damp = 0.0;

      particle.Calculate(PARTICLE_GRAVITATIONAL_ENERGY,                energy_mech_gravity,     r_process_info);
      particle.Calculate(PARTICLE_ELASTIC_ENERGY,                      energy_mech_elastic,     r_process_info);
      particle.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY,      energy_mech_translation, r_process_info);
      particle.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY,         energy_mech_rotation,    r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_FRICTIONAL_ENERGY,         energy_dissip_slid,      r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_ROLLING_RESISTANCE_ENERGY, energy_dissip_roll,      r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_VISCODAMPING_ENERGY,       energy_dissip_damp,      r_process_info);

      const double energy_thermal_slid_pp = std::abs(particle.mGenerationThermalEnergy_slid_particle);
      const double energy_thermal_slid_pw = std::abs(particle.mGenerationThermalEnergy_slid_wall);
      const double energy_thermal_roll_pp = std::abs(particle.mGenerationThermalEnergy_roll_particle);
      const double energy_thermal_roll_pw = std::abs(particle.mGenerationThermalEnergy_roll_wall);
      const double energy_thermal_damp_pp = std::abs(particle.mGenerationThermalEnergy_damp_particle);
      const double energy_thermal_damp_pw = std::abs(particle.mGenerationThermalEnergy_damp_wall);

      const double heat_flux_values_conducdir   = std::abs(particle.mConductionDirectHeatFlux);
      const double heat_flux_values_conducindir = std::abs(particle.mConductionIndirectHeatFlux);
      const double heat_flux_values_rad         = std::abs(particle.mRadiationHeatFlux);
      const double heat_flux_values_gen         = std::abs(particle.mGenerationHeatFlux);
      const double heat_flux_values_conv        = std::abs(particle.mConvectionHeatFlux);
      const double heat_flux_values_prescsurf   = std::abs(particle.mPrescribedHeatFluxSurface);
      const double heat_flux_values_prescvol    = std::abs(particle.mPrescribedHeatFluxVolume);
      const double heat_flux_values             = std::abs(particle.mTotalHeatFlux);

      const double heat_gen_values_slid_pp = std::abs(particle.mGenerationHeatFlux_slid_particle);
      const double heat_gen_values_slid_pw = std::abs(particle.mGenerationHeatFlux_slid_wall);
      const double heat_gen_values_roll_pp = std::abs(particle.mGenerationHeatFlux_roll_particle);
      const double heat_gen_values_roll_pw = std::abs(particle.mGenerationHeatFlux_roll_wall);
      const double heat_gen_values_damp_pp = std::abs(particle.mGenerationHeatFlux_damp_particle);
      const double heat_gen_values_damp_pw = std::abs(particle.mGenerationHeatFlux_damp_wall);
      const double heat_gen_values         = std::abs(particle.mGenerationHeatFlux);

      // Accumulate particle results
      #pragma omp critical
      {
        vol_total += vol;
        if (temp < particle_temp_min) particle_temp_min = temp;
        if (temp > particle_temp_max) particle_temp_max = temp;
        particle_temp_avg     += temp;
        particle_temp_avg_vol += temp * vol;
        particle_temp_dev     += temp * temp;

        energy_mech_gravity_all     += energy_mech_gravity;
        energy_mech_elastic_all     += energy_mech_elastic;
        energy_mech_translation_all += energy_mech_translation;
        energy_mech_rotation_all    += energy_mech_rotation;
        energy_mech_all             += energy_mech_gravity + energy_mech_elastic + energy_mech_translation + energy_mech_rotation;
        
        energy_dissip_slid_all += energy_dissip_slid;
        energy_dissip_roll_all += energy_dissip_roll;
        energy_dissip_damp_all += energy_dissip_damp;
        energy_dissip_all      += energy_dissip_slid + energy_dissip_roll + energy_dissip_damp;

        energy_thermal_slid_pp_all += energy_thermal_slid_pp;
        energy_thermal_slid_pw_all += energy_thermal_slid_pw;
        energy_thermal_roll_pp_all += energy_thermal_roll_pp;
        energy_thermal_roll_pw_all += energy_thermal_roll_pw;
        energy_thermal_damp_pp_all += energy_thermal_damp_pp;
        energy_thermal_damp_pw_all += energy_thermal_damp_pw;
        energy_thermal_all         += energy_thermal_slid_pp + energy_thermal_slid_pw + energy_thermal_roll_pp + energy_thermal_roll_pw + energy_thermal_damp_pp + energy_thermal_damp_pw;

        heat_flux_values_conducdir_all   += heat_flux_values_conducdir;
        heat_flux_values_conducindir_all += heat_flux_values_conducindir;
        heat_flux_values_rad_all         += heat_flux_values_rad;
        heat_flux_values_gen_all         += heat_flux_values_gen;
        heat_flux_values_conv_all        += heat_flux_values_conv;
        heat_flux_values_prescsurf_all   += heat_flux_values_prescsurf;
        heat_flux_values_prescvol_all    += heat_flux_values_prescvol;
        heat_flux_values_all             += heat_flux_values_conducdir + heat_flux_values_conducindir + heat_flux_values_rad + heat_flux_values_gen + heat_flux_values_conv + heat_flux_values_prescsurf + heat_flux_values_prescvol;

        if (heat_flux_values != 0.0) {
          num_particles_flux_contrib++;
          heat_flux_contrib_conducdir   += heat_flux_values_conducdir   / heat_flux_values;
          heat_flux_contrib_conducindir += heat_flux_values_conducindir / heat_flux_values;
          heat_flux_contrib_rad         += heat_flux_values_rad         / heat_flux_values;
          heat_flux_contrib_gen         += heat_flux_values_gen         / heat_flux_values;
          heat_flux_contrib_conv        += heat_flux_values_conv        / heat_flux_values;
          heat_flux_contrib_prescsurf   += heat_flux_values_prescsurf   / heat_flux_values;
          heat_flux_contrib_prescvol    += heat_flux_values_prescvol    / heat_flux_values;
        }

        heat_gen_values_slid_pp_all += heat_gen_values_slid_pp;
        heat_gen_values_slid_pw_all += heat_gen_values_slid_pw;
        heat_gen_values_roll_pp_all += heat_gen_values_roll_pp;
        heat_gen_values_roll_pw_all += heat_gen_values_roll_pw;
        heat_gen_values_damp_pp_all += heat_gen_values_damp_pp;
        heat_gen_values_damp_pw_all += heat_gen_values_damp_pw;
        heat_gen_values_all         += heat_gen_values_slid_pp + heat_gen_values_slid_pw + heat_gen_values_roll_pp + heat_gen_values_roll_pw + heat_gen_values_damp_pp + heat_gen_values_damp_pw;

        if (heat_gen_values != 0.0) {
          num_particles_gen_contrib++;
          heat_gen_contrib_slid_pp += heat_gen_values_slid_pp / heat_gen_values;
          heat_gen_contrib_slid_pw += heat_gen_values_slid_pw / heat_gen_values;
          heat_gen_contrib_roll_pp += heat_gen_values_roll_pp / heat_gen_values;
          heat_gen_contrib_roll_pw += heat_gen_values_roll_pw / heat_gen_values;
          heat_gen_contrib_damp_pp += heat_gen_values_damp_pp / heat_gen_values;
          heat_gen_contrib_damp_pw += heat_gen_values_damp_pw / heat_gen_values;
        }
      }
    }

    // Compute temperature results (avg, dev)
    particle_temp_avg     /= num_particles;
    particle_temp_avg_vol /= vol_total;
    particle_temp_dev      = sqrt(std::max(0.0, particle_temp_dev / num_particles - particle_temp_avg * particle_temp_avg));

    // Compute average of relative contribution of each heat transfer mechanism
    if (num_particles_flux_contrib > 0) {
      heat_flux_contrib_conducdir   /= num_particles_flux_contrib;
      heat_flux_contrib_conducindir /= num_particles_flux_contrib;
      heat_flux_contrib_rad         /= num_particles_flux_contrib;
      heat_flux_contrib_gen         /= num_particles_flux_contrib;
      heat_flux_contrib_conv        /= num_particles_flux_contrib;
      heat_flux_contrib_prescsurf   /= num_particles_flux_contrib;
      heat_flux_contrib_prescvol    /= num_particles_flux_contrib;
    }

    // Compute average of relative contribution of each heat generation mechanism
    if (num_particles_gen_contrib > 0) {
      heat_gen_contrib_slid_pp /= num_particles_gen_contrib;
      heat_gen_contrib_slid_pw /= num_particles_gen_contrib;
      heat_gen_contrib_roll_pp /= num_particles_gen_contrib;
      heat_gen_contrib_roll_pw /= num_particles_gen_contrib;
      heat_gen_contrib_damp_pp /= num_particles_gen_contrib;
      heat_gen_contrib_damp_pw /= num_particles_gen_contrib;
    }

    // Write results to files
    if (mFile_ParticleTempMin.is_open())
      mFile_ParticleTempMin << time_step << " "
                            << time      << " "
                            << particle_temp_min
                            << std::endl;

    if (mFile_ParticleTempMax.is_open())
      mFile_ParticleTempMax << time_step << " "
                            << time      << " "
                            << particle_temp_max
                            << std::endl;

    if (mFile_ParticleTempAvg.is_open())
      mFile_ParticleTempAvg << time_step << " "
                            << time      << " "
                            << particle_temp_avg
                            << std::endl;

    if (mFile_ParticleTempAvgVol.is_open())
      mFile_ParticleTempAvgVol << time_step << " "
                               << time      << " "
                               << particle_temp_avg_vol
                               << std::endl;

    if (mFile_ParticleTempDev.is_open())
      mFile_ParticleTempDev << time_step << " "
                            << time      << " "
                            << particle_temp_dev
                            << std::endl;

    if (mFile_EnergyMechanical.is_open())
      mFile_EnergyMechanical << time_step                   << " "
                             << time                        << " "
                             << energy_mech_gravity_all     << " "
                             << energy_mech_elastic_all     << " "
                             << energy_mech_translation_all << " "
                             << energy_mech_rotation_all    << " "
                             << energy_mech_all
                             << std::endl;

    if (mFile_EnergyDissipated.is_open())
      mFile_EnergyDissipated << time_step              << " "
                             << time                   << " "
                             << energy_dissip_slid_all << " "
                             << energy_dissip_roll_all << " "
                             << energy_dissip_damp_all << " "
                             << energy_dissip_all
                             << std::endl;

    if (mFile_EnergyThermal.is_open())
      mFile_EnergyThermal << time_step                  << " "
                          << time                       << " "
                          << energy_thermal_slid_pp_all << " "
                          << energy_thermal_slid_pw_all << " "
                          << energy_thermal_roll_pp_all << " "
                          << energy_thermal_roll_pw_all << " "
                          << energy_thermal_damp_pp_all << " "
                          << energy_thermal_damp_pw_all << " "
                          << energy_thermal_all        
                          << std::endl;

    if (mFile_HeatFluxContributions.is_open())
      mFile_HeatFluxContributions << time_step                     << " "
                                  << time                          << " "
                                  << heat_flux_contrib_conducdir   << " "
                                  << heat_flux_contrib_conducindir << " "
                                  << heat_flux_contrib_rad         << " "
                                  << heat_flux_contrib_gen         << " "
                                  << heat_flux_contrib_conv        << " "
                                  << heat_flux_contrib_prescsurf   << " "
                                  << heat_flux_contrib_prescvol
                                  << std::endl;

    if (mFile_HeatGenValues.is_open())
      mFile_HeatGenValues << time_step                   << " "
                          << time                        << " "
                          << heat_gen_values_slid_pp_all << " "
                          << heat_gen_values_slid_pw_all << " "
                          << heat_gen_values_roll_pp_all << " "
                          << heat_gen_values_roll_pw_all << " "
                          << heat_gen_values_damp_pp_all << " "
                          << heat_gen_values_damp_pw_all << " "
                          << heat_gen_values_all
                          << std::endl;

    if (mFile_HeatGenContributions.is_open())
      mFile_HeatGenContributions << time_step                << " "
                                 << time                     << " "
                                 << heat_gen_contrib_slid_pp << " "
                                 << heat_gen_contrib_slid_pw << " "
                                 << heat_gen_contrib_roll_pp << " "
                                 << heat_gen_contrib_roll_pw << " "
                                 << heat_gen_contrib_damp_pp << " "
                                 << heat_gen_contrib_damp_pw
                                 << std::endl;

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::CloseFiles(void) {
    if (mFile_ParticleTempMin.is_open())       mFile_ParticleTempMin.close();
    if (mFile_ParticleTempMax.is_open())       mFile_ParticleTempMax.close();
    if (mFile_ParticleTempAvg.is_open())       mFile_ParticleTempAvg.close();
    if (mFile_ParticleTempAvgVol.is_open())    mFile_ParticleTempAvgVol.close();
    if (mFile_ParticleTempDev.is_open())       mFile_ParticleTempDev.close();
    if (mFile_EnergyMechanical.is_open())      mFile_EnergyMechanical.close();
    if (mFile_EnergyDissipated.is_open())      mFile_EnergyDissipated.close();
    if (mFile_EnergyThermal.is_open())         mFile_EnergyThermal.close();
    if (mFile_HeatFluxContributions.is_open()) mFile_HeatFluxContributions.close();
    if (mFile_HeatGenValues.is_open())         mFile_HeatGenValues.close();
    if (mFile_HeatGenContributions.is_open())  mFile_HeatGenContributions.close();
    if (mFile_ParticleTempAll.is_open())               mFile_ParticleTempAll.close();
  }
  
} // namespace Kratos
