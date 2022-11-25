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
  GraphUtilities::GraphUtilities()
  {
    mGraph_ParticleTempMin               = false;
    mGraph_ParticleTempMax               = false;
    mGraph_ParticleTempAvg               = false;
    mGraph_ParticleTempDev               = false;
    mGraph_ModelTempAvg                  = false;
    mGraph_ParticleHeatFluxContributions = false;
    mGraph_ParticleHeatGenContributions  = false;
    mGraph_ParticleEnergyContributions   = false;
  }

  GraphUtilities::~GraphUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteInitialize(bool ParticleTempMin,
                                         bool ParticleTempMax,
                                         bool ParticleTempAvg,
                                         bool ParticleTempDev,
                                         bool ModelTempAvg,
                                         bool ParticleHeatFluxContributions,
                                         bool ParticleHeatGenContributions,
                                         bool ParticleEnergyContributions)
  {
    KRATOS_TRY

    // Set member flags
    mGraph_ParticleTempMin               = ParticleTempMin;
    mGraph_ParticleTempMax               = ParticleTempMax;
    mGraph_ParticleTempAvg               = ParticleTempAvg;
    mGraph_ParticleTempDev               = ParticleTempDev;
    mGraph_ModelTempAvg                  = ModelTempAvg;
    mGraph_ParticleHeatFluxContributions = ParticleHeatFluxContributions;
    mGraph_ParticleHeatGenContributions  = ParticleHeatGenContributions;
    mGraph_ParticleEnergyContributions   = ParticleEnergyContributions;

    // Open files
    if (mGraph_ParticleTempMin) {
      mFile_ParticleTempMin.open("graph_particle_temp_min.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempMin) << "Could not open graph file for minimum particle temperature!" << std::endl;
      mFile_ParticleTempMin << "1 - TIME STEP | ";
      mFile_ParticleTempMin << "2 - TIME | ";
      mFile_ParticleTempMin << "3 - MIN PARTICLE TEMPERATURE";
      mFile_ParticleTempMin << std::endl;
    }
    if (mGraph_ParticleTempMax) {
      mFile_ParticleTempMax.open("graph_particle_temp_max.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempMax) << "Could not open graph file for maximum particle temperature!" << std::endl;
      mFile_ParticleTempMax << "1 - TIME STEP | ";
      mFile_ParticleTempMax << "2 - TIME | ";
      mFile_ParticleTempMax << "3 - MAX PARTICLE TEMPERATURE";
      mFile_ParticleTempMax << std::endl;
    }
    if (mGraph_ParticleTempAvg) {
      mFile_ParticleTempAvg.open("graph_particle_temp_avg.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempAvg) << "Could not open graph file for average particle temperature!" << std::endl;
      mFile_ParticleTempAvg << "1 - TIME STEP | ";
      mFile_ParticleTempAvg << "2 - TIME | ";
      mFile_ParticleTempAvg << "3 - AVERAGE PARTICLE TEMPERATURE";
      mFile_ParticleTempAvg << std::endl;
    }
    if (mGraph_ParticleTempDev) {
      mFile_ParticleTempDev.open("graph_particle_temp_dev.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempDev) << "Could not open graph file for deviation of particle temperature!" << std::endl;
      mFile_ParticleTempDev << "1 - TIME STEP | ";
      mFile_ParticleTempDev << "2 - TIME | ";
      mFile_ParticleTempDev << "3 - PARTICLE TEMPERATURE STANDARD DEVIATION";
      mFile_ParticleTempDev << std::endl;
    }
    if (mGraph_ModelTempAvg) {
      mFile_ModelTempAvg.open("graph_model_temp_avg.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ModelTempAvg) << "Could not open graph file for average model temperature!" << std::endl;
      mFile_ModelTempAvg << "1 - TIME STEP | ";
      mFile_ModelTempAvg << "2 - TIME | ";
      mFile_ModelTempAvg << "3 - AVERAGE MODEL TEMPERATURE";
      mFile_ModelTempAvg << std::endl;
    }
    if (mGraph_ParticleHeatFluxContributions) {
      mFile_ParticleHeatFluxContributions.open("graph_flux_contributions.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleHeatFluxContributions) << "Could not open graph file for heat flux contributions!" << std::endl;
      mFile_ParticleHeatFluxContributions << "1 - TIME STEP | ";
      mFile_ParticleHeatFluxContributions << "2 - TIME | ";
      mFile_ParticleHeatFluxContributions << "3 - DIRECT CONDUCTION | ";
      mFile_ParticleHeatFluxContributions << "4 - INDIRECT CONDUCTION | ";
      mFile_ParticleHeatFluxContributions << "5 - RADIATION | ";
      mFile_ParticleHeatFluxContributions << "6 - GENERATION | ";
      mFile_ParticleHeatFluxContributions << "7 - CONVECTION | ";
      mFile_ParticleHeatFluxContributions << "8 - SURFACE PRESCRIBED | ";
      mFile_ParticleHeatFluxContributions << "9 - VOLUME PRESCRIBED";
      mFile_ParticleHeatFluxContributions << std::endl;
    }
    if (mGraph_ParticleHeatGenContributions) {
      mFile_ParticleHeatGenContributions.open("graph_generation_contributions.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleHeatGenContributions) << "Could not open graph file for heat generation contributions!" << std::endl;
      mFile_ParticleHeatGenContributions << "1 - TIME STEP | ";
      mFile_ParticleHeatGenContributions << "2 - TIME | ";
      mFile_ParticleHeatGenContributions << "3 - SLIDING FRICTION PARTICLE-PARTICLE | ";
      mFile_ParticleHeatGenContributions << "4 - SLIDING FRICTION PARTICLE-WALL | ";
      mFile_ParticleHeatGenContributions << "5 - ROLLING FRICTION PARTICLE-PARTICLE | ";
      mFile_ParticleHeatGenContributions << "6 - ROLLING FRICTION PARTICLE-WALL | ";
      mFile_ParticleHeatGenContributions << "7 - DAMPING FORCE PARTICLE-PARTICLE | ";
      mFile_ParticleHeatGenContributions << "8 - DAMPING FORCE PARTICLE-WALL";
      mFile_ParticleHeatGenContributions << std::endl;
    }
    if (mGraph_ParticleEnergyContributions) {
      mFile_ParticleEnergyContributions.open("graph_energy_contributions.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleEnergyContributions) << "Could not open graph file for energy contributions!" << std::endl;
      mFile_ParticleEnergyContributions << "1 - TIME STEP | ";
      mFile_ParticleEnergyContributions << "2 - TIME | ";
      mFile_ParticleEnergyContributions << "3 - GRAVITATIONAL | ";
      mFile_ParticleEnergyContributions << "4 - ELASTIC | ";
      mFile_ParticleEnergyContributions << "5 - KINETIC TRANSLATION | ";
      mFile_ParticleEnergyContributions << "6 - KINETIC ROTATION | ";
      mFile_ParticleEnergyContributions << "7 - TOTAL ENERGY | ";
      mFile_ParticleEnergyContributions << "8 - ACCUMULATED DISSIP FRICTION | ";
      mFile_ParticleEnergyContributions << "9 - ACCUMULATED DISSIP ROLL | ";
      mFile_ParticleEnergyContributions << "10 - ACCUMULATED DISSIP DAMPING | ";
      mFile_ParticleEnergyContributions << "11 - TOTAL ACCUMULATED DISSIP | ";
      mFile_ParticleEnergyContributions << "12 - TOTAL ENERGY + ACCUMULATED DISSIP";
      mFile_ParticleEnergyContributions << std::endl;
    }

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteFinalizeSolutionStep(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    if (!r_process_info[IS_TIME_TO_PRINT])
      return;

    // Initialize results
    const int num_of_particles                    =  rModelPart.NumberOfElements();
    int       num_ratio_particles_flux            =  0;
    int       num_ratio_particles_gen             =  0;
    double    total_vol                           =  0.0;

    double    particle_temp_min                   =  DBL_MAX;
    double    particle_temp_max                   = -DBL_MAX;
    double    particle_temp_avg                   =  0.0;
    double    particle_temp_dev                   =  0.0;
    double    model_temp_avg                      =  0.0;

    double    particle_flux_conducdir_ratio_avg   =  0.0;
    double    particle_flux_conducindir_ratio_avg =  0.0;
    double    particle_flux_rad_ratio_avg         =  0.0;
    double    particle_flux_gen_ratio_avg         =  0.0;
    double    particle_flux_conv_ratio_avg        =  0.0;
    double    particle_flux_prescsurf_ratio_avg   =  0.0;
    double    particle_flux_prescvol_ratio_avg    =  0.0;

    double    particle_gen_slid_pp_ratio_avg      =  0.0;
    double    particle_gen_slid_pw_ratio_avg      =  0.0;
    double    particle_gen_roll_pp_ratio_avg      =  0.0;
    double    particle_gen_roll_pw_ratio_avg      =  0.0;
    double    particle_gen_damp_pp_ratio_avg      =  0.0;
    double    particle_gen_damp_pw_ratio_avg      =  0.0;

    double    total_energy_potential_gravity      =  0.0;
    double    total_energy_potential_elastic      =  0.0;
    double    total_energy_kinetic_translation    =  0.0;
    double    total_energy_kinetic_rotation       =  0.0;
    double    total_energy                        =  0.0;

    double    total_accum_energy_dissip_slid      =  0.0;
    double    total_accum_energy_dissip_roll      =  0.0;
    double    total_accum_energy_dissip_damp      =  0.0;
    double    total_accum_energy_dissip           =  0.0;

    ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();
    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < num_of_particles; i++) {
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*(it+i));

      // Accumulate values
      const double vol                  = particle.CalculateVolume();
      const double temp                 = particle.GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
      double energy_potential_gravity   = 0.0;
      double energy_potential_elastic   = 0.0;
      double energy_kinetic_translation = 0.0;
      double energy_kinetic_rotation    = 0.0;
      double accum_energy_dissip_slid   = 0.0;
      double accum_energy_dissip_roll   = 0.0;
      double accum_energy_dissip_damp   = 0.0;
      particle.Calculate(PARTICLE_GRAVITATIONAL_ENERGY,                energy_potential_gravity,   r_process_info);
      particle.Calculate(PARTICLE_ELASTIC_ENERGY,                      energy_potential_elastic,   r_process_info);
      particle.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY,      energy_kinetic_translation, r_process_info);
      particle.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY,         energy_kinetic_rotation,    r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_FRICTIONAL_ENERGY,         accum_energy_dissip_slid,   r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_ROLLING_RESISTANCE_ENERGY, accum_energy_dissip_roll,   r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_VISCODAMPING_ENERGY,       accum_energy_dissip_damp,   r_process_info);

      #pragma omp critical
      {
        total_vol += vol;
        if (temp < particle_temp_min) particle_temp_min = temp;
        if (temp > particle_temp_max) particle_temp_max = temp;
        particle_temp_avg                += temp;
        particle_temp_dev                += temp * temp;
        model_temp_avg                   += temp * vol;

        total_energy_potential_gravity   += energy_potential_gravity;
        total_energy_potential_elastic   += energy_potential_elastic;
        total_energy_kinetic_translation += energy_kinetic_translation;
        total_energy_kinetic_rotation    += energy_kinetic_rotation;
        total_energy                     += energy_potential_gravity + energy_potential_elastic + energy_kinetic_translation + energy_kinetic_rotation;
        
        total_accum_energy_dissip_slid   += accum_energy_dissip_slid;
        total_accum_energy_dissip_roll   += accum_energy_dissip_roll;
        total_accum_energy_dissip_damp   += accum_energy_dissip_damp;
        total_accum_energy_dissip        += accum_energy_dissip_slid + accum_energy_dissip_roll + accum_energy_dissip_damp;
      }

      if (mGraph_ParticleHeatFluxContributions) {
        // Get absolute value of particle heat transfer mechanisms
        const double flux_conducdir   = fabs(particle.mConductionDirectHeatFlux);
        const double flux_conducindir = fabs(particle.mConductionIndirectHeatFlux);
        const double flux_rad         = fabs(particle.mRadiationHeatFlux);
        const double flux_gen         = fabs(particle.mGenerationHeatFlux);
        const double flux_conv        = fabs(particle.mConvectionHeatFlux);
        const double flux_prescsurf   = fabs(particle.mPrescribedHeatFluxSurface);
        const double flux_prescvol    = fabs(particle.mPrescribedHeatFluxVolume);
        const double flux_total       = flux_conducdir + flux_conducindir + flux_rad + flux_gen + flux_conv + flux_prescsurf + flux_prescvol;

        // Compute relative contribution of each heat transfer mechanism for current particle
        if (flux_total != 0.0) {
          #pragma omp critical
          {
            num_ratio_particles_flux++;
            particle_flux_conducdir_ratio_avg   += flux_conducdir   / flux_total;
            particle_flux_conducindir_ratio_avg += flux_conducindir / flux_total;
            particle_flux_rad_ratio_avg         += flux_rad         / flux_total;
            particle_flux_gen_ratio_avg         += flux_gen         / flux_total;
            particle_flux_conv_ratio_avg        += flux_conv        / flux_total;
            particle_flux_prescsurf_ratio_avg   += flux_prescsurf   / flux_total;
            particle_flux_prescvol_ratio_avg    += flux_prescvol    / flux_total;
          }
        }
      }

      if (mGraph_ParticleHeatGenContributions) {
        // Get absolute value of particle heat generation mechanisms
        const double gen_slid_pp = fabs(particle.mGenerationHeatFlux_slid_particle);
        const double gen_slid_pw = fabs(particle.mGenerationHeatFlux_slid_wall);
        const double gen_roll_pp = fabs(particle.mGenerationHeatFlux_roll_particle);
        const double gen_roll_pw = fabs(particle.mGenerationHeatFlux_roll_wall);
        const double gen_damp_pp = fabs(particle.mGenerationHeatFlux_damp_particle);
        const double gen_damp_pw = fabs(particle.mGenerationHeatFlux_damp_wall);
        const double gen_total   = fabs(particle.mGenerationHeatFlux);

        // Compute relative contribution of each heat generation mechanism for current particle
        if (gen_total != 0.0) {
          #pragma omp critical
          {
            num_ratio_particles_gen++;
            particle_gen_slid_pp_ratio_avg += gen_slid_pp / gen_total;
            particle_gen_slid_pw_ratio_avg += gen_slid_pw / gen_total;
            particle_gen_roll_pp_ratio_avg += gen_roll_pp / gen_total;
            particle_gen_roll_pw_ratio_avg += gen_roll_pw / gen_total;
            particle_gen_damp_pp_ratio_avg += gen_damp_pp / gen_total;
            particle_gen_damp_pw_ratio_avg += gen_damp_pw / gen_total;
          }
        }
      }
    }

    // Compute temperature results (avg, dev)
    particle_temp_avg /= num_of_particles;
    model_temp_avg    /= total_vol;
    particle_temp_dev  = sqrt(std::max(0.0, particle_temp_dev / num_of_particles - particle_temp_avg * particle_temp_avg));

    // Compute average of relative contribution of each heat transfer mechanism
    if (mGraph_ParticleHeatFluxContributions && num_ratio_particles_flux > 0) {
      particle_flux_conducdir_ratio_avg   /= num_ratio_particles_flux;
      particle_flux_conducindir_ratio_avg /= num_ratio_particles_flux;
      particle_flux_rad_ratio_avg         /= num_ratio_particles_flux;
      particle_flux_gen_ratio_avg         /= num_ratio_particles_flux;
      particle_flux_conv_ratio_avg        /= num_ratio_particles_flux;
      particle_flux_prescsurf_ratio_avg   /= num_ratio_particles_flux;
      particle_flux_prescvol_ratio_avg    /= num_ratio_particles_flux;
    }

    // Compute average of relative contribution of each heat generation mechanism
    if (mGraph_ParticleHeatGenContributions && num_ratio_particles_gen > 0) {
      particle_gen_slid_pp_ratio_avg /= num_ratio_particles_gen;
      particle_gen_slid_pw_ratio_avg /= num_ratio_particles_gen;
      particle_gen_roll_pp_ratio_avg /= num_ratio_particles_gen;
      particle_gen_roll_pw_ratio_avg /= num_ratio_particles_gen;
      particle_gen_damp_pp_ratio_avg /= num_ratio_particles_gen;
      particle_gen_damp_pw_ratio_avg /= num_ratio_particles_gen;
    }

    // Write results to files
    const int    time_step = r_process_info[TIME_STEPS];
    const double time      = r_process_info[TIME];

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

    if (mFile_ParticleTempDev.is_open())
      mFile_ParticleTempDev << time_step << " "
                            << time      << " "
                            << particle_temp_dev
                            << std::endl;

    if (mFile_ModelTempAvg.is_open())
      mFile_ModelTempAvg << time_step << " "
                         << time      << " "
                         << model_temp_avg
                         << std::endl;

    if (mFile_ParticleHeatFluxContributions.is_open())
      mFile_ParticleHeatFluxContributions << time_step                           << " "
                                          << time                                << " "
                                          << particle_flux_conducdir_ratio_avg   << " "
                                          << particle_flux_conducindir_ratio_avg << " "
                                          << particle_flux_rad_ratio_avg         << " "
                                          << particle_flux_gen_ratio_avg         << " "
                                          << particle_flux_conv_ratio_avg        << " "
                                          << particle_flux_prescsurf_ratio_avg   << " "
                                          << particle_flux_prescvol_ratio_avg
                                          << std::endl;

    if (mFile_ParticleHeatGenContributions.is_open())
      mFile_ParticleHeatGenContributions << time_step                      << " "
                                         << time                           << " "
                                         << particle_gen_slid_pp_ratio_avg << " "
                                         << particle_gen_slid_pw_ratio_avg << " "
                                         << particle_gen_roll_pp_ratio_avg << " "
                                         << particle_gen_roll_pw_ratio_avg << " "
                                         << particle_gen_damp_pp_ratio_avg << " "
                                         << particle_gen_damp_pw_ratio_avg
                                         << std::endl;

    if (mFile_ParticleEnergyContributions.is_open())
      mFile_ParticleEnergyContributions << time_step                        << " "
                                        << time                             << " "
                                        << total_energy_potential_gravity   << " "
                                        << total_energy_potential_elastic   << " "
                                        << total_energy_kinetic_translation << " "
                                        << total_energy_kinetic_rotation    << " "
                                        << total_energy                     << " "
                                        << total_accum_energy_dissip_slid   << " "
                                        << total_accum_energy_dissip_roll   << " "
                                        << total_accum_energy_dissip_damp   << " "
                                        << total_accum_energy_dissip        << " "
                                        << total_energy+total_accum_energy_dissip
                                        << std::endl;

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteFinalize(void)
  {
    KRATOS_TRY

    // Close files
    if (mFile_ParticleTempMin.is_open())               mFile_ParticleTempMin.close();
    if (mFile_ParticleTempMax.is_open())               mFile_ParticleTempMax.close();
    if (mFile_ParticleTempAvg.is_open())               mFile_ParticleTempAvg.close();
    if (mFile_ParticleTempDev.is_open())               mFile_ParticleTempDev.close();
    if (mFile_ModelTempAvg.is_open())                  mFile_ModelTempAvg.close();
    if (mFile_ParticleHeatFluxContributions.is_open()) mFile_ParticleHeatFluxContributions.close();
    if (mFile_ParticleHeatGenContributions.is_open())  mFile_ParticleHeatGenContributions.close();
    if (mFile_ParticleEnergyContributions.is_open())   mFile_ParticleEnergyContributions.close();

    KRATOS_CATCH("")
  }

} // namespace Kratos
