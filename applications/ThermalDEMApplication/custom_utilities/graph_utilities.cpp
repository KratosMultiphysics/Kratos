//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
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
  }

  GraphUtilities::~GraphUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteInitialize(bool ParticleTempMin,
                                         bool ParticleTempMax,
                                         bool ParticleTempAvg,
                                         bool ParticleTempDev,
                                         bool ModelTempAvg,
                                         bool ParticleHeatFluxContributions)
  {
    KRATOS_TRY

    // Set member flags
    mGraph_ParticleTempMin               = ParticleTempMin;
    mGraph_ParticleTempMax               = ParticleTempMax;
    mGraph_ParticleTempAvg               = ParticleTempAvg;
    mGraph_ParticleTempDev               = ParticleTempDev;
    mGraph_ModelTempAvg                  = ModelTempAvg;
    mGraph_ParticleHeatFluxContributions = ParticleHeatFluxContributions;

    // Open files
    if (mGraph_ParticleTempMin) {
      mFile_ParticleTempMin.open("graph_particle_temp_min.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempMin) << "Could not open graph file for minimum particle temperature!" << std::endl;
      mFile_ParticleTempMin << "TIME STEP | TIME | MIN PARTICLE TEMPERATURE" << std::endl;
    }
    if (mGraph_ParticleTempMax) {
      mFile_ParticleTempMax.open("graph_particle_temp_max.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempMax) << "Could not open graph file for maximum particle temperature!" << std::endl;
      mFile_ParticleTempMax << "TIME STEP | TIME | MAX PARTICLE TEMPERATURE" << std::endl;
    }
    if (mGraph_ParticleTempAvg) {
      mFile_ParticleTempAvg.open("graph_particle_temp_avg.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempAvg) << "Could not open graph file for average particle temperature!" << std::endl;
      mFile_ParticleTempAvg << "TIME STEP | TIME | AVERAGE PARTICLE TEMPERATURE" << std::endl;
    }
    if (mGraph_ParticleTempDev) {
      mFile_ParticleTempDev.open("graph_particle_temp_dev.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleTempDev) << "Could not open graph file for deviation of particle temperature!" << std::endl;
      mFile_ParticleTempDev << "TIME STEP | TIME | PARTICLE TEMPERATURE STANDARD DEVIATION" << std::endl;
    }
    if (mGraph_ModelTempAvg) {
      mFile_ModelTempAvg.open("graph_model_temp_avg.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ModelTempAvg) << "Could not open graph file for average model temperature!" << std::endl;
      mFile_ModelTempAvg << "TIME STEP | TIME | AVERAGE MODEL TEMPERATURE" << std::endl;
    }
    if (mGraph_ParticleHeatFluxContributions) {
      mFile_ParticleHeatFluxContributions.open("graph_flux_contributions.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ParticleHeatFluxContributions) << "Could not open graph file for heat flux contributions!" << std::endl;
      mFile_ParticleHeatFluxContributions << "TIME STEP | TIME | DIRECT CONDUCTION | INDIRECT CONDUCTION | RADIATION | FRICTION | CONVECTION | SURFACE PRESCRIBED | VOLUME PRESCRIBED" << std::endl;
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
    int       num_ratio_particles                 =  0;
    double    total_vol                           =  0.0;
    double    particle_temp_min                   =  DBL_MAX;
    double    particle_temp_max                   = -DBL_MAX;
    double    particle_temp_avg                   =  0.0;
    double    particle_temp_dev                   =  0.0;
    double    model_temp_avg                      =  0.0;
    double    particle_flux_conducdir_ratio_avg   =  0.0;
    double    particle_flux_conducindir_ratio_avg =  0.0;
    double    particle_flux_rad_ratio_avg         =  0.0;
    double    particle_flux_fric_ratio_avg        =  0.0;
    double    particle_flux_conv_ratio_avg        =  0.0;
    double    particle_flux_prescsurf_ratio_avg   =  0.0;
    double    particle_flux_prescvol_ratio_avg    =  0.0;

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < num_of_particles; i++) {
      ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*it);

      // Accumulate values
      const double vol  = particle.CalculateVolume();
      const double temp = particle.GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
      #pragma omp critical
      {
        total_vol += vol;
        if (temp < particle_temp_min) particle_temp_min = temp;
        if (temp > particle_temp_max) particle_temp_max = temp;
        particle_temp_avg += temp;
        particle_temp_dev += temp * temp;
        model_temp_avg    += temp * vol;
      }

      if (mGraph_ParticleHeatFluxContributions) {
        // Get absolute value of particle heat transfer mechanisms
        const double flux_conducdir   = fabs(particle.mConductionDirectHeatFlux);
        const double flux_conducindir = fabs(particle.mConductionIndirectHeatFlux);
        const double flux_rad         = fabs(particle.mRadiationHeatFlux);
        const double flux_fric        = fabs(particle.mFrictionHeatFlux);
        const double flux_conv        = fabs(particle.mConvectionHeatFlux);
        const double flux_prescsurf   = fabs(particle.mPrescribedHeatFluxSurface);
        const double flux_prescvol    = fabs(particle.mPrescribedHeatFluxVolume);
        const double flux_total       = flux_conducdir + flux_conducindir + flux_rad + flux_fric + flux_conv + flux_prescsurf + flux_prescvol;

        // Compute relative contribution of each heat transfer mechanism for current particle
        if (flux_total != 0.0) {
          #pragma omp critical
          {
            num_ratio_particles++;
            particle_flux_conducdir_ratio_avg   += flux_conducdir   / flux_total;
            particle_flux_conducindir_ratio_avg += flux_conducindir / flux_total;
            particle_flux_rad_ratio_avg         += flux_rad         / flux_total;
            particle_flux_fric_ratio_avg        += flux_fric        / flux_total;
            particle_flux_conv_ratio_avg        += flux_conv        / flux_total;
            particle_flux_prescsurf_ratio_avg   += flux_prescsurf   / flux_total;
            particle_flux_prescvol_ratio_avg    += flux_prescvol    / flux_total;
          }
        }
      }
    }

    // Compute temperature results (avg, dev)
    particle_temp_avg /= num_of_particles;
    model_temp_avg    /= total_vol;
    particle_temp_dev  = sqrt(std::max(0.0, particle_temp_dev / num_of_particles - particle_temp_avg * particle_temp_avg));

    // Compute average of relative contribution of each heat transfer mechanism
    if (mGraph_ParticleHeatFluxContributions && num_ratio_particles > 0) {
      particle_flux_conducdir_ratio_avg   /= num_ratio_particles;
      particle_flux_conducindir_ratio_avg /= num_ratio_particles;
      particle_flux_rad_ratio_avg         /= num_ratio_particles;
      particle_flux_fric_ratio_avg        /= num_ratio_particles;
      particle_flux_conv_ratio_avg        /= num_ratio_particles;
      particle_flux_prescsurf_ratio_avg   /= num_ratio_particles;
      particle_flux_prescvol_ratio_avg    /= num_ratio_particles;
    }

    // Write results to files
    const int    time_step = r_process_info[TIME_STEPS];
    const double time      = r_process_info[TIME];

    if (mFile_ParticleTempMin.is_open())
      mFile_ParticleTempMin << time_step << " " << time << " " << particle_temp_min << std::endl;
    if (mFile_ParticleTempMax.is_open())
      mFile_ParticleTempMax << time_step << " " << time << " " << particle_temp_max << std::endl;
    if (mFile_ParticleTempAvg.is_open())
      mFile_ParticleTempAvg << time_step << " " << time << " " << particle_temp_avg << std::endl;
    if (mFile_ParticleTempDev.is_open())
      mFile_ParticleTempDev << time_step << " " << time << " " << particle_temp_dev << std::endl;
    if (mFile_ModelTempAvg.is_open())
      mFile_ModelTempAvg << time_step << " " << time << " " << model_temp_avg << std::endl;
    if (mFile_ParticleHeatFluxContributions.is_open())
      mFile_ParticleHeatFluxContributions << time_step << " " << time << " " << particle_flux_conducdir_ratio_avg << " " << particle_flux_conducindir_ratio_avg << " " << particle_flux_rad_ratio_avg << " " << particle_flux_fric_ratio_avg << " " << particle_flux_conv_ratio_avg << " " << particle_flux_prescsurf_ratio_avg << " " << particle_flux_prescvol_ratio_avg << std::endl;

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

    KRATOS_CATCH("")
  }

} // namespace Kratos
