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
    mGraph_ParticleTempMin  = false;
    mGraph_ParticleTempMax  = false;
    mGraph_ParticleTempAvg  = false;
    mGraph_ParticleTempDev  = false;
    mGraph_ModelTempAvg     = false;
    mGraph_MechanicalEnergy = false;
    mGraph_ThermalEnergy    = false;
    mGraph_HeatGenValues    = false;
    mGraph_HeatGenContrib   = false;
  }

  GraphUtilities::~GraphUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteInitialize(bool ParticleTempMin,
                                         bool ParticleTempMax,
                                         bool ParticleTempAvg,
                                         bool ParticleTempDev,
                                         bool ModelTempAvg,
                                         bool MechanicalEnergy,
                                         bool ThermalEnergy,
                                         bool HeatGenValues,
                                         bool HeatGenContrib)
  {
    KRATOS_TRY

    // Set member flags
    mGraph_ParticleTempMin  = ParticleTempMin;
    mGraph_ParticleTempMax  = ParticleTempMax;
    mGraph_ParticleTempAvg  = ParticleTempAvg;
    mGraph_ParticleTempDev  = ParticleTempDev;
    mGraph_ModelTempAvg     = ModelTempAvg;
    mGraph_MechanicalEnergy = MechanicalEnergy;
    mGraph_ThermalEnergy    = ThermalEnergy;
    mGraph_HeatGenValues    = HeatGenValues;
    mGraph_HeatGenContrib   = HeatGenContrib;

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
    if (mGraph_MechanicalEnergy) {
      mFile_MechanicalEnergy.open("graph_mechanical_energy.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_MechanicalEnergy) << "Could not open graph file graph_mechanical_energy.txt!" << std::endl;
      mFile_MechanicalEnergy << "1 - TIME STEP | ";
      mFile_MechanicalEnergy << "2 - TIME | ";
      mFile_MechanicalEnergy << "3 - GRAV | ";
      mFile_MechanicalEnergy << "4 - ELAST | ";
      mFile_MechanicalEnergy << "5 - KINET TRANS | ";
      mFile_MechanicalEnergy << "6 - KINET ROT | ";
      mFile_MechanicalEnergy << "7 - TOTAL ENERGY | ";
      mFile_MechanicalEnergy << "8 - DISSIP SLID | ";
      mFile_MechanicalEnergy << "9 - DISSIP ROLL | ";
      mFile_MechanicalEnergy << "10 - DISSIP DAMP | ";
      mFile_MechanicalEnergy << "11 - DISSIP DAMP N | ";
      mFile_MechanicalEnergy << "12 - DISSIP DAMP T | ";
      mFile_MechanicalEnergy << "13 - TOTAL DISSIP | ";
      mFile_MechanicalEnergy << "14 - ENERGY + DISSIP";
      mFile_MechanicalEnergy << std::endl;
    }
    if (mGraph_ThermalEnergy) {
      mFile_ThermalEnergy.open("graph_thermal_energy.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_ThermalEnergy) << "Could not open graph file graph_thermal_energy.txt!" << std::endl;
      mFile_ThermalEnergy << "1 - TIME STEP | ";
      mFile_ThermalEnergy << "2 - TIME | ";
      mFile_ThermalEnergy << "3 - SLID PP | ";
      mFile_ThermalEnergy << "4 - SLID PW | ";
      mFile_ThermalEnergy << "5 - ROLL PP | ";
      mFile_ThermalEnergy << "6 - ROLL PW | ";
      mFile_ThermalEnergy << "7 - DAMP PP | ";
      mFile_ThermalEnergy << "8 - DAMP PW | ";
      mFile_ThermalEnergy << "9 - TOT";
      mFile_ThermalEnergy << std::endl;
    }
    if (mGraph_HeatGenValues) {
      mFile_HeatGenValues.open("graph_heat_generation.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_HeatGenValues) << "Could not open graph file graph_heat_generation.txt!" << std::endl;
      mFile_HeatGenValues << "1 - TIME STEP | ";
      mFile_HeatGenValues << "2 - TIME | ";
      mFile_HeatGenValues << "3 - SLID PP | ";
      mFile_HeatGenValues << "4 - SLID PW | ";
      mFile_HeatGenValues << "5 - ROLL PP | ";
      mFile_HeatGenValues << "6 - ROLL PW | ";
      mFile_HeatGenValues << "7 - DAMP PP | ";
      mFile_HeatGenValues << "8 - DAMP PW | ";
      mFile_HeatGenValues << "9 - TOT";
      mFile_HeatGenValues << std::endl;
    }
    if (mGraph_HeatGenContrib) {
      mFile_HeatGenContrib.open("graph_heat_generation_relative.txt", std::ios::out);
      KRATOS_ERROR_IF_NOT(mFile_HeatGenContrib) << "Could not open graph file graph_heat_generation_relative.txt!" << std::endl;
      mFile_HeatGenContrib << "1 - TIME STEP | ";
      mFile_HeatGenContrib << "2 - TIME | ";
      mFile_HeatGenContrib << "3 - SLID PP | ";
      mFile_HeatGenContrib << "4 - SLID PW | ";
      mFile_HeatGenContrib << "5 - ROLL PP | ";
      mFile_HeatGenContrib << "6 - ROLL PW | ";
      mFile_HeatGenContrib << "7 - DAMP PP | ";
      mFile_HeatGenContrib << "8 - DAMP PW";
      mFile_HeatGenContrib << std::endl;
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

    const int num_particles     =  rModelPart.NumberOfElements();
    int       num_particles_gen =  0;
    double    total_vol         =  0.0;

    double particle_temp_min =  DBL_MAX;
    double particle_temp_max = -DBL_MAX;
    double particle_temp_avg =  0.0;
    double particle_temp_dev =  0.0;
    double model_temp_avg    =  0.0;

    double total_energy_potential_gravity   =  0.0;
    double total_energy_potential_elastic   =  0.0;
    double total_energy_kinetic_translation =  0.0;
    double total_energy_kinetic_rotation    =  0.0;
    double total_energy                     =  0.0;

    double total_dissip_slid         =  0.0;
    double total_dissip_roll         =  0.0;
    double total_dissip_damp         =  0.0;
    double total_dissip_damp_normal  =  0.0;
    double total_dissip_damp_tangent =  0.0;
    double total_dissip              =  0.0;

    double total_thermal_slid_pp = 0.0;
    double total_thermal_slid_pw = 0.0;
    double total_thermal_roll_pp = 0.0;
    double total_thermal_roll_pw = 0.0;
    double total_thermal_damp_pp = 0.0;
    double total_thermal_damp_pw = 0.0;
    double total_thermal         = 0.0;

    double total_heat_gen_slid_pp = 0.0;
    double total_heat_gen_slid_pw = 0.0;
    double total_heat_gen_roll_pp = 0.0;
    double total_heat_gen_roll_pw = 0.0;
    double total_heat_gen_damp_pp = 0.0;
    double total_heat_gen_damp_pw = 0.0;
    double total_heat_gen         = 0.0;

    double contrib_heat_gen_slid_pp = 0.0;
    double contrib_heat_gen_slid_pw = 0.0;
    double contrib_heat_gen_roll_pp = 0.0;
    double contrib_heat_gen_roll_pw = 0.0;
    double contrib_heat_gen_damp_pp = 0.0;
    double contrib_heat_gen_damp_pw = 0.0;

    ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();
    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < num_particles; i++) {
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*(it+i));

      const double vol  = particle.CalculateVolume();
      const double temp = particle.GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);

      double energy_potential_gravity         = 0.0;
      double energy_potential_elastic         = 0.0;
      double energy_kinetic_translation       = 0.0;
      double energy_kinetic_rotation          = 0.0;
      double dissip_slid                      = 0.0;
      double dissip_roll                      = 0.0;
      double dissip_damp                      = 0.0;
      double dissip_damp_normal               = 0.0;
      double dissip_damp_tangent              = 0.0;
      particle.Calculate(PARTICLE_GRAVITATIONAL_ENERGY,                energy_potential_gravity,   r_process_info);
      particle.Calculate(PARTICLE_ELASTIC_ENERGY,                      energy_potential_elastic,   r_process_info);
      particle.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY,      energy_kinetic_translation, r_process_info);
      particle.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY,         energy_kinetic_rotation,    r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_FRICTIONAL_ENERGY,         dissip_slid,                r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_ROLLING_RESISTANCE_ENERGY, dissip_roll,                r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_VISCODAMPING_ENERGY,       dissip_damp,                r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_DAMPING_NORMAL_ENERGY,     dissip_damp_normal,         r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_DAMPING_TANGENT_ENERGY,    dissip_damp_tangent,        r_process_info);

      const double thermal_slid_pp = fabs(particle.mGenerationThermalEnergy_slid_particle);
      const double thermal_slid_pw = fabs(particle.mGenerationThermalEnergy_slid_wall);
      const double thermal_roll_pp = fabs(particle.mGenerationThermalEnergy_roll_particle);
      const double thermal_roll_pw = fabs(particle.mGenerationThermalEnergy_roll_wall);
      const double thermal_damp_pp = fabs(particle.mGenerationThermalEnergy_damp_particle);
      const double thermal_damp_pw = fabs(particle.mGenerationThermalEnergy_damp_wall);

      const double heat_gen_slid_pp = fabs(particle.mGenerationHeatFlux_slid_particle);
      const double heat_gen_slid_pw = fabs(particle.mGenerationHeatFlux_slid_wall);
      const double heat_gen_roll_pp = fabs(particle.mGenerationHeatFlux_roll_particle);
      const double heat_gen_roll_pw = fabs(particle.mGenerationHeatFlux_roll_wall);
      const double heat_gen_damp_pp = fabs(particle.mGenerationHeatFlux_damp_particle);
      const double heat_gen_damp_pw = fabs(particle.mGenerationHeatFlux_damp_wall);
      const double heat_gen_total   = fabs(particle.mGenerationHeatFlux);

      #pragma omp critical
      {
        total_vol += vol;
        if (temp < particle_temp_min) particle_temp_min = temp;
        if (temp > particle_temp_max) particle_temp_max = temp;
        particle_temp_avg += temp;
        particle_temp_dev += temp * temp;
        model_temp_avg    += temp * vol;

        total_energy_potential_gravity   += energy_potential_gravity;
        total_energy_potential_elastic   += energy_potential_elastic;
        total_energy_kinetic_translation += energy_kinetic_translation;
        total_energy_kinetic_rotation    += energy_kinetic_rotation;
        total_energy                     += energy_potential_gravity + energy_potential_elastic + energy_kinetic_translation + energy_kinetic_rotation;
        
        total_dissip_slid         += dissip_slid;
        total_dissip_roll         += dissip_roll;
        total_dissip_damp         += dissip_damp;
        total_dissip_damp_normal  += dissip_damp_normal;
        total_dissip_damp_tangent += dissip_damp_tangent;
        total_dissip              += dissip_slid + dissip_roll + dissip_damp;

        total_thermal_slid_pp += thermal_slid_pp;
        total_thermal_slid_pw += thermal_slid_pw;
        total_thermal_roll_pp += thermal_roll_pp;
        total_thermal_roll_pw += thermal_roll_pw;
        total_thermal_damp_pp += thermal_damp_pp;
        total_thermal_damp_pw += thermal_damp_pw;
        total_thermal         += thermal_slid_pp + thermal_slid_pw + thermal_roll_pp + thermal_roll_pw + thermal_damp_pp + thermal_damp_pw;

        total_heat_gen_slid_pp += heat_gen_slid_pp;
        total_heat_gen_slid_pw += heat_gen_slid_pw;
        total_heat_gen_roll_pp += heat_gen_roll_pp;
        total_heat_gen_roll_pw += heat_gen_roll_pw;
        total_heat_gen_damp_pp += heat_gen_damp_pp;
        total_heat_gen_damp_pw += heat_gen_damp_pw;
        total_heat_gen         += heat_gen_slid_pp + heat_gen_slid_pw + heat_gen_roll_pp + heat_gen_roll_pw + heat_gen_damp_pp + heat_gen_damp_pw;

        if (mGraph_HeatGenContrib && heat_gen_total != 0.0) {
          num_particles_gen++;
          contrib_heat_gen_slid_pp += heat_gen_slid_pp / heat_gen_total;
          contrib_heat_gen_slid_pw += heat_gen_slid_pw / heat_gen_total;
          contrib_heat_gen_roll_pp += heat_gen_roll_pp / heat_gen_total;
          contrib_heat_gen_roll_pw += heat_gen_roll_pw / heat_gen_total;
          contrib_heat_gen_damp_pp += heat_gen_damp_pp / heat_gen_total;
          contrib_heat_gen_damp_pw += heat_gen_damp_pw / heat_gen_total;
        }
      }
    } // Loop over particles

    particle_temp_avg /= num_particles;
    model_temp_avg    /= total_vol;
    particle_temp_dev  = sqrt(std::max(0.0, particle_temp_dev / num_particles - particle_temp_avg * particle_temp_avg));

    if (mGraph_HeatGenContrib && num_particles_gen > 0) {
      contrib_heat_gen_slid_pp /= num_particles_gen;
      contrib_heat_gen_slid_pw /= num_particles_gen;
      contrib_heat_gen_roll_pp /= num_particles_gen;
      contrib_heat_gen_roll_pw /= num_particles_gen;
      contrib_heat_gen_damp_pp /= num_particles_gen;
      contrib_heat_gen_damp_pw /= num_particles_gen;
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

    if (mFile_MechanicalEnergy.is_open())
      mFile_MechanicalEnergy << time_step                        << " "
                             << time                             << " "
                             << total_energy_potential_gravity   << " "
                             << total_energy_potential_elastic   << " "
                             << total_energy_kinetic_translation << " "
                             << total_energy_kinetic_rotation    << " "
                             << total_energy                     << " "
                             << total_dissip_slid                << " "
                             << total_dissip_roll                << " "
                             << total_dissip_damp                << " "
                             << total_dissip_damp_normal         << " "
                             << total_dissip_damp_tangent        << " "
                             << total_dissip                     << " "
                             << total_energy+total_dissip
                             << std::endl;

    if (mFile_ThermalEnergy.is_open())
      mFile_ThermalEnergy << time_step              << " "
                          << time                   << " "
                          << total_thermal_slid_pp << " "
                          << total_thermal_slid_pw << " "
                          << total_thermal_roll_pp << " "
                          << total_thermal_roll_pw << " "
                          << total_thermal_damp_pp << " "
                          << total_thermal_damp_pw << " "
                          << total_thermal
                          << std::endl;

    if (mFile_HeatGenValues.is_open())
      mFile_HeatGenValues << time_step              << " "
                          << time                   << " "
                          << total_heat_gen_slid_pp << " "
                          << total_heat_gen_slid_pw << " "
                          << total_heat_gen_roll_pp << " "
                          << total_heat_gen_roll_pw << " "
                          << total_heat_gen_damp_pp << " "
                          << total_heat_gen_damp_pw << " "
                          << total_heat_gen
                          << std::endl;

    if (mFile_HeatGenContrib.is_open())
      mFile_HeatGenContrib << time_step                << " "
                           << time                     << " "
                           << contrib_heat_gen_slid_pp << " "
                           << contrib_heat_gen_slid_pw << " "
                           << contrib_heat_gen_roll_pp << " "
                           << contrib_heat_gen_roll_pw << " "
                           << contrib_heat_gen_damp_pp << " "
                           << contrib_heat_gen_damp_pw
                           << std::endl;

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void GraphUtilities::ExecuteFinalize(void)
  {
    KRATOS_TRY

    if (mFile_ParticleTempMin.is_open())   mFile_ParticleTempMin.close();
    if (mFile_ParticleTempMax.is_open())   mFile_ParticleTempMax.close();
    if (mFile_ParticleTempAvg.is_open())   mFile_ParticleTempAvg.close();
    if (mFile_ParticleTempDev.is_open())   mFile_ParticleTempDev.close();
    if (mFile_ModelTempAvg.is_open())      mFile_ModelTempAvg.close();
    if (mFile_MechanicalEnergy.is_open())  mFile_MechanicalEnergy.close();
    if (mFile_ThermalEnergy.is_open())     mFile_ThermalEnergy.close();
    if (mFile_HeatGenValues.is_open())     mFile_HeatGenValues.close();
    if (mFile_HeatGenContrib.is_open())    mFile_HeatGenContrib.close();

    KRATOS_CATCH("")
  }

} // namespace Kratos
