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
#include "heat_map_utilities.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  HeatMapUtilities::HeatMapUtilities() {}
  HeatMapUtilities::~HeatMapUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void HeatMapUtilities::ExecuteInitialize(ModelPart& rModelPart) {
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    if (!r_process_info[HEAT_MAP_GENERATION_OPTION])
      return;

    // Store properties
    array_1d<int, 3> subdivisions = r_process_info[HEAT_MAP_SUBDIVISIONS];
    mDimX = subdivisions[0];
    mDimY = subdivisions[1];
    mDimZ = subdivisions[2];

    // Initialize global maps with zeros
    ResetMap(mGlobalHeatMapGenerationDampingPP);
    ResetMap(mGlobalHeatMapGenerationDampingPW);
    ResetMap(mGlobalHeatMapGenerationSlidingPP);
    ResetMap(mGlobalHeatMapGenerationSlidingPW);
    ResetMap(mGlobalHeatMapGenerationRollingPP);
    ResetMap(mGlobalHeatMapGenerationRollingPW);

    // Initialize local maps with zeros
    const int num_of_particles = rModelPart.NumberOfElements();
    ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < num_of_particles; i++) {
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*(it+i));
      ResetMap(particle.mHeatMapGenerationDampingPP);
      ResetMap(particle.mHeatMapGenerationDampingPW);
      ResetMap(particle.mHeatMapGenerationSlidingPP);
      ResetMap(particle.mHeatMapGenerationSlidingPW);
      ResetMap(particle.mHeatMapGenerationRollingPP);
      ResetMap(particle.mHeatMapGenerationRollingPW);
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  // ATTENTION: IF NUMBER OF PARTICLES IS CONSTANT, THIS COULD BE DONE ONLY ONCE IN THE END OF ANALYSIS
  void HeatMapUtilities::ExecuteFinalizeSolutionStep(ModelPart& rModelPart) {
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    if (!r_process_info[IS_TIME_TO_PRINT] || !r_process_info[HEAT_MAP_GENERATION_OPTION])
      return;

    // Merge local maps to global maps and reset local maps with zeros
    const int num_of_particles = rModelPart.NumberOfElements();
    ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < num_of_particles; i++) {
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*(it + i));
      
      // Merge local maps to global maps and reset local maps with zeros
      for (int i = 0; i < mDimX; i++) {
        for (int j = 0; j < mDimY; j++) {
          for (int k = 0; k < mDimZ; k++) {
            mGlobalHeatMapGenerationDampingPP[i][j][k] += particle.mHeatMapGenerationDampingPP[i][j][k];
            mGlobalHeatMapGenerationDampingPW[i][j][k] += particle.mHeatMapGenerationDampingPW[i][j][k];
            mGlobalHeatMapGenerationSlidingPP[i][j][k] += particle.mHeatMapGenerationSlidingPP[i][j][k];
            mGlobalHeatMapGenerationSlidingPW[i][j][k] += particle.mHeatMapGenerationSlidingPW[i][j][k];
            mGlobalHeatMapGenerationRollingPP[i][j][k] += particle.mHeatMapGenerationRollingPP[i][j][k];
            mGlobalHeatMapGenerationRollingPW[i][j][k] += particle.mHeatMapGenerationRollingPW[i][j][k];

            particle.mHeatMapGenerationDampingPP[i][j][k] = 0.0;
            particle.mHeatMapGenerationDampingPW[i][j][k] = 0.0;
            particle.mHeatMapGenerationSlidingPP[i][j][k] = 0.0;
            particle.mHeatMapGenerationSlidingPW[i][j][k] = 0.0;
            particle.mHeatMapGenerationRollingPP[i][j][k] = 0.0;
            particle.mHeatMapGenerationRollingPW[i][j][k] = 0.0;
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void HeatMapUtilities::ExecuteFinalize(ModelPart& rModelPart) {
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    if (!r_process_info[HEAT_MAP_GENERATION_OPTION])
      return;

    const array_1d<double, 3> coords_1     = r_process_info[HEAT_MAP_COORDINATES_1];
    const array_1d<double, 3> coords_2     = r_process_info[HEAT_MAP_COORDINATES_2];
    const array_1d<int, 3>    subdivisions = r_process_info[HEAT_MAP_SUBDIVISIONS];

    // Open files
    std::ofstream file_generation_damping_pp;
    std::ofstream file_generation_damping_pw;
    std::ofstream file_generation_sliding_pp;
    std::ofstream file_generation_sliding_pw;
    std::ofstream file_generation_rolling_pp;
    std::ofstream file_generation_rolling_pw;

    file_generation_damping_pp.open("map_generation_damping_pp.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(file_generation_damping_pp) << "Could not open file map_generation_damping_pp.txt!" << std::endl;
    file_generation_damping_pw.open("map_generation_damping_pw.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(file_generation_damping_pw) << "Could not open file map_generation_damping_pw.txt!" << std::endl;
    file_generation_sliding_pp.open("map_generation_sliding_pp.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(file_generation_sliding_pp) << "Could not open file map_generation_sliding_pp.txt!" << std::endl;
    file_generation_sliding_pw.open("map_generation_sliding_pw.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(file_generation_sliding_pw) << "Could not open file map_generation_sliding_pw.txt!" << std::endl;
    file_generation_rolling_pp.open("map_generation_rolling_pp.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(file_generation_rolling_pp) << "Could not open file map_generation_rolling_pp.txt!" << std::endl;
    file_generation_rolling_pw.open("map_generation_rolling_pw.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(file_generation_rolling_pw) << "Could not open file map_generation_rolling_pw.txt!" << std::endl;

    // Print headers
    file_generation_damping_pp << "MAP - HEAT GENERATION - DAMPING - PARTICLE-PARTICLE" << std::endl;
    file_generation_damping_pp << coords_1[0]     << " " << coords_1[1]     << " " << coords_1[2]     << std::endl;
    file_generation_damping_pp << coords_2[0]     << " " << coords_2[1]     << " " << coords_2[2]     << std::endl;
    file_generation_damping_pp << subdivisions[0] << " " << subdivisions[1] << " " << subdivisions[2] << std::endl;

    file_generation_damping_pw << "MAP - HEAT GENERATION - DAMPING - PARTICLE-WALL" << std::endl;
    file_generation_damping_pw << coords_1[0]     << " " << coords_1[1]     << " " << coords_1[2]     << std::endl;
    file_generation_damping_pw << coords_2[0]     << " " << coords_2[1]     << " " << coords_2[2]     << std::endl;
    file_generation_damping_pw << subdivisions[0] << " " << subdivisions[1] << " " << subdivisions[2] << std::endl;

    file_generation_sliding_pp << "MAP - HEAT GENERATION - SLIDING - PARTICLE-PARTICLE" << std::endl;
    file_generation_sliding_pp << coords_1[0]     << " " << coords_1[1]     << " " << coords_1[2]     << std::endl;
    file_generation_sliding_pp << coords_2[0]     << " " << coords_2[1]     << " " << coords_2[2]     << std::endl;
    file_generation_sliding_pp << subdivisions[0] << " " << subdivisions[1] << " " << subdivisions[2] << std::endl;

    file_generation_sliding_pw << "MAP - HEAT GENERATION - SLIDING - PARTICLE-WALL" << std::endl;
    file_generation_sliding_pw << coords_1[0]     << " " << coords_1[1]     << " " << coords_1[2]     << std::endl;
    file_generation_sliding_pw << coords_2[0]     << " " << coords_2[1]     << " " << coords_2[2]     << std::endl;
    file_generation_sliding_pw << subdivisions[0] << " " << subdivisions[1] << " " << subdivisions[2] << std::endl;

    file_generation_rolling_pp << "MAP - HEAT GENERATION - ROLLING - PARTICLE-PARTICLE" << std::endl;
    file_generation_rolling_pp << coords_1[0]     << " " << coords_1[1]     << " " << coords_1[2]     << std::endl;
    file_generation_rolling_pp << coords_2[0]     << " " << coords_2[1]     << " " << coords_2[2]     << std::endl;
    file_generation_rolling_pp << subdivisions[0] << " " << subdivisions[1] << " " << subdivisions[2] << std::endl;

    file_generation_rolling_pw << "MAP - HEAT GENERATION - ROLLING - PARTICLE-WALL" << std::endl;
    file_generation_rolling_pw << coords_1[0]     << " " << coords_1[1]     << " " << coords_1[2]     << std::endl;
    file_generation_rolling_pw << coords_2[0]     << " " << coords_2[1]     << " " << coords_2[2]     << std::endl;
    file_generation_rolling_pw << subdivisions[0] << " " << subdivisions[1] << " " << subdivisions[2] << std::endl;

    // Print global maps
    file_generation_damping_pp << std::endl;
    file_generation_damping_pw << std::endl;
    file_generation_sliding_pp << std::endl;
    file_generation_sliding_pw << std::endl;
    file_generation_rolling_pp << std::endl;
    file_generation_rolling_pw << std::endl;

    for (int k = 0; k < mDimZ; k++) {
      for (int j = 0; j < mDimY; j++) {
        for (int i = 0; i < mDimX; i++) {
          file_generation_damping_pp << mGlobalHeatMapGenerationDampingPP[i][j][k] << " ";
          file_generation_damping_pw << mGlobalHeatMapGenerationDampingPW[i][j][k] << " ";
          file_generation_sliding_pp << mGlobalHeatMapGenerationSlidingPP[i][j][k] << " ";
          file_generation_sliding_pw << mGlobalHeatMapGenerationSlidingPW[i][j][k] << " ";
          file_generation_rolling_pp << mGlobalHeatMapGenerationRollingPP[i][j][k] << " ";
          file_generation_rolling_pw << mGlobalHeatMapGenerationRollingPW[i][j][k] << " ";
        }
        file_generation_damping_pp << std::endl;
        file_generation_damping_pw << std::endl;
        file_generation_sliding_pp << std::endl;
        file_generation_sliding_pw << std::endl;
        file_generation_rolling_pp << std::endl;
        file_generation_rolling_pw << std::endl;
      }
    }

    // Close files
    if (file_generation_damping_pp.is_open()) file_generation_damping_pp.close();
    if (file_generation_damping_pw.is_open()) file_generation_damping_pw.close();
    if (file_generation_sliding_pp.is_open()) file_generation_sliding_pp.close();
    if (file_generation_sliding_pw.is_open()) file_generation_sliding_pw.close();
    if (file_generation_rolling_pp.is_open()) file_generation_rolling_pp.close();
    if (file_generation_rolling_pw.is_open()) file_generation_rolling_pw.close();
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void HeatMapUtilities::ResetMap(std::vector<std::vector<std::vector<double>>>& map) {
    for (int i = 0; i < mDimX; i++) {
      std::vector<std::vector<double>> vv;
      for (int j = 0; j < mDimY; j++) {
        std::vector<double> v;
        for (int k = 0; k < mDimZ; k++) {
          v.push_back(0.0);
        }
        vv.push_back(v);
      }
      map.push_back(vv);
    }
  }

} // namespace Kratos
