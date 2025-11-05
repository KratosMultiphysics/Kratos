#include "gpu_results_utilities.h"

#define FREQUENCY_TIME 0.01

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  void GPUResultsUtilities::ExecuteInitialize(ModelPart& particlesMP, ModelPart& wallsMP) {
    mTimeToWrite = FREQUENCY_TIME;

    // Open permanent files
    mFile_Global.open("data_global.txt", std::ios::out);
    mFile_Global << "1 - TIME | " << "2 - ENERGY_ELASTIC |" << std::endl;

    // Get radii
    mNumParticles = particlesMP.GetCommunicator().LocalMesh().Elements().size();
    double radius_max = 0.0;
    double radius_min = std::numeric_limits<double>::max();
    for (int i = 0; i < mNumParticles; i++) {
      ModelPart::ElementsContainerType::iterator it = particlesMP.GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
      SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
      const double r = particle.GetRadius();
      radius_max = std::max(radius_max, r);
      radius_min = std::min(radius_min, r);
    }
    mRadiusMean = (radius_max + radius_min) / 2.0;

    // Get number of walls
    mNumWalls = wallsMP.GetCommunicator().LocalMesh().Conditions().size();
  }

  //------------------------------------------------------------------------------------------------------------
  void GPUResultsUtilities::Calculate(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Check if it is time to execute
    const ProcessInfo& r_process_info = particlesMP.GetProcessInfo();
    const double time = r_process_info[TIME];
    if (time < mTimeToWrite) return;
    mTimeToWrite += FREQUENCY_TIME;

    // Initialize variables
    double energy_elast = 0.0;
    
    // Open temporary files
    std::ostringstream particles_filename;
    particles_filename << "data_particles_" << std::fixed << std::setprecision(2) << time << ".txt";
    std::ofstream file_particles_info(particles_filename.str(), std::ios::out);
    file_particles_info << "1 - ID | 2 - RADIUS | 3 - COORD_X | 4 - COORD_Y | 5 - VEL_X | 6 - VEL_Y | 7 - VEL_ROT | 8 - NUM_NEIGHBORS_ALL | 9 - NUM_NEIGHBORS_PART | 10 - NUM_NEIGHBORS_WALL |" << std::endl;

    // Loop over all particles
    ModelPart::ElementsContainerType::iterator itp = particlesMP.GetCommunicator().LocalMesh().Elements().ptr_begin();
    for (int i = 0; i < mNumParticles; i++) {
      // Get particle properties
      SphericParticle& particle = dynamic_cast<SphericParticle&>(*(itp+i));
      const auto& node = particle.GetGeometry()[0];
      const int id = particle.GetId();
      const double radius = particle.GetRadius();
      const array_1d<double, 3>& coords = node.Coordinates();
      const double coord_x = coords[0];
      const double coord_y = coords[1];
      const double coord_z = coords[2];
      const int num_neighbors_p = particle.mNeighbourElements.size();
      const int num_neighbors_w = particle.mNeighbourRigidFaces.size();
      const int num_neighbors = num_neighbors_p + num_neighbors_w;
      const array_1d<double, 3>& vel_trl = node.FastGetSolutionStepValue(VELOCITY);
      const array_1d<double, 3>& vel_rot = node.FastGetSolutionStepValue(ANGULAR_VELOCITY);
      const double vel_trl_x = vel_trl[0];
      const double vel_trl_y = vel_trl[1];
      const double vel_trl_z = vel_trl[2];
      const double vel_rot_x = vel_rot[0];
      const double vel_rot_y = vel_rot[1];
      const double vel_rot_z = vel_rot[2];
      double particle_energy_elast = 0.0;
      particle.Calculate(PARTICLE_ELASTIC_ENERGY, particle_energy_elast, r_process_info);
      energy_elast += particle_energy_elast;
      
      // Write particle information to file
      if (file_particles_info.is_open()) {
        file_particles_info
        << std::fixed << std::setprecision(0)
        << id << " "
        << std::fixed << std::setprecision(15)
        << radius << " " << coord_x << " " << coord_y << " "
        << vel_trl_x << " " << vel_trl_y << " " << vel_rot_z << " "
        << std::fixed << std::setprecision(0)
        << num_neighbors << " " << num_neighbors_p << " " << num_neighbors_w << " "
        << std::endl;
      }
    }

    // Write global results
    if (mFile_Global.is_open()) {
      mFile_Global
      << std::fixed << std::setprecision(6)
      << time << " "
      << std::fixed << std::setprecision(15)
      << energy_elast << " "
      << std::endl;
    }

    // Close temporary files
    if (file_particles_info.is_open()) file_particles_info.close();
  }

  //------------------------------------------------------------------------------------------------------------
  void GPUResultsUtilities::ExecuteFinalize(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Close permanent files
    if (mFile_Global.is_open()) mFile_Global.close();
  }
}