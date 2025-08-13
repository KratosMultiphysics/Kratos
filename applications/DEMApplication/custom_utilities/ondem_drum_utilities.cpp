#include "ondem_drum_utilities.h"

#define FREQUENCY_TIME 0.10
#define FREQUENCY_TOL 1e-9
#define DOMAIN_MIN_X -0.10
#define DOMAIN_MIN_Y -0.10
#define DOMAIN_MIN_Z -0.03
#define DOMAIN_MAX_X 0.10
#define DOMAIN_MAX_Y 0.10
#define DOMAIN_MAX_Z 0.03
#define CELL_SIZE 0.010
#define CELL_MIN_PARTICLES_SMALL 30

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ExecuteInitialize(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Open permanent files
    mFile_Contact.open("data_contact.txt", std::ios::out);
    mFile_Contact << "1 - TIME | " << "2 - ALL | " << "3 - P-P | " << "4 - PL-PL | " << "5 - PS-PS | " << "6 - PL-PS | " << "7 - P-W | " << "8 - PL-W | " << "9 - PS-W | " << "10 - CN_ALL | " << "11 - CN_LARGE | " << "12 - CN_SMALL" << std::endl;

    mFile_Velocity.open("data_velocity.txt", std::ios::out);
    mFile_Velocity << "1 - TIME | " << "2 - AVG TRANSLATION (ALL) | " << "3 - AVG TRANSLATION (LARGE) | " << "4 - AVG TRANSLATION (SMALL) | " << "5 - AVG ROTATION (ALL) | " << "6 - AVG ROTATION (LARGE) | " << "7 - AVG ROTATION (SMALL)" << std::endl;

    mFile_Energy.open("data_energy.txt", std::ios::out);
    mFile_Energy << "1 - TIME | " << "2 - MECHANICAL | " << "3 - KINETIC TRL | " << "4 - KINETIC ROT | " << "5 - ELASTIC | " << "6 - GRAVITY | " << "7 - DISSIP | " << "8 - DISSIP FRICTION | " << "9 - DISSIP DAMPING" << std::endl;

    // Get radii
    mNumParticlesAll = particlesMP.GetCommunicator().LocalMesh().Elements().size();
    double radius_max = 0.0;
    double radius_min = std::numeric_limits<double>::max();
    for (int i = 0; i < mNumParticlesAll; i++) {
      ModelPart::ElementsContainerType::iterator it = particlesMP.GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
      SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
      const double r = particle.GetRadius();
      radius_max = std::max(radius_max, r);
      radius_min = std::min(radius_min, r);
    }
    mRadiusMean = (radius_max + radius_min) / 2.0;

    // Get number of elements
    mNumWalls = wallsMP.GetCommunicator().LocalMesh().Conditions().size();
    mNumParticlesLarge = 0;
    mNumParticlesSmall = 0;
    for (int i = 0; i < mNumParticlesAll; i++) {
      ModelPart::ElementsContainerType::iterator it = particlesMP.GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
      SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
      const double r = particle.GetRadius();
      if (r > mRadiusMean) mNumParticlesLarge++;
      else mNumParticlesSmall++;
    }

    // Initialize maps
    mCellsX = static_cast<int>((DOMAIN_MAX_X - DOMAIN_MIN_X) / CELL_SIZE);
    mCellsY = static_cast<int>((DOMAIN_MAX_Y - DOMAIN_MIN_Y) / CELL_SIZE);
    mCellsZ = static_cast<int>((DOMAIN_MAX_Z - DOMAIN_MIN_Z) / CELL_SIZE);
    InitializeCellsMap(mCountParticlesLarge);
    InitializeCellsMap(mCountParticlesSmall);
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::Calculate(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Check if it is time to execute
    const ProcessInfo& r_process_info = particlesMP.GetProcessInfo();
    const double time = r_process_info[TIME];
    double resid = std::fmod(time, FREQUENCY_TIME);
    if (resid > FREQUENCY_TOL && std::fabs(resid - FREQUENCY_TIME) > FREQUENCY_TOL) return;

    // Initialize variables
    int num_contacts_all   = 0;
    int num_contacts_p_p   = 0;
    int num_contacts_pl_pl = 0;
    int num_contacts_ps_ps = 0;
    int num_contacts_pl_ps = 0;
    int num_contacts_p_w   = 0;
    int num_contacts_pl_w  = 0;
    int num_contacts_ps_w  = 0;
    double coord_num_all   = 0.0;
    double coord_num_pl    = 0.0;
    double coord_num_ps    = 0.0;
    double vel_trl_all_avg = 0.0;
    double vel_trl_pl_avg  = 0.0;
    double vel_trl_ps_avg  = 0.0;
    double vel_rot_all_avg = 0.0;
    double vel_rot_pl_avg  = 0.0;
    double vel_rot_ps_avg  = 0.0;
    double energy_kint     = 0.0;
    double energy_kinr     = 0.0;
    double energy_elast    = 0.0;
    double energy_grav     = 0.0;
    double energy_mech     = 0.0;
    double dissip_fric     = 0.0;
    double dissip_damp     = 0.0;
    double dissip          = 0.0;
    ResetCellsMap(mCountParticlesLarge);
    ResetCellsMap(mCountParticlesSmall);

    // Open temporary files
    std::ostringstream walls_filename;
    walls_filename << "data_walls_" << std::fixed << std::setprecision(6) << time << ".txt";
    std::ofstream file_walls_info(walls_filename.str());
    file_walls_info << "ID X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3" << std::endl;

    std::ostringstream particles_filename;
    particles_filename << "data_particles_" << std::fixed << std::setprecision(6) << time << ".txt";
    std::ofstream file_particles_info(particles_filename.str());
    file_particles_info << "1 - ID | 2 - RADIUS | 3 - X | 4 - Y | 5 - Z | 6 - CONTACTS_P | 7 - CONTACTS_W | 8 - VEL_TRL_X | 9 - VEL_TRL_Y | 10 - VEL_TRL_Z | 11 - VEL_ROT_X | 12 - VEL_ROT_Y | 13 - VEL_ROT_Z | 14 - ENERGY_KINT | 15 - ENERGY_KINR | 16 - ENERGY_ELAST | 17 - ENERGY_GRAV | 18 - DISSIP_FRIC | 19 - DISSIP_DAMP" << std::endl;

    std::ostringstream cells_filename;
    cells_filename << "data_cells_" << std::fixed << std::setprecision(6) << time << ".txt";
    std::ofstream file_cells_info(cells_filename.str());
    file_cells_info << "[XMIN XMAX YMIN YMAX ZMIN ZMAX]: #PL #PS" << std::endl;

    // Get pointers to first elements
    ModelPart::ConditionsContainerType::iterator itw = wallsMP.GetCommunicator().LocalMesh().Conditions().ptr_begin();
    ModelPart::ElementsContainerType::iterator itp = particlesMP.GetCommunicator().LocalMesh().Elements().ptr_begin();

    // Loop over all walls
    for (int i = 0; i < mNumWalls; i++) {
      // Get wall properties
      DEMWall& wall = dynamic_cast<DEMWall&>(*(itw+i));
      const int id = wall.GetId();
      Condition::GeometryType &geom = wall.GetGeometry();
      const double x1 = geom[0][0]; const double y1 = geom[0][1]; const double z1 = geom[0][2];
      const double x2 = geom[1][0]; const double y2 = geom[1][1]; const double z2 = geom[1][2];
      const double x3 = geom[2][0]; const double y3 = geom[2][1]; const double z3 = geom[2][2];

      // Write wall information to file
      if (file_walls_info.is_open()) {
        file_walls_info << std::fixed << std::setprecision(0) << id << " " << std::fixed << std::setprecision(16) << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << " " << x3 << " " << y3 << " " << z3 << std::endl;
      }
    }

    // Loop over all particles
    for (int i = 0; i < mNumParticlesAll; i++) {
      // Get particle properties
      SphericParticle& particle = dynamic_cast<SphericParticle&>(*(itp+i));
      const auto& node1 = particle.GetGeometry()[0];
      const int id1 = particle.GetId();
      const double r1 = particle.GetRadius();
      const array_1d<double, 3>& coords1 = node1.Coordinates();
      const double x1 = coords1[0];
      const double y1 = coords1[1];
      const double z1 = coords1[2];
      const int num_neighbors_p = particle.mNeighbourElements.size();
      const int num_neighbors_w = particle.mNeighbourRigidFaces.size();
      const int num_neighbors = num_neighbors_p + num_neighbors_w;
      const array_1d<double, 3>& vel_trl = node1.FastGetSolutionStepValue(VELOCITY);
      const array_1d<double, 3>& vel_rot = node1.FastGetSolutionStepValue(ANGULAR_VELOCITY);
      const double vel_trl_norm = DEM_MODULUS_3(vel_trl);
      const double vel_rot_norm = DEM_MODULUS_3(vel_rot);
      const double vel_trl_x = vel_trl[0];
      const double vel_trl_y = vel_trl[1];
      const double vel_trl_z = vel_trl[2];
      const double vel_rot_x = vel_rot[0];
      const double vel_rot_y = vel_rot[1];
      const double vel_rot_z = vel_rot[2];
      double particle_energy_kint = 0.0;
      double particle_energy_kinr = 0.0;
      double particle_energy_elast = 0.0;
      double particle_energy_grav = 0.0;
      double particle_dissip_fric = 0.0;
      double particle_dissip_damp = 0.0;
      particle.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY, particle_energy_kint,  r_process_info);
      particle.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY,    particle_energy_kinr,  r_process_info);
      particle.Calculate(PARTICLE_ELASTIC_ENERGY,                 particle_energy_elast, r_process_info);
      particle.Calculate(PARTICLE_GRAVITATIONAL_ENERGY,           particle_energy_grav,  r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_FRICTIONAL_ENERGY,    particle_dissip_fric,  r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_VISCODAMPING_ENERGY,  particle_dissip_damp,  r_process_info);
      
      // Write particle information to file
      if (file_particles_info.is_open()) {
        file_particles_info
        << std::fixed << std::setprecision(0)
        << id1 << " "
        << std::fixed << std::setprecision(16)
        << r1 << " " << x1 << " " << y1 << " " << z1 << " "
        << std::fixed << std::setprecision(0)
        << num_neighbors_p << " " << num_neighbors_w << " "
        << std::fixed << std::setprecision(16)
        << vel_trl_x << " " << vel_trl_y << " " << vel_trl_z << " "
        << vel_rot_x << " " << vel_rot_y << " " << vel_rot_z << " "
        << particle_energy_kint << " " << particle_energy_kinr << " " << particle_energy_elast << " " << particle_energy_grav << " "
        << particle_dissip_fric << " " << particle_dissip_damp
        << std::endl;
      }

      // Update contact statistics
      num_contacts_all += num_neighbors_w;
      num_contacts_p_w += num_neighbors_w;
      coord_num_all += num_neighbors;
      if (r1 > mRadiusMean) {
        num_contacts_pl_w += num_neighbors_w;
        coord_num_pl += num_neighbors;
      }
      else {
        num_contacts_ps_w += num_neighbors_w;
        coord_num_ps += num_neighbors;
      }

      // Loop over particle neighbors
      for (int j = 0; j < num_neighbors_p; j++) {
        // Get particle properties
        SphericParticle& neighbor = dynamic_cast<SphericParticle&>(*particle.mNeighbourElements[j]);
        const int id2 = neighbor.GetId();
        if (id2 <= id1) continue;
        const double r2 = neighbor.GetRadius();
        const auto& node2 = neighbor.GetGeometry()[0];
        const array_1d<double, 3>& coords2 = node2.Coordinates();
        const double x2 = coords2[0];
        const double y2 = coords2[1];
        const double z2 = coords2[2];

        // Check if neighbor is in contact
        const double dist = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        if (dist < (r1 + r2)) {
          // Update contact statistics
          num_contacts_all++;
          num_contacts_p_p++;
          if (r1 > mRadiusMean && r2 > mRadiusMean) 
            num_contacts_pl_pl++;
          else if (r1 < mRadiusMean && r2 < mRadiusMean)
            num_contacts_ps_ps++;
          else
            num_contacts_pl_ps++;
        }
      }
      
      // Accumulate velocities
      vel_trl_all_avg += vel_trl_norm;
      vel_rot_all_avg += vel_rot_norm;
      if (r1 > mRadiusMean) {
        vel_trl_pl_avg += vel_trl_norm;
        vel_rot_pl_avg += vel_rot_norm;
      }
      else {
        vel_trl_ps_avg += vel_trl_norm;
        vel_rot_ps_avg += vel_rot_norm;
      }

      // Accumulate energies
      energy_kint  += particle_energy_kint;
      energy_kinr  += particle_energy_kinr;
      energy_elast += particle_energy_elast;
      energy_grav  += particle_energy_grav;
      energy_mech  += particle_energy_kint + particle_energy_kinr + particle_energy_elast + particle_energy_grav;
      dissip_fric  += particle_dissip_fric;
      dissip_damp  += particle_dissip_damp;
      dissip       += particle_dissip_fric + particle_dissip_damp;

      // Determine cells map indices of current particle
      const int i_x = std::max(0, std::min(mCellsX - 1, static_cast<int>(std::floor((x1 - DOMAIN_MIN_X) / CELL_SIZE))));
      const int i_y = std::max(0, std::min(mCellsY - 1, static_cast<int>(std::floor((y1 - DOMAIN_MIN_Y) / CELL_SIZE))));
      const int i_z = std::max(0, std::min(mCellsZ - 1, static_cast<int>(std::floor((z1 - DOMAIN_MIN_Z) / CELL_SIZE))));

      // Update cells map
      if (r1 > mRadiusMean) mCountParticlesLarge[i_x][i_y][i_z]++;
      else mCountParticlesSmall[i_x][i_y][i_z]++;
    }

    coord_num_all   /= mNumParticlesAll;
    coord_num_pl    /= mNumParticlesLarge;
    coord_num_ps    /= mNumParticlesSmall;
    vel_trl_all_avg /= mNumParticlesAll;
    vel_trl_pl_avg  /= mNumParticlesLarge;
    vel_trl_ps_avg  /= mNumParticlesSmall;
    vel_rot_all_avg /= mNumParticlesAll;
    vel_rot_pl_avg  /= mNumParticlesLarge;
    vel_rot_ps_avg  /= mNumParticlesSmall;

    // Write results
    if (file_cells_info.is_open()) {
      for (int iz = 0; iz < mCellsZ; ++iz) {
        double zmin = DOMAIN_MIN_Z + CELL_SIZE * iz;
        double zmax = DOMAIN_MIN_Z + CELL_SIZE * (iz + 1);
        for (int iy = 0; iy < mCellsY; ++iy) {
          double ymin = DOMAIN_MIN_Y + CELL_SIZE * iy;
          double ymax = DOMAIN_MIN_Y + CELL_SIZE * (iy + 1);
          for (int ix = 0; ix < mCellsX; ++ix) {
            double xmin = DOMAIN_MIN_X + CELL_SIZE * ix;
            double xmax = DOMAIN_MIN_X + CELL_SIZE * (ix + 1);
            file_cells_info << std::fixed << std::setprecision(2);
            file_cells_info << "[" << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << "]: ";
            file_cells_info << std::fixed << std::setprecision(0);
            file_cells_info << mCountParticlesLarge[ix][iy][iz] << " " << mCountParticlesSmall[ix][iy][iz] << std::endl;
          }
        }
      }
    }
    
    if (mFile_Contact.is_open()) {
      mFile_Contact
      << std::fixed << std::setprecision(6)
      << time << " "
      << std::fixed << std::setprecision(0)
      << num_contacts_all   << " "
      << num_contacts_p_p   << " "
      << num_contacts_pl_pl << " "
      << num_contacts_pl_ps << " "
      << num_contacts_ps_ps << " "
      << num_contacts_p_w   << " "
      << num_contacts_pl_w  << " "
      << num_contacts_ps_w  << " "
      << std::fixed << std::setprecision(15)
      << coord_num_all << " "
      << coord_num_pl  << " "
      << coord_num_ps  << " "
      << std::endl;
    }
    if (mFile_Velocity.is_open()) {
      mFile_Velocity 
      << std::fixed << std::setprecision(6)
      << time << " "
      << std::fixed << std::setprecision(15)
      << vel_trl_all_avg << " "
      << vel_trl_pl_avg  << " "
      << vel_trl_ps_avg  << " "
      << vel_rot_all_avg << " "
      << vel_rot_pl_avg  << " "
      << vel_rot_ps_avg  << " "
      << std::endl;
    }
    if (mFile_Energy.is_open()) {
      mFile_Energy
      << std::fixed << std::setprecision(6)
      << time << " "
      << std::fixed << std::setprecision(15)
      << energy_mech  << " "
      << energy_kint  << " "
      << energy_kinr  << " "
      << energy_elast << " "
      << energy_grav  << " "
      << dissip       << " "
      << dissip_fric  << " "
      << dissip_damp  << " "
      << std::endl;
    }

    // Close temporary files
    if (file_walls_info.is_open()) file_walls_info.close();
    if (file_particles_info.is_open()) file_particles_info.close();
    if (file_cells_info.is_open()) file_cells_info.close();
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ExecuteFinalize(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Close permanent files
    if (mFile_Contact.is_open()) mFile_Contact.close();
    if (mFile_Velocity.is_open()) mFile_Velocity.close();
    if (mFile_Energy.is_open()) mFile_Energy.close();
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::InitializeCellsMap(std::vector<std::vector<std::vector<int>>>& map) {
    map.assign(mCellsX, std::vector<std::vector<int>>(mCellsY, std::vector<int>(mCellsZ, 0)));
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ResetCellsMap(std::vector<std::vector<std::vector<int>>>& map) {
    for (auto& plane : map)
      for (auto& row : plane)
        std::fill(row.begin(), row.end(), 0.0);
  }
}