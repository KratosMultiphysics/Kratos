#include "ondem_drum_utilities.h"

#define FREQUENCY_TIME_1 0.10
#define FREQUENCY_TIME_2 0.25
#define FREQUENCY_TOL 1e-9
#define DOMAIN_MIN_X -0.10
#define DOMAIN_MIN_Y -0.10
#define DOMAIN_MIN_Z -0.03
#define DOMAIN_MAX_X 0.10
#define DOMAIN_MAX_Y 0.10
#define DOMAIN_MAX_Z 0.03
#define CELL_SIZE_0100 0.0100
#define CELL_SIZE_0125 0.0125
#define CELL_MIN_PARTICLES_SMALL 30

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ExecuteInitialize(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Open permanent files
    mFile_Contact.open("data_contact.txt", std::ios::out);
    mFile_Contact << "1 - TIME | " << "2 - ALL | " << "3 - P-P | " << "4 - PL-PL | " << "5 - PL-PS | " << "6 - PS-PS | " << "7 - P-W | " << "8 - PL-W | " << "9 - PS-W | " << "10 - CN_ALL | " << "11 - CN_LARGE | " << "12 - CN_SMALL | " << "13 - OVERLAP_RATIO_R | " << "14 - OVERLAP_RATIO_D" << std::endl;

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
    mCellsX_0100 = static_cast<int>(std::ceil((DOMAIN_MAX_X - DOMAIN_MIN_X) / CELL_SIZE_0100));
    mCellsY_0100 = static_cast<int>(std::ceil((DOMAIN_MAX_Y - DOMAIN_MIN_Y) / CELL_SIZE_0100));
    mCellsZ_0100 = static_cast<int>(std::ceil((DOMAIN_MAX_Z - DOMAIN_MIN_Z) / CELL_SIZE_0100));
    mCellsX_0125 = static_cast<int>(std::ceil((DOMAIN_MAX_X - DOMAIN_MIN_X) / CELL_SIZE_0125));
    mCellsY_0125 = static_cast<int>(std::ceil((DOMAIN_MAX_Y - DOMAIN_MIN_Y) / CELL_SIZE_0125));
    mCellsZ_0125 = static_cast<int>(std::ceil((DOMAIN_MAX_Z - DOMAIN_MIN_Z) / CELL_SIZE_0125));
    InitializeCellsMaps_0100();
    InitializeCellsMaps_0125();
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::Calculate(ModelPart& particlesMP, ModelPart& wallsMP, bool force_execute) {
    const ProcessInfo& r_process_info = particlesMP.GetProcessInfo();
    const double time = r_process_info[TIME];
    double resid_1 = std::fmod(time, FREQUENCY_TIME_1);
    double resid_2 = std::fmod(time, FREQUENCY_TIME_2);
    bool execute_1 = !((resid_1 > FREQUENCY_TOL && std::fabs(resid_1 - FREQUENCY_TIME_1) > FREQUENCY_TOL));
    bool execute_2 = !((resid_2 > FREQUENCY_TOL && std::fabs(resid_2 - FREQUENCY_TIME_2) > FREQUENCY_TOL));
    if (force_execute || execute_1 || execute_2)
      ExecuteCalculations(particlesMP, wallsMP);
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ExecuteCalculations(ModelPart& particlesMP, ModelPart& wallsMP) {
    const ProcessInfo& r_process_info = particlesMP.GetProcessInfo();
    const double time = r_process_info[TIME];

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
    double overlap_ratio_r = 0.0;
    double overlap_ratio_d = 0.0;
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

    ResetCellsMap(mCountParticlesLarge_0100);
    ResetCellsMap(mCountParticlesSmall_0100);
    ResetCellsMap(mCountParticlesLarge_0125);
    ResetCellsMap(mCountParticlesSmall_0125);
    
    // Open temporary files
    std::ostringstream particles_filename;
    particles_filename << "data_particles_" << std::fixed << std::setprecision(6) << time << ".txt";
    std::ofstream file_particles_info(particles_filename.str(), std::ios::out);
    file_particles_info << "1 - ID | 2 - RADIUS | 3 - X | 4 - Y | 5 - Z | 6 - CONTACTS_P | 7 - CONTACTS_W | 8 - VEL_TRL_X | 9 - VEL_TRL_Y | 10 - VEL_TRL_Z | 11 - VEL_ROT_X | 12 - VEL_ROT_Y | 13 - VEL_ROT_Z | 14 - ENERGY_KINT | 15 - ENERGY_KINR | 16 - ENERGY_ELAST | 17 - ENERGY_GRAV | 18 - DISSIP_FRIC | 19 - DISSIP_DAMP" << std::endl;

    std::ostringstream walls_filename;
    walls_filename << "data_walls_" << std::fixed << std::setprecision(6) << time << ".txt";
    std::ofstream file_walls_info(walls_filename.str(), std::ios::out);
    file_walls_info << "ID X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3" << std::endl;

    std::ostringstream cells_0100_filename;
    cells_0100_filename << "data_cells_0100_" << std::fixed << std::setprecision(6) << time << ".txt";
    std::ofstream file_cells_0100_info(cells_0100_filename.str(), std::ios::out);
    file_cells_0100_info << "[XMIN XMAX YMIN YMAX ZMIN ZMAX]: #PL #PS" << std::endl;

    std::ostringstream cells_0125_filename;
    cells_0125_filename << "data_cells_0125_" << std::fixed << std::setprecision(6) << time << ".txt";
    std::ofstream file_cells_0125_info(cells_0125_filename.str(), std::ios::out);
    file_cells_0125_info << "[XMIN XMAX YMIN YMAX ZMIN ZMAX]: #PL #PS" << std::endl;
    
    // Loop over all walls
    ModelPart::ConditionsContainerType::iterator itw = wallsMP.GetCommunicator().LocalMesh().Conditions().ptr_begin();
    for (int i = 0; i < mNumWalls; i++) {
      // Write wall information to file
      DEMWall& wall = dynamic_cast<DEMWall&>(*(itw+i));
      const int id = wall.GetId();
      Condition::GeometryType &geom = wall.GetGeometry();
      const double x1 = geom[0][0]; const double y1 = geom[0][1]; const double z1 = geom[0][2];
      const double x2 = geom[1][0]; const double y2 = geom[1][1]; const double z2 = geom[1][2];
      const double x3 = geom[2][0]; const double y3 = geom[2][1]; const double z3 = geom[2][2];
      if (file_walls_info.is_open()) {
        file_walls_info << std::fixed << std::setprecision(0) << id << " " << std::fixed << std::setprecision(16) << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << " " << x3 << " " << y3 << " " << z3 << std::endl;
      }
    }

    // Loop over all particles
    ModelPart::ElementsContainerType::iterator itp = particlesMP.GetCommunicator().LocalMesh().Elements().ptr_begin();
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
      double particle_energy_kint  = 0.0;
      double particle_energy_kinr  = 0.0;
      double particle_energy_elast = 0.0;
      double particle_energy_grav  = 0.0;
      double particle_dissip_fric  = 0.0;
      double particle_dissip_damp  = 0.0;
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
        // Get neighbor properties
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
        const double overlap = (r1 + r2) - dist;
        if (overlap > 0.0) {
          // Update overlap ratios
          const double overlap_ratio_r1 = overlap / r1;
          const double overlap_ratio_d1 = overlap / (2.0 * r1);
          const double overlap_ratio_r2 = overlap / r2;
          const double overlap_ratio_d2 = overlap / (2.0 * r2);
          overlap_ratio_r += overlap_ratio_r1 + overlap_ratio_r2;
          overlap_ratio_d += overlap_ratio_d1 + overlap_ratio_d2;

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

      // Determine cells map indices of current particle and update map
      int i_x, i_y, i_z;
      
      i_x = std::max(0, std::min(mCellsX_0100 - 1, static_cast<int>(std::floor((x1 - DOMAIN_MIN_X) / CELL_SIZE_0100))));
      i_y = std::max(0, std::min(mCellsY_0100 - 1, static_cast<int>(std::floor((y1 - DOMAIN_MIN_Y) / CELL_SIZE_0100))));
      i_z = std::max(0, std::min(mCellsZ_0100 - 1, static_cast<int>(std::floor((z1 - DOMAIN_MIN_Z) / CELL_SIZE_0100))));
      if (r1 > mRadiusMean) mCountParticlesLarge_0100[i_x][i_y][i_z]++;
      else mCountParticlesSmall_0100[i_x][i_y][i_z]++;

      i_x = std::max(0, std::min(mCellsX_0125 - 1, static_cast<int>(std::floor((x1 - DOMAIN_MIN_X) / CELL_SIZE_0125))));
      i_y = std::max(0, std::min(mCellsY_0125 - 1, static_cast<int>(std::floor((y1 - DOMAIN_MIN_Y) / CELL_SIZE_0125))));
      i_z = std::max(0, std::min(mCellsZ_0125 - 1, static_cast<int>(std::floor((z1 - DOMAIN_MIN_Z) / CELL_SIZE_0125))));
      if (r1 > mRadiusMean) mCountParticlesLarge_0125[i_x][i_y][i_z]++;
      else mCountParticlesSmall_0125[i_x][i_y][i_z]++;
    }

    coord_num_all   /= mNumParticlesAll;
    coord_num_pl    /= mNumParticlesLarge;
    coord_num_ps    /= mNumParticlesSmall;
    overlap_ratio_r /= 2.0 * num_contacts_p_p;
    overlap_ratio_d /= 2.0 * num_contacts_p_p;
    vel_trl_all_avg /= mNumParticlesAll;
    vel_trl_pl_avg  /= mNumParticlesLarge;
    vel_trl_ps_avg  /= mNumParticlesSmall;
    vel_rot_all_avg /= mNumParticlesAll;
    vel_rot_pl_avg  /= mNumParticlesLarge;
    vel_rot_ps_avg  /= mNumParticlesSmall;

    // Write results
    if (file_cells_0100_info.is_open()) {
      for (int iz = 0; iz < mCellsZ_0100; ++iz) {
        double zmin = DOMAIN_MIN_Z + CELL_SIZE_0100 * iz;
        double zmax = DOMAIN_MIN_Z + CELL_SIZE_0100 * (iz + 1);
        for (int iy = 0; iy < mCellsY_0100; ++iy) {
          double ymin = DOMAIN_MIN_Y + CELL_SIZE_0100 * iy;
          double ymax = DOMAIN_MIN_Y + CELL_SIZE_0100 * (iy + 1);
          for (int ix = 0; ix < mCellsX_0100; ++ix) {
            double xmin = DOMAIN_MIN_X + CELL_SIZE_0100 * ix;
            double xmax = DOMAIN_MIN_X + CELL_SIZE_0100 * (ix + 1);
            file_cells_0100_info << std::fixed << std::setprecision(2);
            file_cells_0100_info << "[" << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << "]: ";
            file_cells_0100_info << std::fixed << std::setprecision(0);
            file_cells_0100_info << mCountParticlesLarge_0100[ix][iy][iz] << " " << mCountParticlesSmall_0100[ix][iy][iz] << std::endl;
          }
        }
      }
    }

    if (file_cells_0125_info.is_open()) {
      for (int iz = 0; iz < mCellsZ_0125; ++iz) {
        double zmin = DOMAIN_MIN_Z + CELL_SIZE_0125 * iz;
        double zmax = DOMAIN_MIN_Z + CELL_SIZE_0125 * (iz + 1);
        for (int iy = 0; iy < mCellsY_0125; ++iy) {
          double ymin = DOMAIN_MIN_Y + CELL_SIZE_0125 * iy;
          double ymax = DOMAIN_MIN_Y + CELL_SIZE_0125 * (iy + 1);
          for (int ix = 0; ix < mCellsX_0125; ++ix) {
            double xmin = DOMAIN_MIN_X + CELL_SIZE_0125 * ix;
            double xmax = DOMAIN_MIN_X + CELL_SIZE_0125 * (ix + 1);
            file_cells_0125_info << std::fixed << std::setprecision(2);
            file_cells_0125_info << "[" << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << "]: ";
            file_cells_0125_info << std::fixed << std::setprecision(0);
            file_cells_0125_info << mCountParticlesLarge_0125[ix][iy][iz] << " " << mCountParticlesSmall_0125[ix][iy][iz] << std::endl;
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
      << coord_num_all   << " "
      << coord_num_pl    << " "
      << coord_num_ps    << " "
      << overlap_ratio_r << " "
      << overlap_ratio_d << " "
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
    if (file_cells_0100_info.is_open()) file_cells_0100_info.close();
    if (file_cells_0125_info.is_open()) file_cells_0125_info.close();
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ExecuteFinalize(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Close permanent files
    if (mFile_Contact.is_open()) mFile_Contact.close();
    if (mFile_Velocity.is_open()) mFile_Velocity.close();
    if (mFile_Energy.is_open()) mFile_Energy.close();
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::InitializeCellsMaps_0100() {
    mCountParticlesLarge_0100.assign(mCellsX_0100, std::vector<std::vector<int>>(mCellsY_0100, std::vector<int>(mCellsZ_0100, 0)));
    mCountParticlesSmall_0100.assign(mCellsX_0100, std::vector<std::vector<int>>(mCellsY_0100, std::vector<int>(mCellsZ_0100, 0)));
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::InitializeCellsMaps_0125() {
    mCountParticlesLarge_0125.assign(mCellsX_0125, std::vector<std::vector<int>>(mCellsY_0125, std::vector<int>(mCellsZ_0125, 0)));
    mCountParticlesSmall_0125.assign(mCellsX_0125, std::vector<std::vector<int>>(mCellsY_0125, std::vector<int>(mCellsZ_0125, 0)));
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ResetCellsMap(std::vector<std::vector<std::vector<int>>>& map) {
    for (auto& plane : map)
      for (auto& row : plane)
        std::fill(row.begin(), row.end(), 0.0);
  }
}