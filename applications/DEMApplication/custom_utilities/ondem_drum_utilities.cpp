#include "ondem_drum_utilities.h"

#define FREQUENCY_PRINT_LMI 0.010
#define FREQUENCY_PRINT_EMI 0.25
#define FREQUENCY_TOL 1e-9
#define DOMAIN_MIN_X -0.10
#define DOMAIN_MIN_Y -0.10
#define DOMAIN_MIN_Z -0.03
#define DOMAIN_MAX_X 0.10
#define DOMAIN_MAX_Y 0.10
#define DOMAIN_MAX_Z 0.03
#define CELL_SIZE_LMI 0.010
#define CELL_SIZE_EMI 0.125
#define MIN_NUM_PARTICLES_LMI 30

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ExecuteInitialize(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Open permanent files
    mFile_Contact.open("data_contact.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(mFile_Contact) << "Could not open data_contact.txt!" << std::endl;
    mFile_Contact << "1 - TIME | ";
    mFile_Contact << "2 - ALL | ";
    mFile_Contact << "3 - P-P | ";
    mFile_Contact << "4 - PLARGE-PLARGE | ";
    mFile_Contact << "5 - PSMALL-PSMALL | ";
    mFile_Contact << "6 - PLARGE-PSMALL | ";
    mFile_Contact << "7 - P-W | ";
    mFile_Contact << "8 - PLARGE-W | ";
    mFile_Contact << "9 - PSMALL-W | ";
    mFile_Contact << "10 - CN_ALL | ";
    mFile_Contact << "11 - CN_LARGE | ";
    mFile_Contact << "12 - CN_SMALL";
    mFile_Contact << std::endl;

    mFile_Dissipation.open("data_dissipation.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(mFile_Dissipation) << "Could not open data_dissipation.txt!" << std::endl;
    mFile_Dissipation << "1 - TIME | ";
    mFile_Dissipation << "2 - TOTAL | ";
    mFile_Dissipation << "3 - FRICTION | ";
    mFile_Dissipation << "4 - DAMPING";
    mFile_Dissipation << std::endl;

    mFile_Velocity.open("data_velocity.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(mFile_Velocity) << "Could not open data_velocity.txt!" << std::endl;
    mFile_Velocity << "1 - TIME | ";
    mFile_Velocity << "2 - AVG TRANSLATION (ALL) | ";
    mFile_Velocity << "3 - AVG TRANSLATION (LARGE) | ";
    mFile_Velocity << "4 - AVG TRANSLATION (SMALL) | ";
    mFile_Velocity << "5 - AVG ROTATION (ALL) | ";
    mFile_Velocity << "6 - AVG ROTATION (LARGE) | ";
    mFile_Velocity << "7 - AVG ROTATION (SMALL)";
    mFile_Velocity << std::endl;

    mFile_CellsLMI.open("data_cells_lmi.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(mFile_CellsLMI) << "Could not open data_cells_lmi.txt!" << std::endl;
    mFile_CellsLMI << "#TIME" << std::endl;
    mFile_CellsLMI << "[XMIN,XMAX; YMIN,YMAX; ZMIN,ZMAX]: (NUM_PARTICLES_LARGE, NUM_PARTICLES_SMALL)" << std::endl << std::endl;

    mFile_LMI.open("data_lacey.txt", std::ios::out);
    KRATOS_ERROR_IF_NOT(mFile_LMI) << "Could not open data_lacey.txt!" << std::endl;
    mFile_LMI << "1 - TIME | ";
    mFile_LMI << "2 - LMI";
    mFile_LMI << std::endl;

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
    mCellsX = static_cast<int>((DOMAIN_MAX_X - DOMAIN_MIN_X) / CELL_SIZE_LMI);
    mCellsY = static_cast<int>((DOMAIN_MAX_Y - DOMAIN_MIN_Y) / CELL_SIZE_LMI);
    mCellsZ = static_cast<int>((DOMAIN_MAX_Z - DOMAIN_MIN_Z) / CELL_SIZE_LMI);
    InitializeMapLMI(mCountParticlesLarge);
    InitializeMapLMI(mCountParticlesSmall);
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::Calculate(ModelPart& particlesMP, ModelPart& wallsMP) {
    // Check if it is time to execute
    const ProcessInfo& r_process_info = particlesMP.GetProcessInfo();
    const double time = r_process_info[TIME];
    double resid_LMI = std::fmod(time, FREQUENCY_PRINT_LMI);
    if (resid_LMI > FREQUENCY_TOL && std::fabs(resid_LMI - FREQUENCY_PRINT_LMI) > FREQUENCY_TOL) return;

    // Initialize variables
    int num_contacts_all      = 0;
    int num_contacts_p_p      = 0;
    int num_contacts_pl_pl    = 0;
    int num_contacts_ps_ps    = 0;
    int num_contacts_pl_ps    = 0;
    int num_contacts_p_w      = 0;
    int num_contacts_pl_w     = 0;
    int num_contacts_ps_w     = 0;
    double coord_num_all      = 0.0;
    double coord_num_large    = 0.0;
    double coord_num_small    = 0.0;
    double vel_trl_all_avg    = 0.0;
    double vel_trl_large_avg  = 0.0;
    double vel_trl_small_avg  = 0.0;
    double vel_rot_all_avg    = 0.0;
    double vel_rot_large_avg  = 0.0;
    double vel_rot_small_avg  = 0.0;
    double energy_dissip_fric = 0.0;
    double energy_dissip_damp = 0.0;
    double energy_dissip      = 0.0;
    double lmi                = 0.0;
    ResetMapLMI(mCountParticlesLarge);
    ResetMapLMI(mCountParticlesSmall);

    // Open temporary files
    std::ofstream file_walls_info((std::ostringstream() << "data_walls_info_" << std::fixed << std::setprecision(6) << time << ".txt").str());
    KRATOS_ERROR_IF_NOT(file_walls_info) << "Could not open file data_walls_info_" << std::fixed << std::setprecision(6) << time << ".txt!" << std::endl;
    file_walls_info << "ID X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3";
    file_walls_info << std::endl;

    std::ofstream file_particles_info((std::ostringstream() << "data_particles_info_" << std::fixed << std::setprecision(6) << time << ".txt").str());
    KRATOS_ERROR_IF_NOT(file_particles_info) << "Could not open file data_particles_info_" << std::fixed << std::setprecision(6) << time << ".txt!" << std::endl;
    file_particles_info << "1 - ID | 2 - RADIUS | 3 - X | 4 - Y | 5 - Z | 6 - CONTACTS_P | 7 - CONTACTS_W | 8 - VEL_TRL_X | 9 - VEL_TRL_Y | 10 - VEL_TRL_Z | 11 - VEL_ROT_X | 12 - VEL_ROT_Y | 13 - VEL_ROT_Z | 14 - DISSIP_FRIC | 15 - DISSIP_DAMP";
    file_particles_info << std::endl;

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
      const int num_neighbors_particles = particle.mNeighbourElements.size();
      const int num_neighbors_walls = particle.mNeighbourRigidFaces.size();
      const int num_neighbors = num_neighbors_particles + num_neighbors_walls;
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
      double particle_dissip_fric = 0.0;
      double particle_dissip_damp = 0.0;
      particle.Calculate(PARTICLE_INELASTIC_FRICTIONAL_ENERGY, particle_dissip_fric, r_process_info);
      particle.Calculate(PARTICLE_INELASTIC_VISCODAMPING_ENERGY, particle_dissip_damp, r_process_info);

      // Write particle information to file
      if (file_particles_info.is_open()) {
        file_particles_info
        << std::fixed << std::setprecision(0)
        << id1 << " "
        << std::fixed << std::setprecision(16)
        << r1 << " " << x1 << " " << y1 << " " << z1 << " "
        << std::fixed << std::setprecision(0)
        << num_neighbors_particles << " " << num_neighbors_walls << " "
        << std::fixed << std::setprecision(16)
        << vel_trl_x << " " << vel_trl_y << " " << vel_trl_z << " "
        << vel_rot_x << " " << vel_rot_y << " " << vel_rot_z << " "
        << particle_dissip_fric << " " << particle_dissip_damp
        << std::endl;
      }

      // Update contact statistics
      num_contacts_all += num_neighbors_walls;
      num_contacts_p_w += num_neighbors_walls;
      coord_num_all += num_neighbors;
      if (r1 > mRadiusMean) {
        num_contacts_pl_w += num_neighbors_walls;
        coord_num_large += num_neighbors;
      }
      else {
        num_contacts_ps_w += num_neighbors_walls;
        coord_num_small += num_neighbors;
      }

      // Loop over particle neighbors
      if (true) { // Faster: evaluate existing interactions
        for (int j = 0; j < num_neighbors_particles; j++) {
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
      }
      else { // Slower: evaluate all interactions (each interaction only once)
        for (int j = i+1; j < mNumParticlesAll; j++) {
          // Get particle properties
          SphericParticle& neighbor = dynamic_cast<SphericParticle&>(*(itp+j));
          const auto& node2 = neighbor.GetGeometry()[0];
          const double r2 = neighbor.GetRadius();
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
      }
      
      // Accumulate velocities
      vel_trl_all_avg += vel_trl_norm;
      vel_rot_all_avg += vel_rot_norm;
      if (r1 > mRadiusMean) {
        vel_trl_large_avg += vel_trl_norm;
        vel_rot_large_avg += vel_rot_norm;
      }
      else {
        vel_trl_small_avg += vel_trl_norm;
        vel_rot_small_avg += vel_rot_norm;
      }

      // Accumulate energy dissipation
      energy_dissip_fric += particle_dissip_fric;
      energy_dissip_damp += particle_dissip_damp;
      energy_dissip      += particle_dissip_fric + particle_dissip_damp;

      // Determine map indices of current particle
      const int i_x = std::max(0, std::min(mCellsX - 1, static_cast<int>(std::floor((x1 - DOMAIN_MIN_X) / CELL_SIZE_LMI))));
      const int i_y = std::max(0, std::min(mCellsY - 1, static_cast<int>(std::floor((y1 - DOMAIN_MIN_Y) / CELL_SIZE_LMI))));
      const int i_z = std::max(0, std::min(mCellsZ - 1, static_cast<int>(std::floor((z1 - DOMAIN_MIN_Z) / CELL_SIZE_LMI))));

      // Update maps
      if (r1 > mRadiusMean) mCountParticlesLarge[i_x][i_y][i_z]++;
      else mCountParticlesSmall[i_x][i_y][i_z]++;
    }
    
    // Compute Lacey Mixing Index
    std::vector<double> x_list;
    std::vector<double> n_list;
    for (int ix = 0; ix < mCellsX; ++ix) {
      for (int iy = 0; iy < mCellsY; ++iy) {
        for (int iz = 0; iz < mCellsZ; ++iz) {
          int n_small = mCountParticlesSmall[ix][iy][iz];
          int n_large = mCountParticlesLarge[ix][iy][iz];
          int n_tot   = n_small + n_large;
          if (n_small < MIN_NUM_PARTICLES_LMI || n_tot == 0) continue;
          double frac_small = static_cast<double>(n_small) / n_tot;
          x_list.push_back(frac_small);
          n_list.push_back(n_tot);
        }
      }
    }
    const size_t N = x_list.size();
    if (N > 0) {
      // Mean fraction of small particles
      double x_bar = std::accumulate(x_list.begin(), x_list.end(), 0.0) / static_cast<double>(N);
      // Observed variance
      double sigma_obs2 = 0.0;
      for (auto xi : x_list) sigma_obs2 += (xi - x_bar) * (xi - x_bar);
      // Population variance
      sigma_obs2 /= static_cast<double>(N);
      // Mean particles per cell
      double n_bar = std::accumulate(n_list.begin(), n_list.end(), 0.0) / static_cast<double>(N);
      // Random mixing variance (binomial assumption)
      double sigma_rand2 = 0.0;
      if (n_bar > 0.0) sigma_rand2 = x_bar * (1.0 - x_bar) / n_bar;
      // Compute LMI
      if (sigma_rand2 > 0.0) lmi = 1.0 - (sigma_obs2 / sigma_rand2);
    }

    coord_num_all     /= mNumParticlesAll;
    coord_num_large   /= mNumParticlesLarge;
    coord_num_small   /= mNumParticlesSmall;
    vel_trl_all_avg   /= mNumParticlesAll;
    vel_trl_large_avg /= mNumParticlesLarge;
    vel_trl_small_avg /= mNumParticlesSmall;
    vel_rot_all_avg   /= mNumParticlesAll;
    vel_rot_large_avg /= mNumParticlesLarge;
    vel_rot_small_avg /= mNumParticlesSmall;

    // Write results
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
      << coord_num_all      << " "
      << coord_num_large    << " "
      << coord_num_small    << " "
      << std::endl;
    }
    if (mFile_Velocity.is_open()) {
      mFile_Velocity 
      << std::fixed << std::setprecision(6)
      << time << " "
      << std::fixed << std::setprecision(15)
      << vel_trl_all_avg   << " "
      << vel_trl_large_avg << " "
      << vel_trl_small_avg << " "
      << vel_rot_all_avg   << " "
      << vel_rot_large_avg << " "
      << vel_rot_small_avg << " "
      << std::endl;
    }
    if (mFile_Dissipation.is_open()) {
      mFile_Dissipation
      << std::fixed << std::setprecision(6)
      << time << " "
      << std::fixed << std::setprecision(15)
      << energy_dissip      << " "
      << energy_dissip_fric << " "
      << energy_dissip_damp << " "
      << std::endl;
    }
    if (mFile_CellsLMI.is_open()) {
      mFile_CellsLMI << std::fixed << std::setprecision(6) << "#TIME: " << time << std::endl;
      for (int iz = 0; iz < mCellsZ; ++iz) {
        double zmin = DOMAIN_MIN_Z + CELL_SIZE_LMI * iz;
        double zmax = DOMAIN_MIN_Z + CELL_SIZE_LMI * (iz + 1);
        for (int iy = 0; iy < mCellsY; ++iy) {
          double ymin = DOMAIN_MIN_Y + CELL_SIZE_LMI * iy;
          double ymax = DOMAIN_MIN_Y + CELL_SIZE_LMI * (iy + 1);
          for (int ix = 0; ix < mCellsX; ++ix) {
            double xmin = DOMAIN_MIN_X + CELL_SIZE_LMI * ix;
            double xmax = DOMAIN_MIN_X + CELL_SIZE_LMI * (ix + 1);
            mFile_CellsLMI << std::fixed << std::setprecision(2);
            mFile_CellsLMI << "[" << xmin << "," << xmax << "; " << ymin << "," << ymax << "; "  << zmin << "," << zmax << "]: ";
            mFile_CellsLMI << std::fixed << std::setprecision(0);
            mFile_CellsLMI << "(" << mCountParticlesLarge[ix][iy][iz] << ", " << mCountParticlesSmall[ix][iy][iz] << ") " << std::endl;
          }
        }
      }
      mFile_CellsLMI << std::endl;
    }
    if (mFile_LMI.is_open()) {
      mFile_LMI
      << std::fixed << std::setprecision(6)
      << time << " "
      << std::fixed << std::setprecision(15)
      << lmi << " "
      << std::endl;
    }

    // Close temporary files
    if (file_particles_info.is_open()) file_particles_info.close();
    if (file_walls_info.is_open()) file_walls_info.close();
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ExecuteFinalize(ModelPart& particlesMP, ModelPart& wallsMP) {
    if (mFile_Contact.is_open()) mFile_Contact.close();
    if (mFile_Velocity.is_open()) mFile_Velocity.close();
    if (mFile_Dissipation.is_open()) mFile_Dissipation.close();
    if (mFile_CellsLMI.is_open()) mFile_CellsLMI.close();
    if (mFile_LMI.is_open()) mFile_LMI.close();
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::InitializeMapLMI(std::vector<std::vector<std::vector<int>>>& map) {
    map.assign(mCellsX, std::vector<std::vector<int>>(mCellsY, std::vector<int>(mCellsZ, 0)));
  }

  //------------------------------------------------------------------------------------------------------------
  void ONDEMDrumUtilities::ResetMapLMI(std::vector<std::vector<std::vector<int>>>& map) {
    for (auto& plane : map)
      for (auto& row : plane)
        std::fill(row.begin(), row.end(), 0.0);
  }
}