//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "rve_utilities.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RVEUtilities::RVEUtilities() {}
  RVEUtilities::~RVEUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void RVEUtilities::ExecuteInitialize(ModelPart& rDEMModelPart, ModelPart& rFEMModelPart) {
    // Initialize properties
    mCompressionStage = true;
    mFrequency = 1;

    // Assemble vectors of wall elements
    AssembleWallVectors(rFEMModelPart);
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void RVEUtilities::ExecuteFinalizeSolutionStep(ModelPart& rDEMModelPart, ModelPart& rFEMModelPart) {
    ProcessInfo& r_process_info = rDEMModelPart.GetProcessInfo();
    const int time_step = r_process_info[TIME_STEPS];
    if (!mCompressionStage || time_step % mFrequency != 0.0)
      return;

    // Initialize variables
    double volume_solid = 0.0;

    // Compute total volume
    ComputeTotalVolume();
    
    // Loop over all particles
    const int num_of_particles = rDEMModelPart.NumberOfElements();
    ModelPart::ElementsContainerType::iterator it = rDEMModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();
    #pragma omp parallel
    for (int i = 0; i < num_of_particles; i++) {
      SphericParticle& particle = dynamic_cast<SphericParticle&> (*(it+i));

      // Loop over neighbors with higher ID number
      for (unsigned int j = 0; j < particle.mNeighbourElements.size(); j++) {
        if (particle.mNeighbourElements[i] == NULL || particle.GetId() > particle.mNeighbourElements[j]->GetId())
          continue;

      }
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void RVEUtilities::AssembleWallVectors(ModelPart& rFEMModelPart) {
    const double eps = std::numeric_limits<double>::epsilon();

    std::vector<DEMWall*> wall_elems_x;
    std::vector<DEMWall*> wall_elems_y;
    std::vector<DEMWall*> wall_elems_z;

    double xmin =  DBL_MAX;
    double xmax = -DBL_MAX;
    double ymin =  DBL_MAX;
    double ymax = -DBL_MAX;
    double zmin =  DBL_MAX;
    double zmax = -DBL_MAX;

    ModelPart::ConditionsContainerType& r_conditions = rFEMModelPart.GetCommunicator().LocalMesh().Conditions();
    for (int i = 0; i < r_conditions.size(); i++) {
      ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
      DEMWall* p_wall = dynamic_cast<DEMWall*> (&(*it));

      Condition::GeometryType& geom = p_wall->GetGeometry();
      const unsigned int& dim = geom.WorkingSpaceDimension();

      if (dim == 2) {
        auto& node1 = geom[0];
        auto& node2 = geom[1];
        double coords1[2] = {node1[0], node1[1]};
        double coords2[2] = {node2[0], node2[1]};

        if (std::abs(coords1[0] - coords2[0]) < eps) {
          wall_elems_x.push_back(p_wall);
          if (coords1[0] < xmin) xmin = coords1[0];
          if (coords1[0] > xmax) xmax = coords1[0];
        }
        else if (std::abs(coords1[1] - coords2[1]) < eps) {
          wall_elems_y.push_back(p_wall);
          if (coords1[1] < ymin) ymin = coords1[1];
          if (coords1[1] > ymax) ymax = coords1[1];
        }
      }

      else if (dim == 3) {
        auto& node1 = geom[0];
        auto& node2 = geom[1];
        auto& node3 = geom[2];
        double coords1[3] = {node1[0], node1[1], node1[2]};
        double coords2[3] = {node2[0], node2[1], node2[2]};
        double coords3[3] = {node3[0], node3[1], node3[2]};

        if (std::abs(coords1[0] - coords2[0]) < eps &&
            std::abs(coords1[0] - coords3[0]) < eps) {
          wall_elems_x.push_back(p_wall);
          if (coords1[0] < xmin) xmin = coords1[0];
          if (coords1[0] > xmax) xmax = coords1[0];
        }
        else if (std::abs(coords1[1] - coords2[1]) < eps &&
                 std::abs(coords1[1] - coords3[1]) < eps) {
          wall_elems_y.push_back(p_wall);
          if (coords1[1] < ymin) ymin = coords1[1];
          if (coords1[1] > ymax) ymax = coords1[1];
        }
        else if (std::abs(coords1[2] - coords2[2]) < eps &&
                 std::abs(coords1[2] - coords3[2]) < eps) {
          wall_elems_z.push_back(p_wall);
          if (coords1[2] < zmin) zmin = coords1[2];
          if (coords1[2] > zmax) zmax = coords1[2];
        }
      }
    }

    for (int i = 0; i < wall_elems_x.size(); i++) {
      const double x = wall_elems_x[i]->GetGeometry()[0][0];
      if      (std::abs(x - xmin) < eps) mWallXMin.push_back(wall_elems_x[i]);
      else if (std::abs(x - xmax) < eps) mWallXMax.push_back(wall_elems_x[i]);
    }
    for (int i = 0; i < wall_elems_y.size(); i++) {
      const double y = wall_elems_y[i]->GetGeometry()[0][1];
      if      (std::abs(y - ymin) < eps) mWallYMin.push_back(wall_elems_y[i]);
      else if (std::abs(y - ymax) < eps) mWallYMax.push_back(wall_elems_y[i]);
    }
    for (int i = 0; i < wall_elems_z.size(); i++) {
      const double z = wall_elems_z[i]->GetGeometry()[0][2];
      if      (std::abs(z - zmin) < eps) mWallZMin.push_back(wall_elems_z[i]);
      else if (std::abs(z - zmax) < eps) mWallZMax.push_back(wall_elems_z[i]);
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void RVEUtilities::ComputeTotalVolume(void) {
    double dX = 1.0;
    double dY = 1.0;
    double dZ = 1.0;

    if (mWallXMax.size() > 0 && mWallXMin.size() > 0) dX = std::abs(mWallXMax[0]->GetGeometry()[0][0] - mWallXMin[0]->GetGeometry()[0][0]);
    if (mWallYMax.size() > 0 && mWallYMin.size() > 0) dY = std::abs(mWallYMax[0]->GetGeometry()[0][1] - mWallYMin[0]->GetGeometry()[0][1]);
    if (mWallZMax.size() > 0 && mWallZMin.size() > 0) dZ = std::abs(mWallZMax[0]->GetGeometry()[0][2] - mWallZMin[0]->GetGeometry()[0][2]);

    mVolume = dX * dY * dZ;
  }

} // namespace Kratos
