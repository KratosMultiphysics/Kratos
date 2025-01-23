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
#include "tesselation_utilities_3d.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  TesselationUtilities3D::TesselationUtilities3D() {
    mUpdateVoronoi = false;
    mUpdatePorosiy = false;
    mAlphaRadius   = 0.0;
    mSwitches      = "";
  }

  TesselationUtilities3D::~TesselationUtilities3D() {}

  //------------------------------------------------------------------------------------------------------------
  void TesselationUtilities3D::ExecuteInitialize(ModelPart& rModelPart, bool update_voronoi, bool update_porosity) {
    KRATOS_TRY

    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    // Set flags
    mUpdateVoronoi = update_voronoi;
    mUpdatePorosiy = update_porosity;

    // Initialize average porosity
    if (update_porosity)
      r_process_info[AVERAGE_POROSITY] = 1.0;

    // Set switches for TetGen
    if      (update_voronoi)  mSwitches = "JQv";
    else if (update_porosity) mSwitches = "JQ";

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void TesselationUtilities3D::ExecuteInitializeSolutionStep(ModelPart& rModelPart) {
    Tetrahedralization(rModelPart);
  }

  //------------------------------------------------------------------------------------------------------------
  void TesselationUtilities3D::Tetrahedralization(ModelPart& rModelPart) {
    KRATOS_TRY

    const int num_particles = rModelPart.NumberOfElements();
    if (num_particles < 4)
      return;

    // Build input
    struct tetgenio in, out;
    in.firstnumber    = 1;
    in.mesh_dim       = 3;
    in.numberofpoints = num_particles;
    in.pointlist      = (double*)malloc(sizeof(double) * in.numberofpoints * 3);

    ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin();
    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < in.numberofpoints; i++) {
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*(it+i));

      particle.mDelaunayPointListIndex = i;

      const array_1d<double, 3>& coors = particle.GetParticleCoordinates();
      in.pointlist[3 * i + 0] = coors[0];
      in.pointlist[3 * i + 1] = coors[1];
      in.pointlist[3 * i + 2] = coors[2];

      if (mUpdateVoronoi)
        particle.mNeighborVoronoiRadius.clear();
    }

    // Perform tetrahedralization
    int fail = 0;
    try {
      tetrahedralize(&mSwitches[0], &in, &out);
    }
    catch (int error_code) {
      fail = error_code;
    }

    if (fail || out.numberoftetrahedra == 0 || in.numberofpoints != out.numberofpoints) {
      KRATOS_ERROR_IF(rModelPart.GetProcessInfo()[TIME_STEPS] == 1) << "Fail to generate tetrahedralization!" << std::endl;
      KRATOS_WARNING("DEM") << std::endl;
      KRATOS_WARNING("DEM") << "Fail to generate tetrahedralization! Results from previous successful tetrahedralization will be used." << std::endl;
      KRATOS_WARNING("DEM") << std::endl;
      return;
    }

    // Update voronoi diagram
    if (mUpdateVoronoi)
      UpdateVoronoi(rModelPart, out);

    // Update average porosity
    if (mUpdatePorosiy)
      UpdatePorosity(rModelPart, out);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  /*
  * Build a table for each particle with needed information from 3D voronoi diagram:
  * column 1: neighbor particle IDs
  * column 2: voronoi face equivalent radius (for a circle with the same area)
  */
  void TesselationUtilities3D::UpdateVoronoi(ModelPart& rModelPart, struct tetgenio& rOut) {
    ModelPart::ElementsContainerType& rElements = rModelPart.GetCommunicator().LocalMesh().Elements();

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < rOut.numberofvfacets; i++) {
      // Compute area of voronoi face by irradiation of triangles
      tetgenio::vorofacet face = rOut.vfacetlist[i];
      const int num_face_edges = face.elist[0];
      double face_area         = 0.0;

      // Vertices of 1st edge of voronoi face
      tetgenio::voroedge edge = rOut.vedgelist[face.elist[1] - 1];
      int ev1 = edge.v1 - 1;
      int ev2 = edge.v2 - 1;

      // Check if face is bounded
      bool bounded_face = (ev1 >= 0 && ev2 >= 0);

      if (bounded_face) {
        // Irradiating coordinates from a vertex of 1st edge
        const double x0 = rOut.vpointlist[3 * ev1 + 0];
        const double y0 = rOut.vpointlist[3 * ev1 + 1];
        const double z0 = rOut.vpointlist[3 * ev1 + 2];

        // Triangle area irradiation
        for (int j = 2; j <= num_face_edges; j++) {
          // Edge vertices
          edge = rOut.vedgelist[face.elist[j] - 1];
          ev1 = edge.v1 - 1;
          ev2 = edge.v2 - 1;

          // Check if face is bounded
          if (ev1 < 0 || ev2 < 0) {
            bounded_face = false;
            break;
          }

          // Edge vertices coordinates
          const double x1 = rOut.vpointlist[3 * ev1 + 0];
          const double y1 = rOut.vpointlist[3 * ev1 + 1];
          const double z1 = rOut.vpointlist[3 * ev1 + 2];
          const double x2 = rOut.vpointlist[3 * ev2 + 0];
          const double y2 = rOut.vpointlist[3 * ev2 + 1];
          const double z2 = rOut.vpointlist[3 * ev2 + 2];

          // Area of triangle in space
          array_1d<double, 3> AB;
          array_1d<double, 3> AC;
          array_1d<double, 3> ABxAC;
          AB[0] = x1 - x0;
          AB[1] = y1 - y0;
          AB[2] = z1 - z0;
          AC[0] = x2 - x0;
          AC[1] = y2 - y0;
          AC[2] = z2 - z0;
          GeometryFunctions::CrossProduct(AB, AC, ABxAC);
          face_area += DEM_MODULUS_3(ABxAC) / 2.0;
        }
      }

      // Delaunay vertices adjacent to voronoi face
      const int vd1 = face.c1 - 1;
      const int vd2 = face.c2 - 1;

      // Particles corresponding to delaunay vertices
      ModelPart::ElementsContainerType::iterator it1 = rElements.ptr_begin() + vd1;
      ModelPart::ElementsContainerType::iterator it2 = rElements.ptr_begin() + vd2;
      ThermalSphericParticle& particle_1 = dynamic_cast<ThermalSphericParticle&> (*it1);
      ThermalSphericParticle& particle_2 = dynamic_cast<ThermalSphericParticle&> (*it2);

      // Compute equivalent radius
      double radius = 0.0;

      if (bounded_face)
        radius = sqrt(face_area / Globals::Pi);
      else
        // Assumption: radius of unbounded voronoi face is proportional to the particles radii
        radius = 0.75 * (particle_1.GetRadius() + particle_2.GetRadius());

      // Add info to table of both particles
      particle_1.mNeighborVoronoiRadius[vd2] = radius;
      particle_2.mNeighborVoronoiRadius[vd1] = radius;
    }
  }

  //------------------------------------------------------------------------------------------------------------
  /*
  * Compute total volume of delaunay tetahedra and particles to obtain global avarege porosity.
  * Oversized tetahedra may be removed using the alpha-shape method.
  */
  void TesselationUtilities3D::UpdatePorosity(ModelPart& rModelPart, struct tetgenio& rOut) {
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    std::vector<int> addedParticle(rOut.numberofpoints,0);
    double total_volume    = 0.0;
    double particle_volume = 0.0;

    // Compute alpha radius
    // Assumption: Only in the 1st step - keeping alpha radius constant throghout simulation
    if (r_process_info[POROSITY_METHOD_NAME].compare("average_alpha_shape") == 0 && r_process_info[TIME_STEPS] == 1)
      ComputeAlphaRadius(rModelPart, rOut);

    // Accumulate total volume of tetahedra and particles
    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < rOut.numberoftetrahedra; i++) {
      // Get vertices IDs
      const int v1 = rOut.tetrahedronlist[4 * i + 0] - 1;
      const int v2 = rOut.tetrahedronlist[4 * i + 1] - 1;
      const int v3 = rOut.tetrahedronlist[4 * i + 2] - 1;
      const int v4 = rOut.tetrahedronlist[4 * i + 3] - 1;

      // Get vertices coordinates
      const double x1 = rOut.pointlist[3 * v1 + 0];
      const double y1 = rOut.pointlist[3 * v1 + 1];
      const double z1 = rOut.pointlist[3 * v1 + 2];
      const double x2 = rOut.pointlist[3 * v2 + 0];
      const double y2 = rOut.pointlist[3 * v2 + 1];
      const double z2 = rOut.pointlist[3 * v2 + 2];
      const double x3 = rOut.pointlist[3 * v3 + 0];
      const double y3 = rOut.pointlist[3 * v3 + 1];
      const double z3 = rOut.pointlist[3 * v3 + 2];
      const double x4 = rOut.pointlist[3 * v4 + 0];
      const double y4 = rOut.pointlist[3 * v4 + 1];
      const double z4 = rOut.pointlist[3 * v4 + 2];

      // Perform alpha-shape
      if (r_process_info[POROSITY_METHOD_NAME].compare("average_alpha_shape") == 0) {
        std::vector<double> coords = { x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 };
        if (!AlphaShape(coords)) {
          continue;
        }
      }

      // Add tetahedron volume
      array_1d<double, 3> a;
      array_1d<double, 3> b;
      array_1d<double, 3> c;
      array_1d<double, 3> bxc;
      a[0] = x2 - x1;
      a[1] = y2 - y1;
      a[2] = z2 - z1;
      b[0] = x3 - x1;
      b[1] = y3 - y1;
      b[2] = z3 - z1;
      c[0] = x4 - x1;
      c[1] = y4 - y1;
      c[2] = z4 - z1;
      GeometryFunctions::CrossProduct(b, c, bxc);
      total_volume += std::abs(GeometryFunctions::DotProduct(a, bxc)) / 6.0;

      // Add particles volume
      AddParticleVolume(rModelPart, addedParticle, particle_volume, v1);
      AddParticleVolume(rModelPart, addedParticle, particle_volume, v2);
      AddParticleVolume(rModelPart, addedParticle, particle_volume, v3);
      AddParticleVolume(rModelPart, addedParticle, particle_volume, v4);
    }

    // Set new average porosity
    r_process_info[AVERAGE_POROSITY] = 1.0 - particle_volume / total_volume;
  }

  //------------------------------------------------------------------------------------------------------------
  /*
  * Compute mean mesh size for alpha-shape.
  * The mean mesh size is taken as the average of the smallest edge of each tetahedron.
  */
  void TesselationUtilities3D::ComputeAlphaRadius(ModelPart& rModelPart, struct tetgenio& rOut) {
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    double MeanMeshSize = 0.0;

    for (int i = 0; i < rOut.numberoftetrahedra; i++) {
      // Get vertices IDs
      const int v1 = rOut.tetrahedronlist[4 * i + 0] - 1;
      const int v2 = rOut.tetrahedronlist[4 * i + 1] - 1;
      const int v3 = rOut.tetrahedronlist[4 * i + 2] - 1;
      const int v4 = rOut.tetrahedronlist[4 * i + 3] - 1;
      
      // Get vertices coordinates
      const double x1 = rOut.pointlist[3 * v1 + 0];
      const double y1 = rOut.pointlist[3 * v1 + 1];
      const double z1 = rOut.pointlist[3 * v1 + 2];
      const double x2 = rOut.pointlist[3 * v2 + 0];
      const double y2 = rOut.pointlist[3 * v2 + 1];
      const double z2 = rOut.pointlist[3 * v2 + 2];
      const double x3 = rOut.pointlist[3 * v3 + 0];
      const double y3 = rOut.pointlist[3 * v3 + 1];
      const double z3 = rOut.pointlist[3 * v3 + 2];
      const double x4 = rOut.pointlist[3 * v4 + 0];
      const double y4 = rOut.pointlist[3 * v4 + 1];
      const double z4 = rOut.pointlist[3 * v4 + 2];
      
      // Get minimum edge length
      std::vector<double>len;
      len.push_back(sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2)));
      len.push_back(sqrt(pow(x3 - x1, 2) + pow(y3 - y1, 2) + pow(z3 - z1, 2)));
      len.push_back(sqrt(pow(x4 - x1, 2) + pow(y4 - y1, 2) + pow(z4 - z1, 2)));
      len.push_back(sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2) + pow(z3 - z2, 2)));
      len.push_back(sqrt(pow(x4 - x2, 2) + pow(y4 - y2, 2) + pow(z4 - z2, 2)));
      len.push_back(sqrt(pow(x4 - x3, 2) + pow(y4 - y3, 2) + pow(z4 - z3, 2)));
      
      MeanMeshSize += *std::min_element(len.begin(), len.end());
    }

    // Average minimum length
    MeanMeshSize /= rOut.numberoftetrahedra;

    // Alpha radius
    mAlphaRadius = MeanMeshSize * r_process_info[ALPHA_SHAPE_PARAMETER];
  }

  //------------------------------------------------------------------------------------------------------------
  /*
  * Perform alpha-shape verification on a delaunay tetahedron to remove distorted shapes.
  */
  bool TesselationUtilities3D::AlphaShape(std::vector<double>& coords) {
    const double x1 = coords[0];
    const double y1 = coords[1];
    const double z1 = coords[2];
    const double x2 = coords[3];
    const double y2 = coords[4];
    const double z2 = coords[5];
    const double x3 = coords[6];
    const double y3 = coords[7];
    const double z3 = coords[8];
    const double x4 = coords[9];
    const double y4 = coords[10];
    const double z4 = coords[11];

    // Calculate Jacobian
    BoundedMatrix<double, 3, 3> J;
    J(0, 0) = x2 - x1;
    J(1, 0) = x3 - x1;
    J(2, 0) = x4 - x1;
    J(0, 1) = y2 - y1;
    J(1, 1) = y3 - y1;
    J(2, 1) = y4 - y1;
    J(0, 2) = z2 - z1;
    J(1, 2) = z3 - z1;
    J(2, 2) = z4 - z1;

    // Calculate the inverse of the Jacobian
    BoundedMatrix<double, 3, 3> Jinv;
    Jinv(0, 0) =  J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1);
    Jinv(1, 0) = -J(1, 0) * J(2, 2) + J(1, 2) * J(2, 0);
    Jinv(2, 0) =  J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0);
    Jinv(0, 1) = -J(0, 1) * J(2, 2) + J(0, 2) * J(2, 1);
    Jinv(1, 1) =  J(0, 0) * J(2, 2) - J(0, 2) * J(2, 0);
    Jinv(2, 1) = -J(0, 0) * J(2, 1) + J(0, 1) * J(2, 0);
    Jinv(0, 2) =  J(0, 1) * J(1, 2) - J(0, 2) * J(1, 1);
    Jinv(1, 2) = -J(0, 0) * J(1, 2) + J(0, 2) * J(1, 0);
    Jinv(2, 2) =  J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);

    // Calculate the determinant (volume/6)
    const double vol = J(0, 0) * Jinv(0, 0) + J(0, 1) * Jinv(1, 0) + J(0, 2) * Jinv(2, 0);

    // Calculate sphere center
    BoundedVector<double, 3> RHS;
    Vector Center = ZeroVector(3);
    RHS[0]    =  (J(0, 0)  * J(0, 0)  + J(0, 1) * J(0, 1) + J(0, 2) * J(0, 2));
    RHS[1]    =  (J(1, 0)  * J(1, 0)  + J(1, 1) * J(1, 1) + J(1, 2) * J(1, 2));
    RHS[2]    =  (J(2, 0)  * J(2, 0)  + J(2, 1) * J(2, 1) + J(2, 2) * J(2, 2));
    Center[0] =  (+RHS[0]) * (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1));
    Center[1] =  (-RHS[0]) * (J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0));
    Center[2] =  (+RHS[0]) * (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0));
    Center[0] += (+RHS[1]) * (J(2, 1) * J(0, 2) - J(2, 2) * J(0, 1));
    Center[1] += (-RHS[1]) * (J(2, 0) * J(0, 2) - J(2, 2) * J(0, 0));
    Center[2] += (+RHS[1]) * (J(2, 0) * J(0, 1) - J(2, 1) * J(0, 0));
    Center[0] += (+RHS[2]) * (J(0, 1) * J(1, 2) - J(0, 2) * J(1, 1));
    Center[1] += (-RHS[2]) * (J(0, 0) * J(1, 2) - J(0, 2) * J(1, 0));
    Center[2] += (+RHS[2]) * (J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0));
    Center /= (2.0 * vol);

    // Calculate sphere radius
    const double radius = norm_2(Center);

    // Accept or reject tetahedron
    return (radius >= 0 && radius < mAlphaRadius);
  }

  //------------------------------------------------------------------------------------------------------------
  void TesselationUtilities3D::AddParticleVolume(ModelPart& rModelPart, std::vector<int>& addedParticle, double& particle_volume, const int id) {
    if (!addedParticle[id]) {
      ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin() + id;
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*it);
      addedParticle[id] = 1;
      particle_volume += particle.CalculateVolume();
    }
  }

} // namespace Kratos
