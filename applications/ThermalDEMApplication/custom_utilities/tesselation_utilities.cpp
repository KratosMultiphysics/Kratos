//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
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
#include "tesselation_utilities.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  TesselationUtilities::TesselationUtilities()
  {
    mUpdateVoronoi = false;
    mUpdatePorosiy = false;
    mAlphaRadius   = 0.0;
    mSwitches      = "";
  }

  TesselationUtilities::~TesselationUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::ExecuteInitialize(ModelPart& rModelPart, bool update_voronoi, bool update_porosity)
  {
    KRATOS_TRY

    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    // Set flags
    mUpdateVoronoi = update_voronoi;
    mUpdatePorosiy = update_porosity;

    // Initialize average porosity
    if (update_porosity)
      r_process_info[AVERAGE_POROSITY] = 1.0;

    // Set switches
    if (r_process_info[DOMAIN_SIZE] == 2) {
      // Switches for Triangle
      if      (update_voronoi)  mSwitches = "PQev";
      else if (update_porosity) mSwitches = "PQ";
    }
    else if (r_process_info[DOMAIN_SIZE] == 3) {
      // Switches for TetGen
      if      (update_voronoi)  mSwitches = "JQv";
      else if (update_porosity) mSwitches = "JQ";
    }
    else {
      KRATOS_ERROR << "Invalid domain size!";
    }

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::ExecuteInitializeSolutionStep(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    if (r_process_info[DOMAIN_SIZE] == 2)
      Triangulation(rModelPart);
    else if (r_process_info[DOMAIN_SIZE] == 3)
      Tetrahedralization(rModelPart);
    else
      KRATOS_ERROR << "Invalid domain size!";

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::Triangulation(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const int num_particles = rModelPart.NumberOfElements();
    if (num_particles < 3)
      return;

    // Create and clear IO
    struct triangulateio in, out, vorout;
    ClearTriangle(in);
    ClearTriangle(out);
    ClearTriangle(vorout);

    // Build input
    in.numberofpoints = num_particles;
    in.pointlist = (double*)malloc(sizeof(double) * in.numberofpoints * 2);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < in.numberofpoints; i++) {
      ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*it);

      particle.mDelaunayPointListIndex = i;

      const array_1d<double, 3>& coors = particle.GetParticleCoordinates();
      in.pointlist[2 * i + 0] = coors[0];
      in.pointlist[2 * i + 1] = coors[1];

      if (mUpdateVoronoi)
        particle.mNeighborVoronoiRadius.clear();
    }

    // Perform triangulation
    int fail = 0;
    try {
      triangulate(&mSwitches[0], &in, &out, &vorout);
    }
    catch (int error_code) {
      fail = error_code;
    }

    if (fail || out.numberoftriangles == 0 || in.numberofpoints != out.numberofpoints) {
      KRATOS_ERROR_IF(rModelPart.GetProcessInfo()[TIME_STEPS] == 1) << "Fail to generate triangulation!" << std::endl;
      KRATOS_WARNING("DEM") << std::endl;
      KRATOS_WARNING("DEM") << "Fail to generate triangulation! Results from previous successful triangulation will be used." << std::endl;
      KRATOS_WARNING("DEM") << std::endl;
      FreeTriangle(in);
      FreeTriangle(out);
      FreeTriangle(vorout);
      return;
    }

    // Update voronoi diagram
    if (mUpdateVoronoi)
      UpdateVoronoi2D(rModelPart, out, vorout);

    // Update average porosity
    if (mUpdatePorosiy)
      UpdatePorosity2D(rModelPart, out, vorout);

    // Free memory
    FreeTriangle(in);
    FreeTriangle(out);
    FreeTriangle(vorout);

    ClearTriangle(in);
    ClearTriangle(out);
    ClearTriangle(vorout);

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::Tetrahedralization(ModelPart& rModelPart)
  {
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

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < in.numberofpoints; i++) {
      ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*it);

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
      UpdateVoronoi3D(rModelPart, out);

    // Update average porosity
    if (mUpdatePorosiy)
      UpdatePorosity3D(rModelPart, out);

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Build a table for each particle with needed information from 2D voronoi diagram:
  * column 1: neighbor particle IDs
  * column 2: voronoi edge "radius" (length/2)
  */
  void TesselationUtilities::UpdateVoronoi2D(ModelPart& rModelPart, struct triangulateio& rOut, struct triangulateio& rVorOut) {
    ModelPart::ElementsContainerType& rElements = rModelPart.GetCommunicator().LocalMesh().Elements();

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < rOut.numberofpoints; i++) {
      // Particle corresponding to delaunay vertex
      ModelPart::ElementsContainerType::iterator it = rElements.ptr_begin() + i;
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*it);

      for (unsigned int j = 0; j < rOut.numberofedges; j++) {
        // Vertices of delaunay edge
        const int vd1 = rOut.edgelist[2 * j + 0] - 1;
        const int vd2 = rOut.edgelist[2 * j + 1] - 1;

        // Check if delaunay edge contains current point:
        // Only look at one vertex to avoid repeating the process for both vetices of the same edge
        if (vd1 == particle.mDelaunayPointListIndex) {
          // Particle corresponding to neighboring delaunay vertex
          ModelPart::ElementsContainerType::iterator itn = rElements.ptr_begin() + vd2;
          ThermalSphericParticle& neighbor = dynamic_cast<ThermalSphericParticle&> (*itn);

          // Vertices of voronoi edge dual to delaunay edge
          const int vv1 = rVorOut.edgelist[2 * j + 0] - 1;
          const int vv2 = rVorOut.edgelist[2 * j + 1] - 1;

          // Check for bounded or unbounded voronoi cell to compute radius of voronoi edge
          // (unbounded voronoi edge has a negative vertice ID)
          // Assumption: Unbounded edge is only considered when intersected by delaunay edge.
          //             In this case, the radius is the perpendicular distance betweem the
          //             bounded voronoi vertex and the delaunay edge.
          double radius = 0.0;

          if (vv1 >= 0 && vv2 >= 0) { // bounded
            // Coordinates of bounded voronoi edge vertices
            const double xv1 = rVorOut.pointlist[2 * vv1 + 0];
            const double yv1 = rVorOut.pointlist[2 * vv1 + 1];
            const double xv2 = rVorOut.pointlist[2 * vv2 + 0];
            const double yv2 = rVorOut.pointlist[2 * vv2 + 1];
            radius = sqrt(pow(xv2 - xv1, 2) + pow(yv2 - yv1, 2)) / 2.0;
          }
          else { // unbounded
            // Coordiantes of any delaunay edge vertex
            const double xd = rOut.pointlist[2 * vd1 + 0];
            const double yd = rOut.pointlist[2 * vd1 + 1];

            // Coordinates of bounded voronoi vertex
            int vv = vv1;
            if (vv1 < 0) vv = vv2;
            const double xv = rVorOut.pointlist[2 * vv + 0];
            const double yv = rVorOut.pointlist[2 * vv + 1];

            // Direction of unbounded voronoi edge
            const double dirx = rVorOut.normlist[2 * j + 0];
            const double diry = rVorOut.normlist[2 * j + 1];

            // Check if delaunay edge intersects the unbounded voronoi edge
            const double dotProd = dirx * (xd - xv) + diry * (yd - yv);

            if (dotProd > 0.0)
              radius = dotProd / sqrt(dirx * dirx + diry * diry);
          }

          // Add info to table of both particles (current and neighbor to avoid repeating the process for both vertices of the same edge)
          particle.mNeighborVoronoiRadius[vd2] = radius;
          neighbor.mNeighborVoronoiRadius[vd1] = radius;
        }
      }
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Build a table for each particle with needed information from 3D voronoi diagram:
  * column 1: neighbor particle IDs
  * column 2: voronoi face equivalent radius (for a circle with the same area)
  */
  void TesselationUtilities::UpdateVoronoi3D(ModelPart& rModelPart, struct tetgenio& rOut) {
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
        for (unsigned int j = 2; j <= num_face_edges; j++) {
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

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Compute total area of delaunay triangles and particles (circular cross-section) to obtain global avarege porosity.
  * Oversized triangles may be removed using the alpha-shape method.
  */
  void TesselationUtilities::UpdatePorosity2D(ModelPart& rModelPart, struct triangulateio& rOut, struct triangulateio& rVorOut) {
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    std::vector<int> addedParticle(rOut.numberofpoints,0);
    double total_area    = 0.0;
    double particle_area = 0.0;

    // Compute alpha radius
    // Assumption: Only in the 1st step - keeping alpha radius constant throghout simulation
    if (r_process_info[POSORITY_METHOD].compare("average_alpha_shape") == 0 && r_process_info[TIME_STEPS] == 1)
      ComputeAlphaRadius2D(rModelPart, rOut);

    // Accumulate total area of triangles and particles
    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < rOut.numberoftriangles; i++) {
      // Get vertices IDs
      const int v1 = rOut.trianglelist[3 * i + 0] - 1;
      const int v2 = rOut.trianglelist[3 * i + 1] - 1;
      const int v3 = rOut.trianglelist[3 * i + 2] - 1;

      // Get vertices coordinates
      const double x1 = rOut.pointlist[2 * v1 + 0];
      const double y1 = rOut.pointlist[2 * v1 + 1];
      const double x2 = rOut.pointlist[2 * v2 + 0];
      const double y2 = rOut.pointlist[2 * v2 + 1];
      const double x3 = rOut.pointlist[2 * v3 + 0];
      const double y3 = rOut.pointlist[2 * v3 + 1];

      // Perform alpha-shape
      if (r_process_info[POSORITY_METHOD].compare("average_alpha_shape") == 0) {
        std::vector<double> coords = { x1, y1, x2, y2, x3, y3 };
        if (!AlphaShape2D(coords)) {
          continue;
        }
      }

      // Add triangle area
      total_area += fabs(0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)));

      // Add particles area
      AddParticleArea(rModelPart, addedParticle, particle_area, v1);
      AddParticleArea(rModelPart, addedParticle, particle_area, v2);
      AddParticleArea(rModelPart, addedParticle, particle_area, v3);
    }

    // Set new average porosity
    r_process_info[AVERAGE_POROSITY] = 1.0 - particle_area / total_area;
  }

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Compute total volume of delaunay tetahedra and particles to obtain global avarege porosity.
  * Oversized tetahedra may be removed using the alpha-shape method.
  */
  void TesselationUtilities::UpdatePorosity3D(ModelPart& rModelPart, struct tetgenio& rOut) {
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    std::vector<int> addedParticle(rOut.numberofpoints,0);
    double total_volume    = 0.0;
    double particle_volume = 0.0;

    // Compute alpha radius
    // Assumption: Only in the 1st step - keeping alpha radius constant throghout simulation
    if (r_process_info[POSORITY_METHOD].compare("average_alpha_shape") == 0 && r_process_info[TIME_STEPS] == 1)
      ComputeAlphaRadius3D(rModelPart, rOut);

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
      if (r_process_info[POSORITY_METHOD].compare("average_alpha_shape") == 0) {
        std::vector<double> coords = { x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 };
        if (!AlphaShape3D(coords)) {
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
      total_volume += fabs(GeometryFunctions::DotProduct(a, bxc)) / 6.0;

      // Add particles volume
      AddParticleArea(rModelPart, addedParticle, particle_volume, v1);
      AddParticleArea(rModelPart, addedParticle, particle_volume, v2);
      AddParticleArea(rModelPart, addedParticle, particle_volume, v3);
      AddParticleArea(rModelPart, addedParticle, particle_volume, v4);
    }

    // Set new average porosity
    r_process_info[AVERAGE_POROSITY] = 1.0 - particle_volume / total_volume;
  }

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Compute mean mesh size for alpha-shape.
  * The mean mesh size is taken as the average of the smallest side of each triangle.
  */
  void TesselationUtilities::ComputeAlphaRadius2D(ModelPart& rModelPart, struct triangulateio& rOut) {
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    double MeanMeshSize = 0.0;

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < rOut.numberoftriangles; i++) {
      // Get vertices IDs
      const int v1 = rOut.trianglelist[3 * i + 0] - 1;
      const int v2 = rOut.trianglelist[3 * i + 1] - 1;
      const int v3 = rOut.trianglelist[3 * i + 2] - 1;

      // Get vertices coordinates
      const double x1 = rOut.pointlist[2 * v1 + 0];
      const double y1 = rOut.pointlist[2 * v1 + 1];
      const double x2 = rOut.pointlist[2 * v2 + 0];
      const double y2 = rOut.pointlist[2 * v2 + 1];
      const double x3 = rOut.pointlist[2 * v3 + 0];
      const double y3 = rOut.pointlist[2 * v3 + 1];

      // Get minimum edge length
      std::vector<double>len;
      len.push_back(sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)));
      len.push_back(sqrt(pow(x3 - x1, 2) + pow(y3 - y1, 2)));
      len.push_back(sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2)));

      MeanMeshSize += *std::min_element(len.begin(), len.end());
    }

    // Average minimum length
    MeanMeshSize /= rOut.numberoftriangles;

    // Alpha radius
    mAlphaRadius = MeanMeshSize * r_process_info[ALPHA_SHAPE_PARAMETER];
  }

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Compute mean mesh size for alpha-shape.
  * The mean mesh size is taken as the average of the smallest edge of each tetahedron.
  */
  void TesselationUtilities::ComputeAlphaRadius3D(ModelPart& rModelPart, struct tetgenio& rOut) {
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    double MeanMeshSize = 0.0;

    for (unsigned int i = 0; i < rOut.numberoftetrahedra; i++) {
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

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Perform alpha-shape verification on a delaunay triangle to remove distorted shapes.
  */
  bool TesselationUtilities::AlphaShape2D(std::vector<double>& coords) {
    const double x1 = coords[0];
    const double y1 = coords[1];
    const double x2 = coords[2];
    const double y2 = coords[3];
    const double x3 = coords[4];
    const double y3 = coords[5];

    // Calculate Jacobian
    BoundedMatrix<double, 2, 2> J;
    J(0, 0) = x2 - x1;
    J(0, 1) = y2 - y1;
    J(1, 0) = x3 - x1;
    J(1, 1) = y3 - y1;
    J *= 2.0;

    // Calculate the determinant (volume/2)
    const double vol = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);

    // Calculate the inverse of the Jacobian
    BoundedMatrix<double, 2, 2> Jinv;
    Jinv(0, 0) =  J(1, 1);
    Jinv(0, 1) = -J(0, 1);
    Jinv(1, 0) = -J(1, 0);
    Jinv(1, 1) =  J(0, 0);
    Jinv /= vol;

    // Calculate circle center
    Vector Center = ZeroVector(2);
    Center[0] += (x2 * x2);
    Center[0] -= (x1 * x1);
    Center[0] += (y2 * y2);
    Center[0] -= (y1 * y1);
    Center[1] += (x3 * x3);
    Center[1] -= (x1 * x1);
    Center[1] += (y3 * y3);
    Center[1] -= (y1 * y1);
    Center = prod(Jinv, Center);

    // Calculate circle radius
    Center[0] -= x1;
    Center[1] -= y1;
    const double radius = norm_2(Center);

    // Accept or reject triangle
    return (radius >= 0 && radius < mAlphaRadius);
  }

  //-----------------------------------------------------------------------------------------------------------------------
  /*
  * Perform alpha-shape verification on a delaunay tetahedron to remove distorted shapes.
  */
  bool TesselationUtilities::AlphaShape3D(std::vector<double>& coords) {
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

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::AddParticleArea(ModelPart& rModelPart, std::vector<int>& addedParticle, double& particle_area, const int id) {
    if (!addedParticle[id]) {
      ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin() + id;
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*it);
      addedParticle[id] = 1;
      const double r = particle.GetRadius();
      particle_area += Globals::Pi * r * r;
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::AddParticleVolume(ModelPart& rModelPart, std::vector<int>& addedParticle, double& particle_volume, const int id) {
    if (!addedParticle[id]) {
      ModelPart::ElementsContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin() + id;
      ThermalSphericParticle& particle = dynamic_cast<ThermalSphericParticle&> (*it);
      addedParticle[id] = 1;
      particle_volume += particle.CalculateVolume();
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::ClearTriangle(struct triangulateio& rTr)
  {
    KRATOS_TRY

    rTr.pointlist                  = (REAL*) NULL;
    rTr.pointattributelist         = (REAL*) NULL;
    rTr.pointmarkerlist            = (int*)  NULL;
    rTr.numberofpoints             = 0;
    rTr.numberofpointattributes    = 0;

    rTr.trianglelist               = (int*)  NULL;
    rTr.triangleattributelist      = (REAL*) NULL;
    rTr.trianglearealist           = (REAL*) NULL;
    rTr.neighborlist               = (int*)  NULL;
    rTr.numberoftriangles          = 0;
    rTr.numberofcorners            = 3; //for three node triangles
    rTr.numberoftriangleattributes = 0;

    rTr.segmentlist                = (int*) NULL;
    rTr.segmentmarkerlist          = (int*) NULL;
    rTr.numberofsegments           = 0;

    rTr.holelist                   = (REAL*) NULL;
    rTr.numberofholes              = 0;

    rTr.regionlist                 = (REAL*) NULL;
    rTr.numberofregions            = 0;

    rTr.edgelist                   = (int*)  NULL;
    rTr.edgemarkerlist             = (int*)  NULL;
    rTr.normlist                   = (REAL*) NULL;
    rTr.numberofedges              = 0;

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void TesselationUtilities::FreeTriangle(struct triangulateio& rTr)
  {
    KRATOS_TRY

    if (rTr.numberoftriangles) {
      if (rTr.trianglelist)          trifree(rTr.trianglelist);
      if (rTr.triangleattributelist) trifree(rTr.triangleattributelist);
      if (rTr.trianglearealist)      trifree(rTr.trianglearealist);
      if (rTr.neighborlist)          trifree(rTr.neighborlist);
    }
    if (rTr.segmentlist)       trifree(rTr.segmentlist);
    if (rTr.segmentmarkerlist) trifree(rTr.segmentmarkerlist);
    if (rTr.holelist) {
      delete[] rTr.holelist;
      rTr.numberofholes = 0;
    }
    if (rTr.regionlist) {
      delete[] rTr.regionlist;
      rTr.numberofregions = 0;
    }
    if (rTr.edgelist)       trifree(rTr.edgelist);
    if (rTr.edgemarkerlist) trifree(rTr.edgemarkerlist);
    if (rTr.normlist)       trifree(rTr.normlist);
    if (rTr.numberofpoints) {
      if (rTr.pointlist)          trifree(rTr.pointlist);
      if (rTr.pointattributelist) trifree(rTr.pointattributelist);
      if (rTr.pointmarkerlist)    trifree(rTr.pointmarkerlist);
    }

    KRATOS_CATCH("")
  }

} // namespace Kratos