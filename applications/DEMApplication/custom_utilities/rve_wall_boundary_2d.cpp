//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//
//  RVE boundary nomenclature:
//
//               Y+ (Ymax)
//   Vertex 4 --------------- Vertex 3
//           |               |
//           |               |
// X- (Xmin) |      RVE      | X+ (Xmax)
//           |               |
//           |               | 
//   Vertex 1 --------------- Vertex 2
//               Y- (Ymin)
// 
#include "rve_wall_boundary_2d.h"

namespace Kratos
{
    //------------------------------------------------------------------------------------------------------------
    // Create vectors of RVE wall elements.
    // The RVE has four walls (called XMin, XMax, YMin, YMax): the four straight lines forming the boundary.
    // Each of these walls can be made of multiple FEM elements - this function arrange them into four vectors.
    // TODO: THIS FUNCTION IS SPECIFIC FOR SQUARE OR SHEARED RVEs. IT SHOULD BE ADAPTED TO CONSIDER ANY SHAPE.
    void RVEWallBoundary2D::AssembleWallElementVectors(void) {
        ModelPart::ConditionsContainerType &r_conditions = mFemModelPart->GetCommunicator().LocalMesh().Conditions();

        // Determine max wall slope
        double slope_max = -DBL_MAX;

        for (unsigned int i = 0; i < mNumWallElems; i++) {
            ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
            DEMWall *p_wall = dynamic_cast<DEMWall *>(&(*it));
            Condition::GeometryType &geom = p_wall->GetGeometry();

            // Element end coordinates
            const double x1 = geom[0][0]; const double y1 = geom[0][1];
            const double x2 = geom[1][0]; const double y2 = geom[1][1];

            // Calculate element slope
            double slope;
            if (x1 == x2) slope = DBL_MAX; // Vertical element
            else slope = (y2-y1) / (x2-x1);

            // Update max wall slope
            if (std::abs(slope) > slope_max) slope_max = std::abs(slope);
        }

        // Determine X and Y (ie, "vertical" and "horizontal") elements,
        // and determine their minimum intersection projections with global X and Y axes
        std::vector<DEMWall*> wall_elems_x;
        std::vector<DEMWall*> wall_elems_y;
        double intersect_xmin = DBL_MAX;
        double intersect_ymin = DBL_MAX;

        for (unsigned int i = 0; i < mNumWallElems; i++) {
            ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
            DEMWall *p_wall = dynamic_cast<DEMWall *>(&(*it));
            Condition::GeometryType &geom = p_wall->GetGeometry();

            // Element end coordinates
            const double x1 = geom[0][0]; const double y1 = geom[0][1];
            const double x2 = geom[1][0]; const double y2 = geom[1][1];

            // Calculate element slope
            double slope;
            if (x1 == x2) slope = DBL_MAX; // Vertical element
            else slope = (y2-y1) / (x2-x1);

            // Compare with max wall slope
            if (EqualValues(std::abs(slope), slope_max)) { // Element belongs to an X wall (eg, vertical elements)
                wall_elems_x.push_back(p_wall);

                // Projection of element intersection with X axis
                const double intersect_x = x1 - y1 / slope;

                // Update minimum intersection
                if (intersect_x < intersect_xmin) intersect_xmin = intersect_x;
            }
            else { // Element belongs to an Y wall (eg, horizontal elements)
                wall_elems_y.push_back(p_wall);

                // Projection of element intersection with Y axis
                const double intersect_y = y1 - slope * x1;

                // Update minimum intersection
                if (intersect_y < intersect_ymin) intersect_ymin = intersect_y;
            }
        }

        // Determine on which side (negative or positive) of X walls the elements are located
        for (unsigned int i = 0; i < wall_elems_x.size(); i++) {
            Condition::GeometryType &geom = wall_elems_x[i]->GetGeometry();

            // Element end coordinates
            const double x1 = geom[0][0]; const double y1 = geom[0][1];
            const double x2 = geom[1][0]; const double y2 = geom[1][1];

            // Calculate element intersection projection with global X axis
            double slope;
            if (x1 == x2) slope = DBL_MAX;
            else slope = (y2-y1) / (x2-x1);
            const double intersect_x = x1 - y1 / slope;

            // Compare with minimum intersection and insert element into corresponding wall vector
            if (EqualValues(intersect_x, intersect_xmin))
                mWallXMin.push_back(wall_elems_x[i]);
            else
                mWallXMax.push_back(wall_elems_x[i]);
        }

        // Determine on which side (negative or positive) of Y walls the elements are located
        for (unsigned int i = 0; i < wall_elems_y.size(); i++) {
            Condition::GeometryType &geom = wall_elems_y[i]->GetGeometry();

            // Element end coordinates
            const double x1 = geom[0][0]; const double y1 = geom[0][1];
            const double x2 = geom[1][0]; const double y2 = geom[1][1];

            // Calculate element intersection projection with global Y axis
            double slope;
            if (x1 == x2) slope = DBL_MAX;
            else slope = (y2-y1) / (x2-x1);
            const double intersect_y = y1 - slope * x1;

            // Compare with minimum intersection and insert element into corresponding wall vector
            if (EqualValues(intersect_y, intersect_ymin))
                mWallYMin.push_back(wall_elems_y[i]);
            else
                mWallYMax.push_back(wall_elems_y[i]);
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Identify and store the updated coordinates of RVE vertices in counterclockwise order.
    // TODO: THIS FUNCTION ONLY WORKS FOR SQUARED AND SHEARED RVEs WITH INCLINED LATERAL WALLS. IT SHOULD BE GENERALIZED FOR ANY RVE SHAPE.
    void RVEWallBoundary2D::SetVertexCoordinates(void) {
        // Vertex coordinates
        double x1, x2, x3, x4, y1, y2, y3, y4;

        // Coordinates of end nodes (x1,y1 and x2,y2) of any element (in this case, with index [0]) of RVE walls (xmin, xmax, ymin, ymax)
        const double xmin_x1 = mWallXMin[0]->GetGeometry()[0][0]; const double xmin_x2 = mWallXMin[0]->GetGeometry()[1][0];
        const double xmin_y1 = mWallXMin[0]->GetGeometry()[0][1]; const double xmin_y2 = mWallXMin[0]->GetGeometry()[1][1];
        const double xmax_x1 = mWallXMax[0]->GetGeometry()[0][0]; const double xmax_x2 = mWallXMax[0]->GetGeometry()[1][0];
        const double ymin_y1 = mWallYMin[0]->GetGeometry()[0][1]; const double ymin_y2 = mWallYMin[0]->GetGeometry()[1][1];
        const double ymax_y1 = mWallYMax[0]->GetGeometry()[0][1]; const double ymax_y2 = mWallYMax[0]->GetGeometry()[1][1];

        // Check RVE shape
        bool is_square = (EqualValues(xmin_x1,xmin_x2) && EqualValues(xmax_x1,xmax_x2) && EqualValues(ymin_y1,ymin_y2) && EqualValues(ymax_y1,ymax_y2));

        if (is_square) {
            // Set vertex coordinates
            x1 = xmin_x1; y1 = ymin_y1;
            x2 = xmax_x1; y2 = ymin_y1;
            x3 = xmax_x1; y3 = ymax_y1;
            x4 = xmin_x1; y4 = ymax_y1;
        }
        else { // Sheared RVE: Assuming an RVE with inclined lateral walls (Xmin, Xmax), and horizontal Ymin and Ymax walls
            // Determine min/max X coordinates of lower horizontal (Ymin) wall
            double ymin_x_min =  DBL_MAX;
            double ymin_x_max = -DBL_MAX;
            for (unsigned int i = 0; i < mWallYMin.size(); i++) {
                Condition::GeometryType& geom = mWallYMin[i]->GetGeometry();
                const double ymin_x1 = geom[0][0];
                const double ymin_x2 = geom[1][0];
                if (ymin_x1 < ymin_x_min) ymin_x_min = ymin_x1;
                if (ymin_x1 > ymin_x_max) ymin_x_max = ymin_x1;
                if (ymin_x2 < ymin_x_min) ymin_x_min = ymin_x2;
                if (ymin_x2 > ymin_x_max) ymin_x_max = ymin_x2;
            }

            // Determine min/max X coordinates of upper horizontal (Ymax) wall
            const double slope = (xmin_x2-xmin_x1) / (xmin_y2-xmin_y1); // Inclination of lateral walls
            const double height = ymax_y1 - ymin_y1; // Total height of RVE
            const double ymax_x_min = ymin_x_min + slope * height; // Min X coordinates of upper horizontal (Ymax) wall
            const double ymax_x_max = ymin_x_max + slope * height; // Max X coordinates of upper horizontal (Ymax) wall

            // Set vertex coordinates
            x1 = ymin_x_min; y1 = ymin_y1;
            x2 = ymin_x_max; y2 = ymin_y1;
            x3 = ymax_x_max; y3 = ymax_y1;
            x4 = ymax_x_min; y4 = ymax_y1;
        }

        // Ensure that vertices are in counterclockwise order
        SortVerticesCounterClockwise(x1, x2, x3, x4, y1, y2, y3, y4);

        // Set vertex coordinates
        mVertexCoords(0,0) = x1; mVertexCoords(1,0) = y1;
        mVertexCoords(0,1) = x2; mVertexCoords(1,1) = y2;
        mVertexCoords(0,2) = x3; mVertexCoords(1,2) = y3;
        mVertexCoords(0,3) = x4; mVertexCoords(1,3) = y4;
    }

    //------------------------------------------------------------------------------------------------------------
    // TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    void RVEWallBoundary2D::SetVertexCoordinatesInner(void) {
    }

    //------------------------------------------------------------------------------------------------------------
    // Process and sum up global RVE results accumulated from individual particle and interaction contributions.
    void RVEWallBoundary2D::ProcessGlobalResults(void) {
        for (unsigned int i = 0; i < mNumParticles; i++) {
            ModelPart::ElementsContainerType::iterator it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
            SphericParticle &particle = dynamic_cast<SphericParticle&>(*it);

            // Particle properties
            const int id1 = particle.GetId();
            const double r1 = particle.GetRadius();
            const double x1 = particle.GetGeometry()[0][0];
            const double y1 = particle.GetGeometry()[0][1];
            std::vector<double> coords1 = {x1, y1};
            bool is_inner_particle = true;

            // Accumulate particle properties
            mAvgRadius     += r1;
            mVolSolid      += ComputeVolumeParticle(particle);
            mVolSolidInner += ComputeVolumeParticleInner(particle);

            // Loop over contacts with walls
            for (unsigned int j = 0; j < particle.mNeighbourRigidFaces.size(); j++) {
                // Neighbor properties
                const int id2 = particle.mNeighbourRigidFaces[j]->GetId();
                const double x2 = particle.mNeighbourRigidFaces[j]->GetGeometry()[0][0];
                const double y2 = particle.mNeighbourRigidFaces[j]->GetGeometry()[0][1];
                const double x3 = particle.mNeighbourRigidFaces[j]->GetGeometry()[1][0];
                const double y3 = particle.mNeighbourRigidFaces[j]->GetGeometry()[1][1];

                // Check for existing contact
                if (particle.mBallToRigidFaceStoredInfo.find(id2) == particle.mBallToRigidFaceStoredInfo.end()) continue;
                const double indent = particle.mBallToRigidFaceStoredInfo[id2].indentation;
                if (indent <= 0.0) continue;

                // Increment number of contacts
                mAvgCoordNum++;
                mNumContacts++;
                is_inner_particle = false;

                // Normal vector
                const double d  = r1 - indent;
                const double nx = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][1];
                std::vector<double> normal = {nx, ny};
                std::vector<double> branch = {d * nx, d * ny};

                // Update rose diagram
                AddContactToRoseDiagram(mRoseDiagram, normal);

                // Applied wall force (normal component)
                const double fx = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[0];
                const double fy = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[1];
                std::vector<double> global_contact_force = {fx, fy};
                mWallForces += std::abs(fx * nx + fy * ny);

                // Contact results
                std::vector<double> chain{x1, y1, 0.0, x1+branch[0], y1+branch[1], 0.0, fx, fy, 0.0};
                mContactChain.insert(mContactChain.end(), chain.begin(), chain.end());

                // Tensors
                for (int k = 0; k < mDim; k++) {
                    for (int l = 0; l < mDim; l++) {
                        mFabricTensor(k,l) += normal[k] * normal[l];
                        mStressTensor(k,l) += branch[k] * global_contact_force[l];
                    }
                }
            }
            if (is_inner_particle) {
                mNumParticlesInner++;
            }

            // Loop over contacts with particles
            for (unsigned int j = 0; j < particle.mNeighbourElements.size(); j++) {
                // Neighbor properties
                const int id2 = particle.mNeighbourElements[j]->GetId();
                const double r2 = particle.mNeighbourElements[j]->GetRadius();
                const double x2 = particle.mNeighbourElements[j]->GetGeometry()[0][0];
                const double y2 = particle.mNeighbourElements[j]->GetGeometry()[0][1];
                std::vector<double> coords2 = {x2, y2};

                // Check for existing contact
                if (particle.mBallToBallStoredInfo.find(id2) == particle.mBallToBallStoredInfo.end()) continue;
                const double indent = particle.mBallToBallStoredInfo[id2].indentation;
                if (indent <= 0.0) continue;

                // Increment number of contacts
                mAvgCoordNum++;
                if (is_inner_particle) mAvgCoordNumInner++;

                // Normal vector
                const double d  = r1 + r2 - indent;
                const double nx = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][1];
                std::vector<double> normal = {nx, ny};
                std::vector<double> branch = {d * nx, d * ny};

                // Update rose diagram
                AddContactToRoseDiagram(mRoseDiagram, normal);
                if (is_inner_particle) AddContactToRoseDiagram(mRoseDiagramInner, normal);

                // Unique contacts (each binary contact evaluated only once)
                if (id1 < id2) {
                    // Check for inner contact
                    const double inner_contact_len = ComputeBranchLengthInner(coords1, coords2);
                    const double inner_contact_ratio = inner_contact_len / d;
                    bool is_inner_contact = (inner_contact_len != 0.0);

                    // Increment number of unique contacts
                    mNumContacts++;
                    if (is_inner_contact) mNumContactsInner++;

                    // Contact results
                    const double fx = particle.mBallToBallStoredInfo[id2].global_contact_force[0];
                    const double fy = particle.mBallToBallStoredInfo[id2].global_contact_force[1];
                    std::vector<double> global_contact_force = {fx, fy};
                    std::vector<double> chain{x1, y1, 0.0, x2, y2, 0.0, fx, fy, 0.0};
                    mContactChain.insert(mContactChain.end(), chain.begin(), chain.end());

                    // Tensors
                    for (unsigned int k = 0; k < mDim; k++) {
                        for (unsigned int l = 0; l < mDim; l++) {
                            mFabricTensor(k,l) += normal[k] * normal[l];
                            mStressTensor(k,l) += branch[k] * global_contact_force[l];
                            if (is_inner_contact) {
                                mFabricTensorInner(k,l) += normal[k] * normal[l];
                                mStressTensorInner(k,l) += branch[k] * global_contact_force[l] * inner_contact_ratio;
                            }
                        }
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double RVEWallBoundary2D::ComputeBranchLengthInner(std::vector<double>& coords1, std::vector<double>& coords2) {
        return 0.0;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute particle volume, currently considering circle area in 2D.
    double RVEWallBoundary2D::ComputeVolumeParticle(SphericParticle& particle) {
        return Globals::Pi * particle.GetRadius() * particle.GetRadius();
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute particle volume inside RVE inner volume (areas in 2D).
    // TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double RVEWallBoundary2D::ComputeVolumeParticleInner(SphericParticle& particle) {
        return 0.0;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the total area of RVE.
    double RVEWallBoundary2D::ComputeVolumeRVE(void) {
        const double x1 = mVertexCoords(0,0); const double y1 = mVertexCoords(1,0);
        const double x2 = mVertexCoords(0,1); const double y2 = mVertexCoords(1,1);
        const double x3 = mVertexCoords(0,2); const double y3 = mVertexCoords(1,2);
        const double x4 = mVertexCoords(0,3); const double y4 = mVertexCoords(1,3);
        return ComputeAreaFromVertices(x1, x2, x3, x4, y1, y2, y3, y4);
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the inner area of RVE, considered as the effective area without wall boundary effects.
    // METHOD 1: Applies an offset proportional to the particle radius to the RVE walls to obtain an inner quadrilateral.
    //           The optimal offset value was determined as 0.8R.
    //           This simplistic method was used in the simulations of the reference works.
    // METHOD 2: Perform a Delaunay triangulation over the center of the inner particles to obtain the convex hull.
    //           (This has been removed!)
    double RVEWallBoundary2D::ComputeVolumeRVEInner(void) {
        // Offset with respect to walls
        const double offset = mInnerVolOffset * mAvgRadius;

        // Vertex coordinates
        double x1 = mVertexCoords(0,0); double y1 = mVertexCoords(1,0);
        double x2 = mVertexCoords(0,1); double y2 = mVertexCoords(1,1);
        double x3 = mVertexCoords(0,2); double y3 = mVertexCoords(1,2);
        double x4 = mVertexCoords(0,3); double y4 = mVertexCoords(1,3);

        // Wall direction vectors
        double w1x = x2-x1; double w1y = y2-y1;
        double w2x = x3-x2; double w2y = y3-y2;
        double w3x = x4-x3; double w3y = y4-y3;
        double w4x = x1-x4; double w4y = y1-y4;
        
        // Wall lengths
        double len1 = std::sqrt(w1x*w1x + w1y*w1y);
        double len2 = std::sqrt(w2x*w2x + w2y*w2y);
        double len3 = std::sqrt(w3x*w3x + w3y*w3y);
        double len4 = std::sqrt(w4x*w4x + w4y*w4y);

        // Inward vectors perpendicular to the walls (normalized and scaled by offset)
        double n1x = -(w1y/len1)*offset; double n1y = (w1x/len1)*offset;
        double n2x = -(w2y/len2)*offset; double n2y = (w2x/len2)*offset;
        double n3x = -(w3y/len3)*offset; double n3y = (w3x/len3)*offset;
        double n4x = -(w4y/len4)*offset; double n4y = (w4x/len4)*offset;

        // Vertex coordinates of inner quadrilateral
        double x1i = x1 + n1x + n4x; double y1i = y1 + n1y + n4y;
        double x2i = x2 + n1x + n2x; double y2i = y2 + n1y + n2y;
        double x3i = x3 + n2x + n3x; double y3i = y3 + n2y + n3y;
        double x4i = x4 + n3x + n4x; double y4i = y4 + n3y + n4y;

        // Compute area of inner quadrilateral
        return ComputeAreaFromVertices(x1, x2, x3, x4, y1, y2, y3, y4);
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute RVE porosity considering its inner volume.
    double RVEWallBoundary2D::ComputePorosityInner(void) {
        return RVEUtilities::ComputePorosityInner();

        // OLD ALTERNATIVE METHOD: ONLY FOR SQUARE RVEs (NEEDS TO BE CHECKED)
        /*
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        // Vertex coordinates
        double x1 = mVertexCoords(0,0); double y1 = mVertexCoords(1,0);
        double x2 = mVertexCoords(0,1); double y2 = mVertexCoords(1,1);
        double x3 = mVertexCoords(0,2); double y3 = mVertexCoords(1,2);
        double x4 = mVertexCoords(0,3); double y4 = mVertexCoords(1,3);

        // Check RVE shape
        bool is_square = (EqualValues(y1,y2) && EqualValues(y3,y4) && EqualValues(x1,x4) && EqualValues(x2,x3));
        if (!is_square)
            return 0.0;
    
        // Inner boundaries and corners
        const double offset = mInnerVolOffset * mRVE_MeanRadius;

        const double xmin = mRVE_WallXMin[0]->GetGeometry()[0][0] + offset;
        const double xmax = mRVE_WallXMax[0]->GetGeometry()[0][0] - offset;
        const double ymin = mRVE_WallYMin[0]->GetGeometry()[0][1] + offset;
        const double ymax = mRVE_WallYMax[0]->GetGeometry()[0][1] - offset;

        const double corner_1[3] = { xmin, ymin, 0.0 };
        const double corner_2[3] = { xmax, ymin, 0.0 };
        const double corner_3[3] = { xmin, ymax, 0.0 };
        const double corner_4[3] = { xmax, ymax, 0.0 };

        // Inner volume
        const double inner_vol = (xmax - xmin) * (ymax - ymin);

        // Compute solid volume
        mRVE_VolSolidInner = 0.0;

      // Loop over all particles
      for (int i = 0; i < mListOfSphericParticles.size(); i++) {
        const double r = mListOfSphericParticles[i]->GetRadius();
        auto& central_node = mListOfSphericParticles[i]->GetGeometry()[0];
        const double x = central_node[0];
        const double y = central_node[1];
        const double z = central_node[2];
        const double coords[3] = { x, y, z };

        // Distances from particle center to corners
        const double dist_corner_1 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_1);
        const double dist_corner_2 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_2);
        const double dist_corner_3 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_3);
        const double dist_corner_4 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_4);

        // Area inside inner region
        double area_inside = 0.0;
        
        // 1 - Particle inside region
        if (x + r <= xmax && x - r >= xmin && y + r <= ymax && y - r >= ymin) {
          area_inside = Globals::Pi * r * r;
        }

        // 2A - Cell corner (bottom-left) inside circle
        else if (dist_corner_1 < r) {
          const double dx = x - xmin;
          const double dy = y - ymin;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx + dx;
          const double ly = sy + dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // 2B - Cell corner (bottom-right) inside circle
        else if (dist_corner_2 < r) {
          const double dx = x - xmax;
          const double dy = y - ymin;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx - dx;
          const double ly = sy + dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // 2C - Cell corner (top-left) inside circle
        else if (dist_corner_3 < r) {
          const double dx = x - xmin;
          const double dy = y - ymax;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx + dx;
          const double ly = sy - dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // 2D - Cell corner (top-right) inside circle
        else if (dist_corner_4 < r) {
          const double dx = x - xmax;
          const double dy = y - ymax;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx - dx;
          const double ly = sy - dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // Particle crossed by edges
        else {
          bool is_inside = false;

          // Left edge
          if (std::abs(x - xmin) < r && y > ymin && y < ymax) {
            const double ident = r - std::abs(x - xmin);
            const double ovlp  = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (x > xmin) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Right edge
          if (std::abs(x - xmax) < r && y > ymin && y < ymax) {
            const double ident = r - std::abs(x - xmax);
            const double ovlp  = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (x < xmax) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Lower edge
          if (std::abs(y - ymin) < r && x > xmin && x < xmax) {
            const double ident = r - std::abs(y - ymin);
            const double ovlp = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (y > ymin) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Top edge
          if (std::abs(y - ymax) < r && x > xmin && x < xmax) {
            const double ident = r - std::abs(y - ymax);
            const double ovlp = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (y < ymax) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Add particle area if it is inside cell
          if (is_inside)
            area_inside = area_inside + Globals::Pi * r * r;
        }

        // Remove overlaps: loop over possible particle neighbors
        for (int j = i+1; j < mListOfSphericParticles.size(); j++) {
          const double rn = mListOfSphericParticles[j]->GetRadius();
          const double xn = mListOfSphericParticles[j]->GetGeometry()[0][0];
          const double yn = mListOfSphericParticles[j]->GetGeometry()[0][1];
          const double zn = mListOfSphericParticles[j]->GetGeometry()[0][2];

          // Identation
          auto& neighbour_node = mListOfSphericParticles[j]->GetGeometry()[0];
          array_1d<double, 3> v;
          noalias(v) = central_node.Coordinates() - neighbour_node.Coordinates();
          const double d = DEM_MODULUS_3(v);
          const double ident = r + rn - d;

          if (ident <= 0.0)
            continue;

          // Overlap area
          double area_ovlp = r*r * std::acos((d*d + r*r - rn*rn) / (2.0 * d * r)) + rn*rn * std::acos((d*d + rn*rn - r*r) / (2.0 * d * rn)) - 0.5 * sqrt((d + r + rn) * (-d + r + rn) * (d + r - rn) * (d - r + rn));

          // Unit vector
          array_1d<double, 3> vu = v;
          GeometryFunctions::normalize(vu);

          // Contact radius
          const double rc = sqrt(r*r - pow((r*r - rn* rn + d*d)/(2*d),2.0));

          // Coordinates of overlap center
          const double d_1O = sqrt(r*r - rc*rc);
          array_1d<double,3> vu_d1O = vu;
          DEM_MULTIPLY_BY_SCALAR_3(vu_d1O, d_1O);
          array_1d<double, 3> pO;
          noalias(pO) = central_node.Coordinates() + vu_d1O;

          // Coordiantes of overlap ends
          array_1d<double, 3> vu_OA;
          array_1d<double, 3> vu_OB;
          array_1d<double, 3> positive_z = ZeroVector(3);
          positive_z[2] = 1.0;
          GeometryFunctions::CrossProduct(positive_z, vu, vu_OA);
          GeometryFunctions::CrossProduct(vu, positive_z, vu_OB);

          array_1d<double, 3> v_OA = vu_OA;
          array_1d<double, 3> v_OB = vu_OB;
          DEM_MULTIPLY_BY_SCALAR_3(v_OA, rc);
          DEM_MULTIPLY_BY_SCALAR_3(v_OB, rc);

          array_1d<double, 3> pA;
          array_1d<double, 3> pB;
          noalias(pA) = pO + v_OA;
          noalias(pB) = pO + v_OB;

          // Check if overlap area is inside region or cut by edges
          const double x1 = GeometryFunctions::min(pA[0], pB[0]);
          const double x2 = GeometryFunctions::max(pA[0], pB[0]);
          const double y1 = GeometryFunctions::min(pA[1], pB[1]);
          const double y2 = GeometryFunctions::max(pA[1], pB[1]);
          double ratio_in = 0.0;

          if      (x1 > xmin && x2 < xmax && y1 > ymin && y2 < ymax) ratio_in = 1.0;                     // Inside
          else if (x1 < xmin && x2 > xmin && y1 > ymin && y2 < ymax) ratio_in = (x2 - xmin) / (x2 - x1); // Left edge
          else if (x1 < xmax && x2 > xmax && y1 > ymin && y2 < ymax) ratio_in = (xmax - x1) / (x2 - x1); // Right edge
          else if (y1 < ymin && y2 > ymin && x1 > xmin && x2 < xmax) ratio_in = (y2 - ymin) / (y2 - y1); // Bottom edge
          else if (y1 < ymax && y2 > ymax && x1 > xmin && x2 < xmax) ratio_in = (ymax - y1) / (y2 - y1); // Top edge

          area_ovlp *= ratio_in;

          // Discount overlap
          area_inside -= area_ovlp;
        }

        // Add particle area to solid volume
        mRVE_VolSolidInner += area_inside;
      }

      // Compute porosity
      mRVE_PorosityInner = 1.0 - mRVE_VolSolidInner / inner_vol;
      */
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the boundary perimeter of RVE.
    double RVEWallBoundary2D::ComputeSurfaceArea(void) {
        // Vertex coordinates
        const double x1 = mVertexCoords(0,0); const double y1 = mVertexCoords(1,0);
        const double x2 = mVertexCoords(0,1); const double y2 = mVertexCoords(1,1);
        const double x3 = mVertexCoords(0,2); const double y3 = mVertexCoords(1,2);
        const double x4 = mVertexCoords(0,3); const double y4 = mVertexCoords(1,3);

        // Wall direction vectors
        const double w1x = x2-x1; const double w1y = y2-y1;
        const double w2x = x3-x2; const double w2y = y3-y2;
        const double w3x = x4-x3; const double w3y = y4-y3;
        const double w4x = x1-x4; const double w4y = y1-y4;

        // Wall lengths
        const double len1 = std::sqrt(w1x*w1x + w1y*w1y);
        const double len2 = std::sqrt(w2x*w2x + w2y*w2y);
        const double len3 = std::sqrt(w3x*w3x + w3y*w3y);
        const double len4 = std::sqrt(w4x*w4x + w4y*w4y);

        return len1 + len2 + len3 + len4;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute representative index of fabric tensor.
    double RVEWallBoundary2D::ComputeFabricIndex(Matrix fabric_tensor) {
        return fabric_tensor(0,0) - fabric_tensor(1,1);
    }

    //------------------------------------------------------------------------------------------------------------
    // Increment rose diagram with a contact in a given normal direction.
    void RVEWallBoundary2D::AddContactToRoseDiagram(std::vector<int>& rose_diagram, std::vector<double>& normal) {
        const double bin_length = 360.0 / rose_diagram.size();
        const double bin_half = bin_length / 2.0; // to centralize bin at 0 degrees
        double angle = atan2(normal[1], normal[0]) * 180.0 / Globals::Pi;
        if (angle < -bin_half) angle += 360.0;
        const int idx = std::abs((angle + bin_half) / bin_length);
        rose_diagram[idx] += 1;
    }

    //------------------------------------------------------------------------------------------------------------
    // Evaluate the uniformity of the rose diagram according to selected statistical metrics (currently, the standard deviation).
    void RVEWallBoundary2D::EvaluateRoseUniformity(void) {
        mRoseUnif      = ComputeRoseDiagramStdDev(mRoseDiagram);
        mRoseUnifInner = ComputeRoseDiagramStdDev(mRoseDiagramInner);
    }

    //------------------------------------------------------------------------------------------------------------
    // Write file headers for selected results.
    void RVEWallBoundary2D::WriteFileHeaders(void) {
        if (mFileGlobalResults.is_open()) {
            mFileGlobalResults << "R = ROW, C = COLUMN" << std::endl;
            mFileGlobalResults << "Ri -> C1: STEP | C2: TIME | C3: IS_MOVING | C4: EQUILIB_STEPS | C5-C12: [X1 Y1 X2 Y2 X3 Y3 X4 Y4] (VERTEX) | C13: #PARTICLES | C14: #PARTICLES_INN | C15: #CONTACTS | C16: #CONTACTS_INN | C17: AVG_COORD | C18: AVG_COORD_INN | C19: AVG_RADIUS_ALL | C20: VOL_TOT | C21: VOL_TOT_INN | C22: VOLD_SOLID | C23: VOLD_SOLID_INN | C24: POROSITY | C25: POROSITY_INN" << std::endl;
        }

        if (mFileParticleResults.is_open()) {
            mFileParticleResults << "R = ROW, C = COLUMN" << std::endl;
            mFileParticleResults << "Ri -> C1: STEP | C2: TIME | C3: ID | C4: RADIUS | C5-C6: [X Y] | C7: #CONTACTS_PARTICLE | C8: #CONTACTS_WALLS | C9-Cn: [FX FY] (EACH PARTICLE NEIGHBOR) | Cn+1-Cp: [FX FY] (EACH WALL NEIGHBOR)" << std::endl;
        }

        if (mFileContactResults.is_open()) {
            mFileContactResults << "R = ROW, C = COLUMN" << std::endl;
            mFileContactResults << "Ri -> C1: STEP | C2: TIME | C3-C8: X1 Y1 X2 Y2 FX FY" << std::endl;
        }

        if (mFileTensorResults.is_open()) {
            mFileTensorResults << "R = ROW, C = COLUMN, F = FABRIC TENSOR, S = STRESS TENSOR" << std::endl;
            mFileTensorResults << "Ri -> C1: STEP | C2: TIME" << std::endl;
            mFileTensorResults << "Rj -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY" << std::endl;
            mFileTensorResults << "Rk -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY" << std::endl;
            mFileTensorResults << "Rl -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS  | C7: WALL STRESS (ALL CONTACTS)" << std::endl;
            mFileTensorResults << "Rm -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS (INN CONTACTS)" << std::endl;
        }

        if (mFileRoseDiagram.is_open()) {
            mFileRoseDiagram << "R = ROW, C = COLUMN" << std::endl;
            mFileRoseDiagram << "Ri -> C1: STEP | C2: TIME" << std::endl;
            mFileRoseDiagram << "Rj -> C1-Cn: [BIN VALUES] | Cn+1: ROSE UNIFORMITY (ALL CONTACTS)" << std::endl;
            mFileRoseDiagram << "Rk -> C1-Cn: [BIN VALUES] | Cn+1: ROSE UNIFORMITY (INN CONTACTS)" << std::endl;
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Write selected results to opened files.
    void RVEWallBoundary2D::WriteResultFiles(void) {
        if (!IsTimeToPrintResults(mDemModelPart->GetProcessInfo()[TIME_STEPS])) return;

        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();
        const int time_step = r_process_info[TIME_STEPS];
        const double time   = r_process_info[TIME];

        // Define column widths
        const int wid_default = 10; // Column width for general numbers
        const int wid_float15 = 20; // Column width for fixed precision of 12 decimal places
        const int wid_float12 = 18; // Column width for fixed precision of 12 decimal places
        const int wid_float10 = 16; // Column width for fixed precision of 10 decimal places
        const int wid_float6  = 12; // Column width for fixed precision of 6 decimal places
        const int wid_float5  = 10; // Column width for fixed precision of 5 decimal places
        const int wid_float3  = 8;  // Column width for fixed precision of 3 decimal places

        if (mFileGlobalResults.is_open()) {
            mFileGlobalResults << std::setw(wid_default) << std::defaultfloat << time_step << " "
                               << std::setw(wid_default) << std::defaultfloat << time      << " ";

            mFileGlobalResults << std::setw(wid_default) << std::defaultfloat                   << mIsMoving          << " "
                               << std::setw(wid_default) << std::defaultfloat                   << mEquilibriumSteps  << " "
                               << "[ "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(0,0) << " "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(1,0) << " "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(0,1) << " "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(1,1) << " "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(0,2) << " "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(1,2) << " "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(0,3) << " "
                               << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mVertexCoords(1,3) << " "
                               << "] "
                               << std::setw(wid_default) << std::defaultfloat                   << mNumParticles      << " "
                               << std::setw(wid_default) << std::defaultfloat                   << mNumParticlesInner << " "
                               << std::setw(wid_default) << std::defaultfloat                   << mNumContacts       << " "
                               << std::setw(wid_default) << std::defaultfloat                   << mNumContactsInner  << " "
                               << std::setw(wid_float3)  << std::fixed << std::setprecision(3)  << mAvgCoordNum       << " "
                               << std::setw(wid_float3)  << std::fixed << std::setprecision(3)  << mAvgCoordNumInner  << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mAvgRadius         << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mVolTotal          << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mVolInner          << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mVolSolid          << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mVolSolidInner     << " "
                               << std::setw(wid_float5)  << std::fixed << std::setprecision(5)  << mPorosity          << " "
                               << std::setw(wid_float5)  << std::fixed << std::setprecision(5)  << mPorosityInner     << " ";

            mFileGlobalResults << std::endl; 
        }

        if (mFileParticleResults.is_open()) {
            for (unsigned int i = 0; i < mNumParticles; i++) {
                ModelPart::ElementsContainerType::iterator it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
                SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
                const int neighbors_p = particle.mNeighbourElements.size();
                const int neighbors_w = particle.mNeighbourRigidFaces.size();

                mFileParticleResults << std::setw(wid_default) << std::defaultfloat << time_step << " "
                                     << std::setw(wid_default) << std::defaultfloat << time      << " ";

                mFileParticleResults << std::setw(wid_default) << std::defaultfloat                   << particle.Id()                << " "
                                     << std::setw(wid_float12) << std::fixed << std::setprecision(12) << particle.GetRadius()         << " "
                                     << "[ "
                                     << std::setw(wid_float12) << std::fixed << std::setprecision(12) << particle.GetGeometry()[0][0] << " "
                                     << std::setw(wid_float12) << std::fixed << std::setprecision(12) << particle.GetGeometry()[0][1] << " "
                                     << "] "
                                     << std::setw(wid_default) << std::defaultfloat                   << neighbors_p                  << " "
                                     << std::setw(wid_default) << std::defaultfloat                   << neighbors_w                  << " ";

                for (unsigned int j = 0; j < neighbors_p; j++) {
                    array_1d<double, 3> force = particle.mNeighbourElasticContactForces[j];
                    mFileParticleResults << std::setw(wid_float15) << std::fixed << std::setprecision(15) << force[0] << " "
                                         << std::setw(wid_float15) << std::fixed << std::setprecision(15) << force[1] << " ";
                }
                for (unsigned int j = 0; j < neighbors_w; j++) {
                    array_1d<double, 3> force = particle.mNeighbourRigidFacesElasticContactForce[j];
                    mFileParticleResults << std::setw(wid_float15) << std::fixed << std::setprecision(15) << force[0] << " "
                                         << std::setw(wid_float15) << std::fixed << std::setprecision(15) << force[1] << " ";
                }
                mFileParticleResults << std::endl;
            }
            mFileParticleResults << std::endl;
        }

        if (mFileContactResults.is_open()) {
            for (int i = 0; i < mContactChain.size(); i+=9) {
                mFileContactResults << std::setw(wid_default) << std::defaultfloat << time_step << " "
                                    << std::setw(wid_default) << std::defaultfloat << time      << " ";

                mFileContactResults << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mContactChain[i+0] << " " // X1
                                    << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mContactChain[i+1] << " " // Y1
                                    << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mContactChain[i+3] << " " // X2
                                    << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mContactChain[i+4] << " " // Y2
                                    << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mContactChain[i+6] << " " // FX
                                    << std::setw(wid_float12) << std::fixed << std::setprecision(12) << mContactChain[i+7] << " " // FY
                                    << std::endl;
            }
            mFileContactResults << std::endl;
        }

        if (mFileTensorResults.is_open()) {
            mFileTensorResults << std::setw(wid_default) << std::defaultfloat << time_step << " "
                               << std::setw(wid_default) << std::defaultfloat << time      << " "
                               << std::endl;

            mFileTensorResults << "[ " 
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensor(0,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensor(0,1) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensor(1,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensor(1,1) << " "
                               << "] "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mFidx              << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mAnisotropy        << " "
                               << std::endl;

            mFileTensorResults << "[ " 
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,1) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,1) << " "
                               << "] "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mFidxInner              << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mAnisotropyInner        << " "
                               << std::endl;

            mFileTensorResults << "[ " 
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensor(0,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensor(0,1) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensor(1,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensor(1,1) << " "
                               << "] "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mEffStress         << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mDevStress         << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mWallStress        << " "
                               << std::endl;

            mFileTensorResults << "[ " 
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,1) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,0) << " "
                               << std::setw(wid_float10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,1) << " "
                               << "] "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mEffStressInner         << " "
                               << std::setw(wid_float6)  << std::fixed << std::setprecision(6)  << mDevStressInner         << " "
                               << std::endl;

            mFileTensorResults << std::endl;
        }

        if (mFileRoseDiagram.is_open()) {
            mFileRoseDiagram << std::setw(wid_default) << std::defaultfloat << time_step << " "
                             << std::setw(wid_default) << std::defaultfloat << time      << " "
                             << std::endl;

            mFileRoseDiagram << "[ ";
            for (int i = 0; i < mRoseDiagram.size(); i++) mFileRoseDiagram << std::setw(wid_default) << std::defaultfloat << mRoseDiagram[i] << " ";
            mFileRoseDiagram << "] ";
            mFileRoseDiagram << std::setw(wid_float6) << std::fixed << std::setprecision(6) << mRoseUnif;
            mFileRoseDiagram << std::endl;

            mFileRoseDiagram << "[ ";
            for (int i = 0; i < mRoseDiagramInner.size(); i++) mFileRoseDiagram << std::setw(wid_default) << std::defaultfloat << mRoseDiagramInner[i] << " ";
            mFileRoseDiagram << "] ";
            mFileRoseDiagram << std::setw(wid_float6) << std::fixed << std::setprecision(6) << mRoseUnifInner;
            mFileRoseDiagram << std::endl;

            mFileRoseDiagram << std::endl;
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Sort the vertices (x1,y1), (x2,y2), (x3,y3), (x4,y4) of a quadrilateral in counterclockwise order by
    // making (x1,y1) the "lower-left" vertex, ie, with its angle closest to 45° (facing northeast).
    void RVEWallBoundary2D::SortVerticesCounterClockwise(double& x1, double& x2, double& x3, double& x4, double& y1, double& y2, double& y3, double& y4) {
        // RVE Centroid
        double cx = (x1 + x2 + x3 + x4) / 4.0;
        double cy = (y1 + y2 + y3 + y4) / 4.0;

        // Store vertices in an array
        struct Point {double x, y;};
        Point points[4] = {{x1, y1}, {x2, y2}, {x3, y3}, {x4, y4}};

        // Find the "lower-left" vertex (northeast-facing angle closest to 45°)
        int lower_left_idx = 0;
        double closest_angle = std::atan2(points[0].y - cy, points[0].x - cx);
        for (int i = 1; i < 4; ++i) {
            double angle = std::atan2(points[i].y - cy, points[i].x - cx);
            if (std::abs(angle - Globals::Pi/4) < std::abs(closest_angle - Globals::Pi/4)) {
                lower_left_idx = i;
                closest_angle = angle;
            }
        }

        // Make the "lower-left" vertex as V1
        Point temp = points[lower_left_idx];
        points[lower_left_idx] = points[0];
        points[0] = temp;

        // Sort the other vertices by their angle with respect to the centroid
        std::sort(points+1, points+4, [cx,cy](const Point& a, const Point& b) {
            return std::atan2(a.y-cy, a.x-cx) < std::atan2(b.y-cy, b.x-cx);
        });

        // Assign the sorted vertices back to the original coordinates
        x1 = points[0].x; y1 = points[0].y;
        x2 = points[1].x; y2 = points[1].y;
        x3 = points[2].x; y3 = points[2].y;
        x4 = points[3].x; y4 = points[3].y;
    }

    //------------------------------------------------------------------------------------------------------------
    // Uses the shoelace formula to compute the area of a quadrilateral with vertices (x1,y1), (x2,y2), (x3,y3), (x4,y4)
    // listed in counterclockwise or clockwise order.
    double RVEWallBoundary2D::ComputeAreaFromVertices(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
        return 0.5 * std::abs(x1*y2 + x2*y3 + x3*y4 + x4*y1 - x2*y1 - x3*y2 - x4*y3 - x1*y4);
    }
    
    //------------------------------------------------------------------------------------------------------------
    // Compute the standard deviation of the bin values of a given rose diagram.
    double RVEWallBoundary2D::ComputeRoseDiagramStdDev(std::vector<int> rose_diagram) {
        // Number of bins
        int len = rose_diagram.size();
        if (len == 0) return 0.0;

        // Summation
        double sum = 0.0;
        for (int i = 0; i < len; i++) sum += rose_diagram[i];
        if (sum == 0.0) return 0.0;

        // Normalize
        std::vector<double> norm(len, 0.0);
        for (int i = 0; i < len; i++) norm[i] = rose_diagram[i] / sum;

        // Mean  (from normalized values)
        double mean = 0.0;
        for (int i = 0; i < len; i++) mean += norm[i];
        mean /= len;

        // Variance  (from normalized values)
        double var = 0.0;
        for (int i = 0; i < len; i++) var += pow(norm[i] - mean, 2.0);
        var /= len;

        // Standard deviation
        return sqrt(var);
    }
}
