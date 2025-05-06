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
    // TODO: THIS FUNCTION ONLY WORKS FOR SQUARE AND SHEARED RVEs WITH INCLINED LATERAL WALLS. IT SHOULD BE ADAPTED TO CONSIDER ANY SHAPE.
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
        SortVerticesCCW(x1, x2, x3, x4, y1, y2, y3, y4);

        // Set vertex coordinates
        mVertexCoords(0,0) = x1; mVertexCoords(1,0) = y1;
        mVertexCoords(0,1) = x2; mVertexCoords(1,1) = y2;
        mVertexCoords(0,2) = x3; mVertexCoords(1,2) = y3;
        mVertexCoords(0,3) = x4; mVertexCoords(1,3) = y4;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute and store the coordinates of RVE inner vertices based on original wall vertices.
    // Inner volume (area in 2D):
    //    Considered as the effective area without wall boundary effects.
    //    Applies an offset proportional to the particle radius to the RVE walls to obtain an inner quadrilateral.
    //    The optimal offset value was determined as 0.8R.
    //    This simplistic method was used in the simulations of the reference works.
    // TODO: CHECK IT FOR NON-ORTHOGONAL QUADRILATERALS!
    void RVEWallBoundary2D::SetVertexCoordinatesInner(void) {
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

        // Set inner vertex coordinates
        mVertexCoordsInner(0,0) = x1i; mVertexCoordsInner(1,0) = y1i;
        mVertexCoordsInner(0,1) = x2i; mVertexCoordsInner(1,1) = y2i;
        mVertexCoordsInner(0,2) = x3i; mVertexCoordsInner(1,2) = y3i;
        mVertexCoordsInner(0,3) = x4i; mVertexCoordsInner(1,3) = y4i;
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

                // Check for valid existing contact
                if (particle.mBallToRigidFaceStoredInfo.find(id2) == particle.mBallToRigidFaceStoredInfo.end() ||
                    particle.mBallToRigidFaceStoredInfo[id2].indentation <= 0.0)
                    continue;

                // Increment number of contacts
                mAvgCoordNum++;
                mNumContacts++;
                is_inner_particle = false;

                // Normal vector
                const double d  = r1 - particle.mBallToRigidFaceStoredInfo[id2].indentation;
                const double nx = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToRigidFaceStoredInfo[id2].local_coord_system[2][1];
                std::vector<double> normal = {nx, ny};
                std::vector<double> branch = {d * nx, d * ny};

                // Update rose diagram
                AddContactToRoseDiagram(mRoseDiagram, normal);

                // Applied wall force (normal component)
                const double fx = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[0];
                const double fy = particle.mBallToRigidFaceStoredInfo[id2].global_contact_force[1];
                std::vector<double> force = {fx, fy};
                mWallForces += std::abs(fx * nx + fy * ny);

                // Contact results
                std::vector<double> chain{x1, y1, 0.0, x1+branch[0], y1+branch[1], 0.0, fx, fy, 0.0};
                mContactChain.insert(mContactChain.end(), chain.begin(), chain.end());

                // Tensors
                for (unsigned int k = 0; k < mDim; k++) {
                    for (unsigned int l = 0; l < mDim; l++) {
                        mFabricTensor(k,l) += normal[k] * normal[l];
                        mStressTensor(k,l) += branch[k] * force[l];
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

                // Check for valid existing contact
                if (particle.mBallToBallStoredInfo.find(id2) == particle.mBallToBallStoredInfo.end()) continue;
                const double indent = particle.mBallToBallStoredInfo[id2].indentation;
                if (indent <= 0.0) continue;

                // Normal vector
                const double d  = r1 + r2 - indent;
                const double nx = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][0];
                const double ny = -particle.mBallToBallStoredInfo[id2].local_coord_system[2][1];
                std::vector<double> normal = {nx, ny};
                std::vector<double> branch = {d * nx, d * ny};

                // Check for inner contact
                const double inner_contact_len = ComputeBranchLengthInner(coords1, coords2);
                const double inner_contact_ratio = inner_contact_len / d;
                bool is_inner_contact = (inner_contact_ratio != 0.0);

                // Increment number of contacts
                mAvgCoordNum++;
                if (is_inner_particle) mAvgCoordNumInner++;

                // Update rose diagram
                AddContactToRoseDiagram(mRoseDiagram, normal);
                if (is_inner_contact) AddContactToRoseDiagram(mRoseDiagramInner, normal);

                // Unique contacts (each binary contact evaluated only once)
                if (id1 < id2) {
                    // Increment number of unique contacts
                    mNumContacts++;
                    if (is_inner_contact) mNumContactsInner++;

                    // Contact results
                    const double fx = particle.mBallToBallStoredInfo[id2].global_contact_force[0];
                    const double fy = particle.mBallToBallStoredInfo[id2].global_contact_force[1];
                    std::vector<double> force = {fx, fy};
                    std::vector<double> chain{x1, y1, 0.0, x2, y2, 0.0, fx, fy, 0.0};
                    mContactChain.insert(mContactChain.end(), chain.begin(), chain.end());

                    // Tensors
                    for (unsigned int k = 0; k < mDim; k++) {
                        for (unsigned int l = 0; l < mDim; l++) {
                            mFabricTensor(k,l) += normal[k] * normal[l];
                            mStressTensor(k,l) += branch[k] * force[l];
                            if (is_inner_contact) {
                                mFabricTensorInner(k,l) += normal[k] * normal[l];
                                mStressTensorInner(k,l) += branch[k] * force[l] * inner_contact_ratio;
                            }
                        }
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the length of a branch vector inside RVE inner volume.
    double RVEWallBoundary2D::ComputeBranchLengthInner(std::vector<double>& coords1, std::vector<double>& coords2) {
        // Branch vector end coordinates
        double x1 = coords1[0]; double y1 = coords1[1];
        double x2 = coords2[0]; double y2 = coords2[1];

        // Inner vertex coordinates
        double vx1 = mVertexCoordsInner(0,0); double vy1 = mVertexCoordsInner(1,0);
        double vx2 = mVertexCoordsInner(0,1); double vy2 = mVertexCoordsInner(1,1);
        double vx3 = mVertexCoordsInner(0,2); double vy3 = mVertexCoordsInner(1,2);
        double vx4 = mVertexCoordsInner(0,3); double vy4 = mVertexCoordsInner(1,3);
        double edges[4][4] = {{vx1,vy1,vx2,vy2}, {vx2,vy2,vx3,vy3}, {vx3,vy3,vx4,vy4}, {vx4,vy4,vx1,vy1}};

        // Check if the endpoints of the line are inside the quadrilateral
        std::array<std::pair<double, double>, 4> intersections;
        int count = 0;
        if (PointInsideQuadrilateral(x1,y1,vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4)) intersections[count++] = {x1,y1};
        if (PointInsideQuadrilateral(x2,y2,vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4)) intersections[count++] = {x2,y2};

        // Compute intersections of the line with the quadrilateral edges
        for (int i = 0; i < 4; i++) {
            double ex1 = edges[i][0], ey1 = edges[i][1], ex2 = edges[i][2], ey2 = edges[i][3];
            double denom = (x1-x2) * (ey1-ey2) - (y1-y2) * (ex1-ex2);
            if (denom == 0)
                continue; // Parallel lines, no intersection
            double t = ((x1-ex1) * (ey1-ey2) - (y1-ey1) * (ex1-ex2)) / denom;
            double u = ((x1-ex1) * (y1-y2)   - (y1-ey1) * (x1-x2))   / denom;
            if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
                double ix = x1 + t * (x2-x1);
                double iy = y1 + t * (y2-y1);
                intersections[count++] = {ix,iy};
            }
        }

        // Compute the maximum segment length inside the quadrilateral
        if (count < 2) { // Zero length
            return 0.0;
        }
        else if (count == 2) {
            return std::hypot(intersections[1].first - intersections[0].first, intersections[1].second - intersections[0].second);
        }
        else {
            double max_len = 0.0;
            for (int i = 0; i < count; i++) {
                for (int j = i+1; j < count; j++) {
                    double len = std::hypot(intersections[j].first - intersections[i].first, intersections[j].second - intersections[i].second);
                    max_len = std::max(max_len, len);
                }
            }
            return max_len;
        }  
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute particle volume, currently considering circle area in 2D.
    double RVEWallBoundary2D::ComputeVolumeParticle(SphericParticle& particle) {
        return Globals::Pi * particle.GetRadius() * particle.GetRadius();
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute particle volume inside RVE inner volume (areas in 2D).
    double RVEWallBoundary2D::ComputeVolumeParticleInner(SphericParticle& particle) {
        // Particle properties
        double r = particle.GetRadius();
        double x = particle.GetGeometry()[0][0];
        double y = particle.GetGeometry()[0][1];

        // Inner vertex coordinates
        double x1 = mVertexCoordsInner(0,0); double y1 = mVertexCoordsInner(1,0);
        double x2 = mVertexCoordsInner(0,1); double y2 = mVertexCoordsInner(1,1);
        double x3 = mVertexCoordsInner(0,2); double y3 = mVertexCoordsInner(1,2);
        double x4 = mVertexCoordsInner(0,3); double y4 = mVertexCoordsInner(1,3);
        double xv[] = {x1,x2,x3,x4};         double yv[] = {y1,y2,y3,y4};

        // Distances from circle center to each wall
        double de1 = ComputeDistancePointToSegment(x,y,x1,y1,x2,y2);
        double de2 = ComputeDistancePointToSegment(x,y,x2,y2,x3,y3);
        double de3 = ComputeDistancePointToSegment(x,y,x3,y3,x4,y4);
        double de4 = ComputeDistancePointToSegment(x,y,x4,y4,x1,y1);

        // Distances from circle center to each vertex
        double dv1 = std::sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
        double dv2 = std::sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2));
        double dv3 = std::sqrt((x-x3)*(x-x3)+(y-y3)*(y-y3));
        double dv4 = std::sqrt((x-x4)*(x-x4)+(y-y4)*(y-y4));

        // Check is circle center is inside quadrilateral
        bool is_inside = PointInsideQuadrilateral(x,y,x1,y1,x2,y2,x3,y3,x4,y4);

        // Compute particle volume for each case of interaction with quadrilateral edges/vertices
        double vol = 0.0; // default: circle completely outside

        if (de1 >= r && de2 >= r && de3 >= r && de4 >= r) {
            // Circle completely inside
            if (is_inside) vol = ComputeVolumeParticle(particle);
        }
        else {
            // Circle cut by one or more edges, but with no vertex inside
            if (dv1 >= r && dv2 >= r && dv3 >= r && dv4 >= r) {
                double a1 = (de1 < r) ? ComputeAreaCircleSector(r,de1) : 0.0;
                double a2 = (de2 < r) ? ComputeAreaCircleSector(r,de2) : 0.0;
                double a3 = (de3 < r) ? ComputeAreaCircleSector(r,de3) : 0.0;
                double a4 = (de4 < r) ? ComputeAreaCircleSector(r,de4) : 0.0;
                if (is_inside)
                    vol = ComputeVolumeParticle(particle) - (a1+a2+a3+a4);
                else
                    vol = std::max({a1,a2,a3,a4}) - std::max(std::min({a1,a2,a3,a4}), 0.0);
            }
            // Circle with vertex inside (no more than 1)
            else {
                double dv[] = {dv1,dv2,dv3,dv4};
                for (int i = 0; i < 4; i++) {
                    if (dv[i] < r) {
                        // Vertex (common end to the two edges)
                        double xa = xv[i];
                        double ya = yv[i];

                        // Other edge end coordinates
                        double xe1, ye1;
                        double xe2, ye2;
                        if (i == 1) {
                            xe1 = xv[2]; ye1 = yv[2];
                            xe2 = xv[4]; ye2 = yv[4];
                        }
                        else if (i == 4) {
                            xe1 = xv[1]; ye1 = yv[1];
                            xe2 = xv[3]; ye2 = yv[3];
                        }
                        else {
                            xe1 = xv[i-1]; ye1 = yv[i-1];
                            xe2 = xv[i+1]; ye2 = yv[i+1];
                        }

                        // Intersection between edges and circle
                        double xb, yb;
                        double xc, yc;
                        ComputeIntersectionCircleSegment(r,x,y,xa,ya,xe1,ye1,xb,yb);
                        ComputeIntersectionCircleSegment(r,x,y,xa,ya,xe2,ye2,xc,yc);

                        // Area of triangle formed by the vertex and intersection points
                        double area_triangle = ComputeAreaTriangle(xa,xb,xc,ya,yb,yc);

                        // Area of circle sector formed by the intersection points and circumference
                        double xt = (xb+xc)/2.0, yt = (yb+yc)/2.0;            // Middle coordinates of the chord between intersection points (xb,yb)->(xc,yc)
                        double d1 = std::sqrt((xt-x)*(xt-x) + (yt-y)*(yt-y)); // Distance: circle center - chord middle
                        double d2 = std::sqrt((xa-x)*(xa-x) + (ya-y)*(ya-y)); // Distance: circle center - vertex
                        double dot = (xa-x)*(xt-x) + (ya-y)*(yt-y);           // Dot product between colinear vectors (circle center)->(vertex) and (circle center)->(chord middle)

                        double area_sector;
                        if (dot > 0.0 && d2 > d1)
                            area_sector = ComputeVolumeParticle(particle) - ComputeAreaCircleSector(r,d1);
                        else
                            area_sector = ComputeAreaCircleSector(r,d1);

                        // Total area delimited by the edges
                        vol = area_triangle + area_sector;
                        break;
                    }
                }
            }
        }
        return vol;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the total area of RVE.
    double RVEWallBoundary2D::ComputeVolumeRVE(void) {
        const double x1 = mVertexCoords(0,0); const double y1 = mVertexCoords(1,0);
        const double x2 = mVertexCoords(0,1); const double y2 = mVertexCoords(1,1);
        const double x3 = mVertexCoords(0,2); const double y3 = mVertexCoords(1,2);
        const double x4 = mVertexCoords(0,3); const double y4 = mVertexCoords(1,3);
        return ComputeAreaQuadrilateral(x1, x2, x3, x4, y1, y2, y3, y4);
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the inner volume (area in 2D) of RVE.
    double RVEWallBoundary2D::ComputeVolumeRVEInner(void) {
        double x1 = mVertexCoordsInner(0,0); double y1 = mVertexCoordsInner(1,0);
        double x2 = mVertexCoordsInner(0,1); double y2 = mVertexCoordsInner(1,1);
        double x3 = mVertexCoordsInner(0,2); double y3 = mVertexCoordsInner(1,2);
        double x4 = mVertexCoordsInner(0,3); double y4 = mVertexCoordsInner(1,3);
        return ComputeAreaQuadrilateral(x1, x2, x3, x4, y1, y2, y3, y4);
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
    void RVEWallBoundary2D::WriteFileHeadersGlobalResults(void) {
        mFileGlobalResults << "R = ROW, C = COLUMN" << std::endl;
        mFileGlobalResults << "Ri -> C1: STEP | C2: TIME | C3: IS_MOVING | C4: EQUILIB_STEPS | C5-C12: [X1 Y1 X2 Y2 X3 Y3 X4 Y4] (VERTEX) | C13: #PARTICLES | C14: #PARTICLES_INN | C15: #CONTACTS | C16: #CONTACTS_INN | C17: AVG_COORD | C18: AVG_COORD_INN | C19: AVG_RADIUS_ALL | C20: VOL_TOT | C21: VOL_TOT_INN | C22: VOLD_SOLID | C23: VOLD_SOLID_INN | C24: POROSITY | C25: POROSITY_INN" << std::endl;
    }
    void RVEWallBoundary2D::WriteFileHeadersParticleResults(void) {
        mFileParticleResults << "R = ROW, C = COLUMN" << std::endl;
        mFileParticleResults << "Ri -> C1: STEP | C2: TIME | C3: ID | C4: RADIUS | C5-C6: [X Y] | C7: #CONTACTS_PARTICLE | C8: #CONTACTS_WALLS | C9-Cn: [FX FY] (EACH PARTICLE NEIGHBOR) | Cn+1-Cp: [FX FY] (EACH WALL NEIGHBOR)" << std::endl;
    }
    void RVEWallBoundary2D::WriteFileHeadersContactResults(void) {
        mFileContactResults << "R = ROW, C = COLUMN" << std::endl;
        mFileContactResults << "Ri -> C1: STEP | C2: TIME | C3-C8: X1 Y1 X2 Y2 FX FY" << std::endl;
    }
    void RVEWallBoundary2D::WriteFileHeadersTensorResults(void) {
        mFileTensorResults << "R = ROW, C = COLUMN, F = FABRIC TENSOR, S = STRESS TENSOR" << std::endl;
        mFileTensorResults << "Ri -> C1: STEP | C2: TIME" << std::endl;
        mFileTensorResults << "Rj -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Rk -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY (INN CONTACTS)" << std::endl;
        mFileTensorResults << "Rl -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS  | C7: WALL STRESS (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Rm -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS (INN CONTACTS)" << std::endl;
    }
    void RVEWallBoundary2D::WriteFileHeadersRoseDiagram(void) {
        mFileRoseDiagram << "R = ROW, C = COLUMN" << std::endl;
        mFileRoseDiagram << "Ri -> C1: STEP | C2: TIME" << std::endl;
        mFileRoseDiagram << "Rj -> C1-Cn: [BIN VALUES] | Cn+1: ROSE UNIFORMITY (ALL CONTACTS)" << std::endl;
        mFileRoseDiagram << "Rk -> C1-Cn: [BIN VALUES] | Cn+1: ROSE UNIFORMITY (INN CONTACTS)" << std::endl;
    }

    //------------------------------------------------------------------------------------------------------------
    // Write selected results to opened files.
    void RVEWallBoundary2D::WriteResultFilesGlobalResults(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        mFileGlobalResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME_STEPS] << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME]       << " ";

        mFileGlobalResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat << mIsMoving         << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat << mEquilibriumSteps << " "
                           << "[ "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(0,0) << " "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(1,0) << " "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(0,1) << " "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(1,1) << " "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(0,2) << " "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(1,2) << " "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(0,3) << " "
                           << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mVertexCoords(1,3) << " "
                           << "] "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat                   << mNumParticles      << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat                   << mNumParticlesInner << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat                   << mNumContacts       << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat                   << mNumContactsInner  << " "
                           << std::setw(WIDTH_FLOAT03) << std::fixed << std::setprecision(3)  << mAvgCoordNum       << " "
                           << std::setw(WIDTH_FLOAT03) << std::fixed << std::setprecision(3)  << mAvgCoordNumInner  << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mAvgRadius         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mVolTotal          << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mVolInner          << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mVolSolid          << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mVolSolidInner     << " "
                           << std::setw(WIDTH_FLOAT05) << std::fixed << std::setprecision(5)  << mPorosity          << " "
                           << std::setw(WIDTH_FLOAT05) << std::fixed << std::setprecision(5)  << mPorosityInner     << " ";

        mFileGlobalResults << std::endl; 
    }
    void RVEWallBoundary2D::WriteResultFilesParticleResults(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        for (unsigned int i = 0; i < mNumParticles; i++) {
            ModelPart::ElementsContainerType::iterator it = mDemModelPart->GetCommunicator().LocalMesh().Elements().ptr_begin() + i;
            SphericParticle& particle = dynamic_cast<SphericParticle&>(*it);
            const int neighbors_p = particle.mNeighbourElements.size();
            const int neighbors_w = particle.mNeighbourRigidFaces.size();

            mFileParticleResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME_STEPS] << " "
                                 << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME]       << " ";

            mFileParticleResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat                   << particle.Id()        << " "
                                 << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << particle.GetRadius() << " "
                                 << "[ "
                                 << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << particle.GetGeometry()[0][0] << " "
                                 << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << particle.GetGeometry()[0][1] << " "
                                 << "] "
                                 << std::setw(WIDTH_DEFAULT) << std::defaultfloat << neighbors_p << " "
                                 << std::setw(WIDTH_DEFAULT) << std::defaultfloat << neighbors_w << " ";

            for (unsigned int j = 0; j < neighbors_p; j++) {
                array_1d<double, 3> force = particle.mNeighbourElasticContactForces[j];
                mFileParticleResults << std::setw(WIDTH_FLOAT15) << std::fixed << std::setprecision(15) << force[0] << " "
                                     << std::setw(WIDTH_FLOAT15) << std::fixed << std::setprecision(15) << force[1] << " ";
            }
            for (unsigned int j = 0; j < neighbors_w; j++) {
                array_1d<double, 3> force = particle.mNeighbourRigidFacesElasticContactForce[j];
                mFileParticleResults << std::setw(WIDTH_FLOAT15) << std::fixed << std::setprecision(15) << force[0] << " "
                                     << std::setw(WIDTH_FLOAT15) << std::fixed << std::setprecision(15) << force[1] << " ";
            }
            mFileParticleResults << std::endl;
        }
        mFileParticleResults << std::endl;
    }
    void RVEWallBoundary2D::WriteResultFilesContactResults(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        for (int i = 0; i < mContactChain.size(); i+=9) {
            mFileContactResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME_STEPS] << " "
                                << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME]       << " ";

            mFileContactResults << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mContactChain[i+0] << " " // X1
                                << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mContactChain[i+1] << " " // Y1
                                << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mContactChain[i+3] << " " // X2
                                << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mContactChain[i+4] << " " // Y2
                                << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mContactChain[i+6] << " " // FX
                                << std::setw(WIDTH_FLOAT12) << std::fixed << std::setprecision(12) << mContactChain[i+7] << " " // FY
                                << std::endl;
        }
        mFileContactResults << std::endl;
    }
    void RVEWallBoundary2D::WriteResultFilesTensorResults(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        mFileTensorResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME_STEPS] << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME]       << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mFidx              << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mAnisotropy        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mFidxInner              << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mAnisotropyInner        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mEffStress         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mDevStress         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mWallStress        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mEffStressInner         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mDevStressInner         << " "
                           << std::endl;

        mFileTensorResults << std::endl;
    }
    void RVEWallBoundary2D::WriteResultFilesRoseDiagram(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        mFileRoseDiagram << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME_STEPS] << " "
                         << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME]       << " "
                         << std::endl;

        mFileRoseDiagram << "[ ";
        for (int i = 0; i < mRoseDiagram.size(); i++) mFileRoseDiagram << std::setw(WIDTH_DEFAULT) << std::defaultfloat << mRoseDiagram[i] << " ";
        mFileRoseDiagram << "] ";
        mFileRoseDiagram << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6) << mRoseUnif;
        mFileRoseDiagram << std::endl;

        mFileRoseDiagram << "[ ";
        for (int i = 0; i < mRoseDiagramInner.size(); i++) mFileRoseDiagram << std::setw(WIDTH_DEFAULT) << std::defaultfloat << mRoseDiagramInner[i] << " ";
        mFileRoseDiagram << "] ";
        mFileRoseDiagram << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6) << mRoseUnifInner;
        mFileRoseDiagram << std::endl;

        mFileRoseDiagram << std::endl;
    }

    //------------------------------------------------------------------------------------------------------------
    // Sort the vertices (x1,y1), (x2,y2), (x3,y3), (x4,y4) of a quadrilateral in counterclockwise order by
    // making (x1,y1) the "lower-left" vertex, ie, with its angle closest to 45° (facing northeast).
    void RVEWallBoundary2D::SortVerticesCCW(double& x1, double& x2, double& x3, double& x4, double& y1, double& y2, double& y3, double& y4) {
        // RVE Centroid
        double cx = (x1 + x2 + x3 + x4) / 4.0;
        double cy = (y1 + y2 + y3 + y4) / 4.0;

        // Store vertices in an array
        struct Point {double x, y;};
        Point points[4] = {{x1, y1}, {x2, y2}, {x3, y3}, {x4, y4}};

        // Find the "lower-left" vertex (northeast-facing angle closest to 45°)
        int lower_left_idx = 0;
        double closest_angle = std::atan2(cy - points[0].y, cx - points[0].x);
        for (int i = 1; i < 4; ++i) {
            double angle = std::atan2(cy - points[i].y, cx - points[i].x);
            if (std::abs(angle - Globals::Pi/4) < std::abs(closest_angle - Globals::Pi/4)) {
                lower_left_idx = i;
                closest_angle = angle;
            }
        }

        // Make the "lower-left" vertex as V1
        Point temp = points[lower_left_idx];
        points[lower_left_idx] = points[0];
        points[0] = temp;

        // Sort the other vertices by their angle with respect to the "lower-left" vertex
        std::sort(points+1, points+4, [&](const Point& p1, const Point& p2) {
            return atan2(p1.y - points[0].y, p1.x - points[0].x) < atan2(p2.y - points[0].y, p2.x - points[0].x);
        });

        // Assign the sorted vertices back to the original coordinates
        x1 = points[0].x; y1 = points[0].y;
        x2 = points[1].x; y2 = points[1].y;
        x3 = points[2].x; y3 = points[2].y;
        x4 = points[3].x; y4 = points[3].y;
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the area of a quadrilateral with vertices (x1,y1), (x2,y2), (x3,y3), (x4,y4)
    // listed in counterclockwise or clockwise order, using the shoelace formula.
    double RVEWallBoundary2D::ComputeAreaQuadrilateral(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
        return 0.5 * std::abs(x1*y2 + x2*y3 + x3*y4 + x4*y1 - x2*y1 - x3*y2 - x4*y3 - x1*y4);
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the area of a triangle with vertices (x1,y1), (x2,y2), (x3,y3), listed in counterclockwise or clockwise order.
    double RVEWallBoundary2D::ComputeAreaTriangle(double x1, double x2, double x3, double y1, double y2, double y3) {
        return 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the area of a sector of a circle with radius r cut by a line whose distance from the circle center is d.
    double RVEWallBoundary2D::ComputeAreaCircleSector(double r, double d) {
        return (r*r * acos(d/r)) - (d * sqrt(r*r - d*d));
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute the normal distance of a point (x,y) to a line segment edge with vertices (xa,ya)-(xb,yb).
    // If the point projection falls outside the segment, it returns the closest end point (xa,ya) or (xb,yb).
    double RVEWallBoundary2D::ComputeDistancePointToSegment(double x, double y, double xa, double ya, double xb, double yb) {
        double dx = xb - xa;
        double dy = yb - ya;
        double l2 = dx * dx + dy * dy;

        // Avoid division by zero in case of degenerate segment
        if (l2 == 0.0) return std::sqrt((x - xa) * (x - xa) + (y - ya) * (y - ya));

        // Compute projection parameter t and clamp t to the segment range [0,1]
        double t = ((x - xa) * dx + (y - ya) * dy) / l2;
        t = std::max(0.0, std::min(1.0, t));

        // Compute closest point on the segment
        double near_x = xa + t * dx;
        double near_y = ya + t * dy;

        // Compute Euclidean distance
        return std::sqrt((x - near_x) * (x - near_x) + (y - near_y) * (y - near_y));
    }

    //------------------------------------------------------------------------------------------------------------
    // Compute intersection coordinates of a circle with radius r and center (x,y) and a line segment with vertices (xa,ya)-(xb,yb).
    // It assumes that (xa,ya) is inside the circle and (xb,yb) is outside.
    void RVEWallBoundary2D::ComputeIntersectionCircleSegment(double r, double x, double y, double xa, double ya, double xb, double yb, double& xc, double& yc) {
        double dx = xb - xa;
        double dy = yb - ya;
        double A = dx*dx + dy*dy;
        double B = 2.0 * ((xa-x)*dx + (ya-y)*dy);
        double C = (xa-x) * (xa-x) + (ya-y) * (ya-y) - r*r;
        double t = (-B - std::sqrt(B*B-4*A*C)) / (2*A);
        xc = xa + t * dx;
        yc = ya + t * dy;
    }

    //------------------------------------------------------------------------------------------------------------
    // Check if a point (x,y) is inside a quadrilateral with vertices (x1,y1), (x2,y2), (x3,y3), (x4,y4),
    // listed in a consistent order: clockwise or counterclockwise.
    bool RVEWallBoundary2D::PointInsideQuadrilateral(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
        int count = 0;
        double vertices[4][2] = {{x1,y1},{x2,y2},{x3,y3},{x4,y4}};
        for (int i = 0; i < 4; i++) {
            double xa = vertices[i][0];
            double ya = vertices[i][1];
            double xb = vertices[(i+1)%4][0];
            double yb = vertices[(i+1)%4][1];
            if (((ya>y) != (yb>y)) && ( x < (xb-xa)*(y-ya)/(yb-ya)+xa)) count++;
        }
        return count % 2 == 1;
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
