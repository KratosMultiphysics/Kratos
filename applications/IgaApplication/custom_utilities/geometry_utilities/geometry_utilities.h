//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

#if !defined(INSIDE_OUTSIDE_GEOMETRY_UTILITIES_H_INCLUDED)
#define INSIDE_OUTSIDE_GEOMETRY_UTILITIES_H_INCLUDED

// System includes
#include <limits>
#include <algorithm>
#include <cmath>

// External includes
#include <delaunator-cpp/delaunator.hpp>

namespace Kratos
{

    namespace GeometryUtilities
    {
        static void GetLocalCoordinatesTriangle2D(
            const double& rLocationX,
            const double& rLocationY,
            const double& rNode1X,
            const double& rNode1Y,
            const double& rNode2X,
            const double& rNode2Y,
            const double& rNode3X,
            const double& rNode3Y,
            double& rLocalCoordinateX,
            double& rLocalCoordinateY
        )
        {
            double J_00 = rNode2X - rNode1X;
            double J_01 = rNode3X - rNode1X;
            double J_10 = rNode2Y - rNode1Y;
            double J_11 = rNode3Y - rNode1Y;
            const double det_J = J_00 * J_11 - J_01 * J_10;

            // Compute eta and xi
            const double eta = (J_10*(rNode1X - rLocationX) +
                J_00 * (rLocationY - rNode1Y)) / det_J;
            const double xi = (J_11*(rLocationX - rNode1X) +
                J_01 * (rNode1Y - rLocationY)) / det_J;

            array_1d<double, 2> new_location;
            rLocalCoordinateX = xi;
            rLocalCoordinateY = eta;
        }

        static bool IsInsideTriangle2D(
            const double& rLocationX,
            const double& rLocationY,
            const double& rNode1X,
            const double& rNode1Y,
            const double& rNode2X,
            const double& rNode2Y,
            const double& rNode3X,
            const double& rNode3Y,
            const double Tolerance = std::numeric_limits<double>::epsilon()
        )
        {
            double local_coordinate_x = 0.0;
            double local_coordinate_y = 0.0;
            GetLocalCoordinatesTriangle2D(rLocationX, rLocationY,
                rNode1X, rNode1Y, rNode2X, rNode2Y, rNode3X, rNode3Y,
                local_coordinate_x, local_coordinate_y);

            if ((local_coordinate_x >= (0.0 - Tolerance)) && (local_coordinate_x <= (1.0 + Tolerance)))
            {
                if ((local_coordinate_y >= (0.0 - Tolerance)) && (local_coordinate_y <= (1.0 + Tolerance)))
                {
                    if ((local_coordinate_x + local_coordinate_y) <= (1.0 + Tolerance))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        static bool IsInside2D(
            const double& rLocationX,
            const double& rLocationY,
            const std::vector<std::vector<double>>& outer_loops,
            const std::vector<std::vector<double>>& inner_loops
        )
        {
            for (int loop_i = 0; loop_i < inner_loops.size(); ++loop_i)
            {
                delaunator::Delaunator delaunators(inner_loops[loop_i]);

                for (size_t i = 0; i < delaunators.triangles.size(); i += 3) {
                    bool check_is_inside = IsInsideTriangle2D(
                        rLocationX,
                        rLocationY,
                        inner_loops[loop_i][2 * delaunators.triangles[i]],
                        inner_loops[loop_i][2 * delaunators.triangles[i] + 1],
                        inner_loops[loop_i][2 * delaunators.triangles[i + 1]],
                        inner_loops[loop_i][2 * delaunators.triangles[i + 1] + 1],
                        inner_loops[loop_i][2 * delaunators.triangles[i + 2]],
                        inner_loops[loop_i][2 * delaunators.triangles[i + 2] + 1]
                    );
                    if (check_is_inside)
                        return false;
                }
            }

            for (int loop_i = 0; loop_i < outer_loops.size(); ++loop_i)
            {
                delaunator::Delaunator delaunators(outer_loops[loop_i]);

                for (size_t i = 0; i < delaunators.triangles.size(); i += 3) {
                    bool check_is_inside = IsInsideTriangle2D(
                        rLocationX,
                        rLocationY,
                        outer_loops[loop_i][2 * delaunators.triangles[i]],
                        outer_loops[loop_i][2 * delaunators.triangles[i] + 1],
                        outer_loops[loop_i][2 * delaunators.triangles[i + 1]],
                        outer_loops[loop_i][2 * delaunators.triangles[i + 1] + 1],
                        outer_loops[loop_i][2 * delaunators.triangles[i + 2]],
                        outer_loops[loop_i][2 * delaunators.triangles[i + 2] + 1]
                    );

                    if (check_is_inside)
                        return true;
                }
            }

            return false;
        }



        static int CrossProductTest(
            const double& rLocationX,
            const double& rLocationY,
            const double& rAX,
            const double& rAY,
            const double& rBX,
            const double& rBY)
        {
            double AX = rAX;
            double AY = rAY;
            double BX = rBX;
            double BY = rBY;

            if (rAY > rBY)
            {
                AX = rBX;
                AY = rBY;
                BX = rAX;
                BY = rAY;
            }

            if (rLocationY <= AY || rLocationY > BY)
                return 1;

            double delta = (AX - rLocationX) * (BY - rLocationY) - (AY - rLocationY) * (BX - rLocationX);

            if (delta > 0)
                return 1;
            if (delta < 0)
                return -1;
            else
                return 0;
        }

        // -1 = outside
        //  0 = on the boundary
        //  1 = inside
        static int IsInsidePolygon2D(
            const double& rLocationX,
            const double& rLocationY,
            const std::vector<std::vector<double>>& outer_loops,
            const std::vector<std::vector<double>>& inner_loops
        )
        {
            int t;
            for (int loop_i = 0; loop_i < inner_loops.size(); ++loop_i)
            {
                t = -1;
                for (int index_node = 0; index_node < inner_loops[loop_i].size() - 2; index_node+=2)
                {
                    t = t * CrossProductTest(
                        rLocationX,
                        rLocationY,
                        inner_loops[loop_i][index_node],
                        inner_loops[loop_i][index_node + 1],
                        inner_loops[loop_i][index_node + 2],
                        inner_loops[loop_i][index_node + 3]);

                    if (t == 0)
                        return 0;
                }
                if (inner_loops[loop_i][inner_loops[loop_i].size() - 2] != inner_loops[loop_i][0]
                    && inner_loops[loop_i][inner_loops[loop_i].size() - 1] != inner_loops[loop_i][1])
                {
                    t = t * CrossProductTest(
                        rLocationX,
                        rLocationY,
                        inner_loops[loop_i][inner_loops[loop_i].size() - 2],
                        inner_loops[loop_i][inner_loops[loop_i].size() - 1],
                        inner_loops[loop_i][0],
                        inner_loops[loop_i][1]);
                }
                if (t > 0)
                    return -1;
            }

            for (int loop_i = 0; loop_i < outer_loops.size(); ++loop_i)
            {
                t = -1;
                for (int index_node = 0; index_node < outer_loops[loop_i].size() - 2; index_node += 2)
                {
                    t = t * CrossProductTest(
                        rLocationX,
                        rLocationY,
                        outer_loops[loop_i][index_node],
                        outer_loops[loop_i][index_node + 1],
                        outer_loops[loop_i][index_node + 2],
                        outer_loops[loop_i][index_node + 3]);

                    if (t == 0)
                        return 0;
                }
                if (outer_loops[loop_i][outer_loops[loop_i].size() - 2] != outer_loops[loop_i][0]
                    && outer_loops[loop_i][outer_loops[loop_i].size() - 1] != outer_loops[loop_i][1])
                {
                    t = t * CrossProductTest(
                        rLocationX,
                        rLocationY,
                        outer_loops[loop_i][outer_loops[loop_i].size() - 2],
                        outer_loops[loop_i][outer_loops[loop_i].size() - 1],
                        outer_loops[loop_i][0],
                        outer_loops[loop_i][1]);
                }
                if (t > 0)
                    return 1;
            }
            return t;
        }
    }

} // namespace Kratos

#endif // INSIDE_OUTSIDE_GEOMETRY_UTILITIES_H_INCLUDED
