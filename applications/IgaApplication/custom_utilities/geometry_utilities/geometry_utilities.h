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

#if !defined(GEOMETRY_UTILITIES_H_INCLUDED)
#define GEOMETRY_UTILITIES_H_INCLUDED

#include "includes/model_part.h"

// external library include
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
                rNode1X, rNode1Y, rNode2X, rNode2Y, rNode3X, rNode3Y, local_coordinate_x, local_coordinate_y);

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
            //for (int i = 0; i < inner_loops.size(); ++i)
            //{
            //    delaunator::Delaunator delaunators(inner_loops[i]);
            //}
 /*           }

            for (int i = 0; i < inner_delaunators.size(); ++i)
            {*/
            //    for (size_t i = 0; i < delaunators.triangles.size(); i += 3) {
            //        bool check_is_inside = IsInsideTriangle2D(
            //            rLocationX,
            //            rLocationY,
            //            inner_loops[i][2 * delaunators.triangles[i]],
            //            inner_loops[i][2 * delaunators.triangles[i] + 1],
            //            inner_loops[i][2 * delaunators.triangles[i + 1]],
            //            inner_loops[i][2 * delaunators.triangles[i + 1] + 1],
            //            inner_loops[i][2 * delaunators.triangles[i + 2]],
            //            inner_loops[i][2 * delaunators.triangles[i + 2] + 1]
            //        );
            //        if (check_is_inside)
            //            return false;
            //    }
            //}
 //           //std::vector<delaunator::Delaunator> outer_delaunators(outer_loops.size());
 //           for (int i = 0; i < outer_loops.size(); ++i)
 //           {
 //               delaunator::Delaunator delaunators(outer_loops[i]);
 //           //}
 //           //for (int i = 0; i < outer_delaunators.size(); ++i)
 //           //{
 //               for (size_t i = 0; i < delaunators.triangles.size(); i += 3) {
 //                   bool check_is_inside = IsInsideTriangle2D(
 //                       rLocationX,
 //                       rLocationY,
 //                       outer_loops[i][2 * delaunators.triangles[i]],
 //                       outer_loops[i][2 * delaunators.triangles[i] + 1],
 //                       outer_loops[i][2 * delaunators.triangles[i + 1]],
 //                       outer_loops[i][2 * delaunators.triangles[i + 1] + 1],
 //                       outer_loops[i][2 * delaunators.triangles[i + 2]],
 //                       outer_loops[i][2 * delaunators.triangles[i + 2] + 1]
 //                   );
 //                   if (check_is_inside)
 //                       return true;
 //               }
 //           }

            return false;
        }
    }



} // namespace Kratos

#endif // GEOMETRY_UTILITIES_H_INCLUDED
