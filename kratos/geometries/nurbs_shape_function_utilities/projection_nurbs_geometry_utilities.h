//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//                   Tobias Teschemacher
//

#if !defined(KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED)
#define KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
// #include "geometries/nurbs_curve_geometry.h"
// #include "geometries/nurbs_surface_geometry.h"

namespace Kratos
{
    template<int TDimension, class TPointType> class NurbsCurveGeometry;
    template<int TDimension, class TPointType> class NurbsSurfaceGeometry;
    namespace ProjectionNurbsGeometryUtilities
    {
        typedef array_1d<double, 3> CoordinatesArrayType;

        /*
        * @brief Returns the projection of a point onto a Nurbs curve
        *        geometry using the Newton-Rapshon iterative method
        * @param rInitialGuessParameter Intial guess for the Newton-Rapshon 
        *        algorithm
        * @param rPoint The point to be projected onto the Nurbs curve geometry
        *        This is overwritten by the Cartesian coordinates of the projected
        *        point in case the projection is successful 
        * @param rResult The projection onto the Nurbs curve geometry
        * @param rNurbsCurve The Nurbs curve geometry onto which the point is 
        *        to be projected
        * @param MaxIterations Maximum number of iterations for the Newton-Rapshon 
        *        algorithm
        * @param ModelTolerance Tolerance of the CAD model
        * @param Accuracy Accuracy for the the Newton-Rapshon algorithm
        */
        template <int TDimension, class TPointType>
        static bool NewtonRaphsonCurve(
            const double rInitialGuessParameter,
            CoordinatesArrayType& rPoint,
            CoordinatesArrayType& rResult,
            const NurbsCurveGeometry<TDimension, TPointType>& rNurbsCurve,
            const int MaxIterations = 20,
            const double ModelTolerance = 1e-6,
            const double Accuracy = 1e-6)
        {
            rResult[0] = rInitialGuessParameter;

            // Loop over all Newton-Raphson iterations
            for (int i = 0; i < MaxIterations; i++) 
            {
                // Compute the position, the base and the acceleration vector
                auto derivatives = rNurbsCurve.GlobalDerivatives(
                    rResult,
                    2);

                // Compute the distance vector between the point and its 
                // projection on the curve
                array_1d<double, 3> distance_vector = derivatives[0] - rPoint;

                // Compute the distance between the point and its projection
                // on the curve
                if (norm_2(distance_vector) < std::max(ModelTolerance,Accuracy)) {
                    KRATOS_WATCH("Inside the first check")
                    rPoint = derivatives[0];
                    return true;
                }

                // Compute the residual
                double res = inner_prod(distance_vector, derivatives[1]);
                if (std::abs(res) < Accuracy) {
                    rPoint = derivatives[0];
                    return true;
                }

                // Compute the increment
                double delta_t = res / (inner_prod(derivatives[2], distance_vector) + pow(norm_2(derivatives[1]), 2));

                // Increment the parametric coordinate
                rResult[0] -= delta_t;

                // Check if the increment is too small and if yes return true
                if (norm_2(delta_t*derivatives[1]) < Accuracy) {
                    KRATOS_WATCH("Inside the third check")
                    rPoint = derivatives[0];
                    return true;
                }

                // Check if the parameter gets out of its interval of definition and if so clamp it 
                // back to the boundaries
                rNurbsCurve.DomainInterval().IsInside(rResult[0]);
            }

            // Return false if the Newton-Raphson iterations did not converge
            return false;
    }

    template <int TDimension, class TPointType>
    static bool NewtonRaphsonSurface(
        const CoordinatesArrayType rInitialGuessParameter,
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rResult,
        const NurbsSurfaceGeometry<TDimension, TPointType>& rNurbsSurface,
        const int MaxIterations,
        const double ModelTolerance,
        const double Accuracy)
    {
        rPoint[0] = rInitialGuessParameter[0];
        rPoint[1] = rInitialGuessParameter[1];

        for (int i = 0; i < MaxIterations; i++) {
            const auto s = rNurbsSurface->GlobalDerivatives(rPoint[0], rPoint[1], 2);

            const array_1d<double, 3> distance = s[0] - rPoint;

            const double c1v = norm_2(distance);

            if (c1v < Accuracy) {
                rResult[0] = rPoint[0];
                rResult[1] = rPoint[1];

                return true;
            }

            double s1_l = norm_2(s[1]);
            double s2_l = norm_2(s[2]);

            double c2an = std::abs(inner_prod(s[1], distance));
            double c2ad = s1_l * c1v;

            double c2bn = std::abs(inner_prod(s[2], distance));
            double c2bd = s2_l * c1v;

            double c2av = c2an / c2ad;
            double c2bv = c2bn / c2bd;

            if (c2av < Accuracy && c2bv < Accuracy) {
                rResult[0] = rPoint[0];
                rResult[1] = rPoint[1];

                return true;
            }

            double f = inner_prod(s[1], distance);
            double g = inner_prod(s[2], distance);

            double J_00 = inner_prod(s[1], s[1]) + inner_prod(s[3], distance);
            double J_01 = inner_prod(s[1], s[2]) + inner_prod(s[4], distance);
            double J_11 = inner_prod(s[2], s[2]) + inner_prod(s[5], distance);

            double det_J = J_00 * J_11 - J_01 * J_01;

            double d_u = (g * J_01 - f * J_11) / det_J;
            double d_v = (f * J_01 - g * J_00) / det_J;

            rResult[0] += d_u;
            rResult[1] += d_v;
        }

        rResult[0] = rPoint[0];
        rResult[1] = rPoint[1];
        return false;
    }
    }
} // namespace Kratos

#endif // KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED
