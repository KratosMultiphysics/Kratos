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
#include "geometries/geometry.h"
#include "containers/array_1d.h"

namespace Kratos
{
template<int TDimension, class TPointType> class NurbsSurfaceGeometry;
class ProjectionNurbsGeometryUtilities
{
public:

    typedef array_1d<double, 3> CoordinatesArrayType;

    /*
    * @brief Returns the projection of a point onto a Nurbs curve
    *        geometry using the Newton-Rapshon iterative method
    * @param rProjectedPointLoaclCoordinates Intial guess for the Newton-Rapshon algorithm
    *        overwritten by the local coordinates of the projected point onto
    *        the Nurbs curve geometry
    * @param rPoint The point to be projected onto the Nurbs curve geometry
    *        This is overwritten by the Cartesian coordinates of the projected
    *        point in case the projection is successful
    * @param rProjectedPointGlobalCoordinates The projection onto the Nurbs curve geometry
    * @param rNurbsCurve The Nurbs curve geometry onto which the point is
    *        to be projected
    * @param MaxIterations Maximum number of iterations for the Newton-Rapshon
    *        algorithm
    * @param Accuracy Accuracy for the the Newton-Rapshon algorithm
    */
    template <class TPointType>
    static bool NewtonRaphsonCurve(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinatesCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const Geometry<TPointType>& rGeometry,
        const int MaxIterations = 20,
        const double Accuracy = 1e-6)
    {
        // Intialize variables
        double residual, delta_t;

        std::vector<array_1d<double, 3>> derivatives(3);
        array_1d<double, 3> distance_vector;

        bool projection_reset_to_boundary = false;

        // Loop over all Newton-Raphson iterations
        for (int i = 0; i < MaxIterations; ++i)
        {
            // Compute the position, the base and the acceleration vector
            rGeometry.GlobalSpaceDerivatives(
                derivatives,
                rProjectedPointLocalCoordinates,
                2);
            rProjectedPointGlobalCoordinates = derivatives[0];

            // Compute the distance vector between the point and its
            // projection on the curve
            distance_vector = rProjectedPointGlobalCoordinates - rPointGlobalCoordinatesCoordinates;
            if (norm_2(distance_vector) < Accuracy)
                return true;

            // Compute the residual
            residual = inner_prod(distance_vector, derivatives[1]);
            if (std::abs(residual) < Accuracy)
                return true;

            // Compute the increment
            delta_t = residual / (inner_prod(derivatives[2], distance_vector) + pow(norm_2(derivatives[1]), 2));
            
            // Compute the update vector and correct direction if needed
            const array_1d<double, 3> update_vector = delta_t * derivatives[1];
            const double alignment = inner_prod(update_vector, distance_vector);
            
            // If the update vector is not aligned with the distance vector, invert the sign to ensure moving towards the right point
            if (alignment > 0.0) {
                delta_t *= -1.0;
            }

            // Increment the parametric coordinate
            rProjectedPointLocalCoordinates[0] += delta_t;

            // Check if the increment is too small and if yes return true
            if (norm_2(delta_t * derivatives[1]) < Accuracy)
                return true;

            // Check if the parameter gets out of its interval of definition and if so clamp it
            // back to the boundaries
            int check = rGeometry.ClosestPointLocalToLocalSpace(
                rProjectedPointLocalCoordinates, rProjectedPointLocalCoordinates);
            if (check == 0) {
                if (projection_reset_to_boundary)
                    return false;
                else
                    projection_reset_to_boundary = true;
            }
        }

        // Return false if the Newton-Raphson iterations did not converge
        return false;
    }

    /*
    * @brief Returns the projection of a point onto a Nurbs surface
    *        geometry using the Newton-Rapshon iterative method
    * @param rProjectedPointLocalCoordinates Intial guess for the Newton-Rapshon algorithm
    *        overwritten by the local coordinates of the projected point onto
    *        the Nurbs surface geometry
    * @param rPoint The point to be projected onto the Nurbs surface geometry
    *        This is overwritten by the Cartesian coordinates of the projected
    *        point in case the projection is successful
    * @param rResult The projection onto the Nurbs surface geometry
    * @param rNurbsCurve The Nurbs curve geometry onto which the point is
    *        to be projected
    * @param MaxIterations Maximum number of iterations for the Newton-Rapshon
    *        algorithm
    * @param Accuracy Accuracy for the the Newton-Rapshon algorithm
    */
    template <int TDimension, class TPointType>
    static bool NewtonRaphsonSurface(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const NurbsSurfaceGeometry<TDimension, TPointType>& rNurbsSurface,
        const int MaxIterations = 20,
        const double Accuracy = 1e-6)
    {
        // Initialize variables
        bool is_first_row_zero, is_second_row_zero, is_first_column_zero, is_second_column_zero, is_system_invertible;
        double d_u = 0.0;
        double d_v = 0.0;
        double xi_cos, eta_cos, residual_u, residual_v, j_00, j_01, j_11, det_j;

        // Loop over all the Newton-Raphson iterations
        for (int i = 0; i < MaxIterations; i++) {

            // Compute the position, the base and the acceleration vectors
            std::vector<array_1d<double, 3>> s;
            rNurbsSurface.GlobalSpaceDerivatives(s, rProjectedPointLocalCoordinates, 2);
            rProjectedPointGlobalCoordinates = s[0];

            // Compute the distance vector
            const array_1d<double, 3> distance_vector = s[0] - rPointGlobalCoordinates;

            // Compute the distance
            const double distance = norm_2(distance_vector);
            if (distance < Accuracy)
                return true;

            // Compute the residuals along both parametric directions
            residual_u = -inner_prod(s[1], distance_vector);
            residual_v = -inner_prod(s[2], distance_vector);

            // Compute the cosine with respect to the u-parametric coordinate
            xi_cos = std::abs(residual_u) / norm_2(s[1]) / norm_2(distance_vector);

            // Compute the cosine with respect to the v-parametric coordinate
            eta_cos = std::abs(residual_v) / norm_2(s[2]) / norm_2(distance_vector);

            // Check the orthogonality condition
            if (xi_cos < Accuracy && eta_cos < Accuracy)
                return true;

            // Compute the Jacobian of the nonlinear system
            j_00 = inner_prod(s[1], s[1]) + inner_prod(s[3], distance_vector);
            j_01 = inner_prod(s[1], s[2]) + inner_prod(s[4], distance_vector);
            j_11 = inner_prod(s[2], s[2]) + inner_prod(s[5], distance_vector);

            // Check for singularities otherwise update the parametric coordinates as usual
            is_first_row_zero = false;
            if ((std::abs(j_00) < Accuracy && std::abs(j_01) < Accuracy)) {
                is_first_row_zero = true;
            }
            is_second_row_zero = false;
            if (std::abs(j_01) < Accuracy && fabs(j_11) < Accuracy) {
                is_second_row_zero = true;
            }
            is_first_column_zero = false;
            if ((std::abs(j_00) < Accuracy && std::abs(j_01) < Accuracy)) {
                is_first_column_zero = true;
            }
            is_second_column_zero = false;
            if ((std::abs(j_01) < Accuracy && std::abs(j_11) < Accuracy)) {
                is_second_column_zero = true;
            }

            // Check if the system is solvable by checking the condition of the diagonal entries
            is_system_invertible = true;
            if (is_first_row_zero || is_second_row_zero || is_first_column_zero || is_second_column_zero) {
                is_system_invertible = false;
            }

            // Solve the 2x2 linear equation system and take into account special cases where singularities occur
            if (is_system_invertible) {
                det_j = j_00 * j_11 - j_01 * j_01;
                d_u = -(residual_v * j_01 - residual_u * j_11) / det_j;
                d_v = -(residual_u * j_01 - residual_v * j_00) / det_j;
            }
            else {
                if (is_first_row_zero) {
                    d_u = residual_v / j_11;
                    d_v = 0.0;
                }
                else if (is_second_row_zero) {
                    d_u = residual_u / j_00;
                    d_v = 0.0;
                }
                else if (is_first_column_zero) {
                    d_v = (residual_u + residual_v) / (j_01 + j_11);
                    d_u = 0.0;
                }
                else if (is_second_column_zero) {
                    d_u = (residual_u + residual_v) / (j_00 + j_01);
                    d_v = 0.0;
                }
            }

            // Check if the step size is too small
            if (norm_2(d_u * s[1] + d_v * s[2]) < Accuracy)
                return true;

            // Update the parametric coordinates
            rProjectedPointLocalCoordinates[0] += d_u;
            rProjectedPointLocalCoordinates[1] += d_v;

            // Check if the parametric coordinates get out of their interval of definition
            // and if so clamp them back to their boundaries
            rNurbsSurface.DomainIntervalU().IsInside(rProjectedPointLocalCoordinates[0]);
            rNurbsSurface.DomainIntervalV().IsInside(rProjectedPointLocalCoordinates[1]);
        }

        return false;
    }
};
} // namespace Kratos

#endif // KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED
