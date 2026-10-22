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

    enum class ProjectionAlgorithm
    {
        NewtonRaphson,
        LevenbergMarquardt
    };

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
    * @brief Returns the closest-point projection of a point onto a Nurbs curve
    *        geometry using a Levenberg-Marquardt iterative method.
    * @param rProjectedPointLocalCoordinates Initial guess for the iterative
    *        algorithm. This is overwritten by the local coordinate of the
    *        projected point onto the Nurbs curve geometry.
    * @param rPointGlobalCoordinates The global coordinates of the point to be
    *        projected onto the Nurbs curve geometry.
    * @param rProjectedPointGlobalCoordinates The global coordinates of the
    *        projected point on the Nurbs curve geometry. This is overwritten
    *        in case the projection is successful.
    * @param rGeometry The Nurbs curve geometry onto which the point is to be
    *        projected.
    * @param MaxIterations Maximum number of iterations for the
    *        Levenberg-Marquardt algorithm.
    * @param Accuracy Accuracy used for the distance, orthogonality, and step-size
    *        convergence checks.
    */
    template <class TPointType>
    static bool LevenbergMarquardtCurve(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const Geometry<TPointType>& rGeometry,
        const int MaxIterations = 30,
        const double Accuracy = 1e-6)
    {
        double lambda = 1e-8;
        const double lambda_factor = 10.0;
        const double max_lambda = 1e12;

        bool projection_reset_to_boundary = false;

        for (int i = 0; i < MaxIterations; ++i) {

            // Compute position and first derivative only
            std::vector<array_1d<double, 3>> derivatives(2);

            rGeometry.GlobalSpaceDerivatives(derivatives, rProjectedPointLocalCoordinates, 1);

            rProjectedPointGlobalCoordinates = derivatives[0];

            const array_1d<double, 3> distance_vector = rProjectedPointGlobalCoordinates - rPointGlobalCoordinates;

            const double distance = norm_2(distance_vector);

            if (distance < Accuracy) {
                return true;
            }

            const double tangent_norm = norm_2(derivatives[1]);

            if (tangent_norm < Accuracy) {
                return false;
            }

            // Gradient of 1/2 ||C(t)-x||^2
            const double gradient = inner_prod(derivatives[1], distance_vector);

            const double tangent_cos = std::abs(gradient) / (tangent_norm * distance);

            // Orthogonality condition
            if (tangent_cos < Accuracy) {
                return true;
            }

            // 1D Levenberg-Marquardt system:
            // (C_t · C_t + lambda) delta_t = - C_t · (C - x)
            const double a = inner_prod(derivatives[1], derivatives[1]);

            bool step_accepted = false;

            for (int lm_it = 0; lm_it < 10; ++lm_it) {

                const double denominator = a + lambda;

                if (std::abs(denominator) < Accuracy) {
                    lambda *= lambda_factor;
                    continue;
                }

                const double delta_t = -gradient / denominator;

                if (norm_2(delta_t * derivatives[1]) < Accuracy) {
                    return tangent_cos < Accuracy;
                }

                CoordinatesArrayType trial_local_coordinates = rProjectedPointLocalCoordinates;

                trial_local_coordinates[0] += delta_t;

                int check = rGeometry.ClosestPointLocalToLocalSpace(
                    trial_local_coordinates,
                    trial_local_coordinates);

                if (check == 0) {
                    if (projection_reset_to_boundary) {
                        return false;
                    } else {
                        projection_reset_to_boundary = true;
                    }
                }

                std::vector<array_1d<double, 3>> trial_derivatives(1);

                rGeometry.GlobalSpaceDerivatives(trial_derivatives, trial_local_coordinates, 0);

                const array_1d<double, 3> trial_distance_vector = trial_derivatives[0] - rPointGlobalCoordinates;

                const double trial_distance = norm_2(trial_distance_vector);

                if (trial_distance < distance) {
                    rProjectedPointLocalCoordinates = trial_local_coordinates;
                    rProjectedPointGlobalCoordinates = trial_derivatives[0];

                    lambda /= lambda_factor;
                    lambda = std::max(lambda, 1e-14);

                    step_accepted = true;
                    break;
                } else {
                    lambda *= lambda_factor;

                    if (lambda > max_lambda) {
                        return false;
                    }
                }
            }

            if (!step_accepted) {
                return false;
            }
        }

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

    /*
    * @brief Returns the closest-point projection of a point onto a Nurbs surface
    *        geometry using a Levenberg-Marquardt iterative method.
    * @param rProjectedPointLocalCoordinates Initial guess for the iterative algorithm.
    *        This is overwritten by the local coordinates of the projected point
    *        onto the Nurbs surface geometry.
    * @param rPointGlobalCoordinates The global coordinates of the point to be
    *        projected onto the Nurbs surface geometry.
    * @param rProjectedPointGlobalCoordinates The global coordinates of the
    *        projected point on the Nurbs surface geometry. This is overwritten
    *        in case the projection is successful.
    * @param rNurbsSurface The Nurbs surface geometry onto which the point is
    *        to be projected.
    * @param MaxIterations Maximum number of iterations for the
    *        Levenberg-Marquardt algorithm.
    * @param Accuracy Accuracy used for the distance, orthogonality, and step-size
    *        convergence checks.
    */
    template <int TDimension, class TPointType>
    static bool LevenbergMarquardtSurface(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const NurbsSurfaceGeometry<TDimension, TPointType>& rNurbsSurface,
        const int MaxIterations = 30,
        const double Accuracy = 1e-6)
    {
        double lambda = 1e-8;
        const double lambda_factor = 10.0;
        const double max_lambda = 1e12;
        // Orthogonality condition
        const double orthogonality_tolerance = 1e-4;

        for (int i = 0; i < MaxIterations; ++i) {

            // Compute position and first derivatives only
            std::vector<array_1d<double, 3>> s;
            rNurbsSurface.GlobalSpaceDerivatives(s, rProjectedPointLocalCoordinates, 1);

            rProjectedPointGlobalCoordinates = s[0];

            const array_1d<double, 3> distance_vector = s[0] - rPointGlobalCoordinates;

            const double distance = norm_2(distance_vector);

            if (distance < Accuracy) {
                return true;
            }

            const double norm_su = norm_2(s[1]);
            const double norm_sv = norm_2(s[2]);

            if (norm_su < Accuracy || norm_sv < Accuracy) {
                return false;
            }

            // Gradient of 1/2 ||S(u,v)-x||^2
            const double g_u = inner_prod(s[1], distance_vector);
            const double g_v = inner_prod(s[2], distance_vector);

            const double xi_cos = std::abs(g_u) / (norm_su * distance);

            const double eta_cos = std::abs(g_v) / (norm_sv * distance);

            if (xi_cos < orthogonality_tolerance && eta_cos < orthogonality_tolerance) {
                return true;
            }

            // Gauss-Newton matrix J^T J with LM damping
            const double a_00 = inner_prod(s[1], s[1]);
            const double a_01 = inner_prod(s[1], s[2]);
            const double a_11 = inner_prod(s[2], s[2]);

            bool step_accepted = false;

            for (int lm_it = 0; lm_it < 10; ++lm_it) {

                const double j_00 = a_00 + lambda;
                const double j_01 = a_01;
                const double j_11 = a_11 + lambda;

                const double det_j = j_00 * j_11 - j_01 * j_01;

                if (std::abs(det_j) < Accuracy * (std::abs(j_00 * j_11) + std::abs(j_01 * j_01))) {
                    lambda *= lambda_factor;
                    continue;
                }

                // Solve:
                // (J^T J + lambda I) d = -g
                const double d_u =
                    (-g_u * j_11 + g_v * j_01) / det_j;

                const double d_v =
                    ( g_u * j_01 - g_v * j_00) / det_j;

                if (norm_2(d_u * s[1] + d_v * s[2]) < Accuracy) {
                    return xi_cos < orthogonality_tolerance && eta_cos < orthogonality_tolerance;
                }

                CoordinatesArrayType trial_local_coordinates =
                    rProjectedPointLocalCoordinates;

                trial_local_coordinates[0] += d_u;
                trial_local_coordinates[1] += d_v;

                rNurbsSurface.DomainIntervalU().IsInside(trial_local_coordinates[0]);
                rNurbsSurface.DomainIntervalV().IsInside(trial_local_coordinates[1]);

                std::vector<array_1d<double, 3>> trial_s;
                rNurbsSurface.GlobalSpaceDerivatives(trial_s, trial_local_coordinates, 0);

                const array_1d<double, 3> trial_distance_vector = trial_s[0] - rPointGlobalCoordinates;

                const double trial_distance = norm_2(trial_distance_vector);

                if (trial_distance < distance) {
                    rProjectedPointLocalCoordinates = trial_local_coordinates;
                    rProjectedPointGlobalCoordinates = trial_s[0];

                    lambda /= lambda_factor;
                    lambda = std::max(lambda, 1e-14);

                    step_accepted = true;
                    break;
                } else {
                    lambda *= lambda_factor;

                    if (lambda > max_lambda) {
                        return false;
                    }
                }
            }

            if (!step_accepted) {
                return false;
            }
        }

        return false;
    }
};
} // namespace Kratos

#endif // KRATOS_PROJECTION_NURBS_GEOMETRY_UTILITIES_H_INCLUDED
