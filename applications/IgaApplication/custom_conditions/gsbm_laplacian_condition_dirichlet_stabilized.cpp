//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

// System includes
#include <vector>

// External includes

// Project includes
#include "custom_conditions/gsbm_laplacian_condition_dirichlet_stabilized.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"
#include "geometries/nurbs_surface_geometry.h"
#include "iga_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace
{

using IndexType = GsbmLaplacianConditionDirichletStabilized::IndexType;
using SizeType = GsbmLaplacianConditionDirichletStabilized::SizeType;
using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;
using Directions = std::vector<array_1d<double, 3>>;

/**
 * @brief Coefficients e_k of
 * \nabla^n u[v_1, ..., v_n] = \sum_k e_k \partial_x^(n-k) \partial_y^k u.
 *
 * @param rDirections The directions v_1, ..., v_n.
 * @return [e_0, ..., e_n].
 */
Vector ComputeContractionCoefficients(const Directions& rDirections)
{
    Vector coefficients(1, 1.0);
    for (const auto& r_direction : rDirections) {
        Vector product = ZeroVector(coefficients.size() + 1);
        for (IndexType k = 0; k < coefficients.size(); ++k) {
            product[k] += r_direction[0] * coefficients[k];
            product[k + 1] += r_direction[1] * coefficients[k];
        }
        coefficients = product;
    }
    return coefficients;
}

/**
 * @brief Coefficients of t^k D_k, with
 * D_k = \sum_i (-1)^i \partial_\perp^(2i) \partial_e^(k-2-2i).
 *
 * @param rDistance Distance vector d from the facet point to the skin.
 * @param K The order k >= 2.
 * @return |d|^2 \sum_i (-1)^i \nabla^(k-2)[(d_\perp)^(2i), (-d)^(k-2-2i)].
 */
Vector ComputeResidualOperatorCoefficients(
    const array_1d<double, 3>& rDistance,
    const int K)
{
    // t e = -d, t e_\perp = d_\perp
    const array_1d<double, 3> minus_distance = -rDistance;
    array_1d<double, 3> distance_perp = ZeroVector(3);
    distance_perp[0] = -rDistance[1];
    distance_perp[1] = rDistance[0];

    Vector coefficients = ZeroVector(K - 1);
    for (int i = 0; 2 * i <= K - 2; ++i) {
        Directions directions(2 * i, distance_perp);
        directions.insert(directions.end(), K - 2 - 2 * i, minus_distance);
        const double sign = i % 2 == 0 ? 1.0 : -1.0;
        noalias(coefficients) += sign * ComputeContractionCoefficients(directions);
    }
    return inner_prod(rDistance, rDistance) * coefficients;
}

/// \rho^m is added unless USE_RESIDUAL_TERMS is false
bool UsesResidualTerms(const Properties& rProperties)
{
    return rProperties.Has(USE_RESIDUAL_TERMS) ? rProperties[USE_RESIDUAL_TERMS] : true;
}

} // namespace

Vector GsbmLaplacianConditionDirichletStabilized::CalculateShiftedShapeFunctions(
    const GeometryType& rSide,
    const Node& rProjectionNode,
    const array_1d<double, 3>& rDistance,
    const SizeType Order) const
{
    const int m = static_cast<int>(Order);
    const array_1d<double, 3> minus_distance = -rDistance;  // t e = -d

    // Weights w of \partial_x^(n-k) \partial_y^k N(x_\Gamma)
    std::vector<Vector> weights(m + 1);
    for (int k = 1; k <= m; ++k) {
        // Series from x_\Gamma: -(t^k / k!) \partial_e^k N
        weights[k] = -ComputeContractionCoefficients(Directions(k, minus_distance))
            / MathUtils<double>::Factorial(k);
    }
    if (UsesResidualTerms(GetProperties())) {
        for (int k = 2; k <= m; ++k) {
            // \rho^m: (t^k / k!) D_k \Delta N
            const Vector coefficients =
                ComputeResidualOperatorCoefficients(rDistance, k) / MathUtils<double>::Factorial(k);
            for (IndexType l = 0; l < coefficients.size(); ++l) {
                weights[k][l] += coefficients[l];      // \partial_x^2
                weights[k][l + 2] += coefficients[l];  // \partial_y^2
            }
        }
    }

    // \partial^\alpha N(x_\Gamma) of the side's span polynomial (uv = xy)
    const auto p_surface = dynamic_cast<const NurbsSurfaceType*>(&rSide.GetGeometryParent(0));
    KRATOS_ERROR_IF(p_surface == nullptr) << Info() << ": the parent of a facet side is "
        << "not a 2D NURBS surface." << std::endl;
    const SizeType p_u = p_surface->PolynomialDegree(0);
    const SizeType p_v = p_surface->PolynomialDegree(1);
    const auto& r_spans = rSide.GetValue(KNOT_SPAN_INDICES);
    NurbsSurfaceShapeFunction shape_functions(p_u, p_v, Order);
    shape_functions.ComputeBSplineShapeFunctionValuesAtSpan(
        p_surface->KnotsU(), p_surface->KnotsV(),
        r_spans[0] + p_u - 1, r_spans[1] + p_v - 1,
        rProjectionNode.X(), rProjectionNode.Y()
    );

    // Shifted N = N(\tilde{x}) + \sum w_\alpha \partial^\alpha N(x_\Gamma)
    Vector shifted_shape_functions = row(rSide.ShapeFunctionsValues(), 0);
    for (int n = 1; n <= m; ++n) {
        for (int k = 0; k <= n; ++k) {
            const IndexType derivative_row =
                NurbsSurfaceShapeFunction::IndexOfShapeFunctionRow(n - k, k);
            for (IndexType i = 0; i < shifted_shape_functions.size(); ++i) {
                shifted_shape_functions[i] += weights[n][k] * shape_functions(i, derivative_row);
            }
        }
    }

    return shifted_shape_functions;
}

double GsbmLaplacianConditionDirichletStabilized::CalculateShiftedBoundaryValue(
    const Node& rProjectionNode,
    const array_1d<double, 3>& rDistance,
    const SizeType Order,
    const Variable<double>& rUnknownVariable) const
{
    const double boundary_value = rProjectionNode.GetValue(rUnknownVariable);
    const int m = static_cast<int>(Order);
    if (m < 2 || !UsesResidualTerms(GetProperties())) {
        return boundary_value;  // G = g
    }

    // \partial_x^(n-k) \partial_y^k f(x_\Gamma), n = 0 to m - 2
    const Vector& r_f = rProjectionNode.GetValue(VOLUME_SOURCE_DERIVATIVES);
    const SizeType f_size = static_cast<SizeType>((m - 1) * m / 2);
    KRATOS_ERROR_IF(r_f.size() < f_size) << Info() << ": projection node #"
        << rProjectionNode.Id() << " needs VOLUME_SOURCE_DERIVATIVES, f and its Cartesian "
        << "derivatives of orders 0 to " << m - 2 << " (size " << f_size << "), but it "
        << "has size " << r_f.size() << "." << std::endl;

    // G = g - \sum_k (t^k / k!) D_k f
    double value = boundary_value;
    for (int k = 2; k <= m; ++k) {
        const Vector coefficients =
            ComputeResidualOperatorCoefficients(rDistance, k) / MathUtils<double>::Factorial(k);
        const IndexType offset = static_cast<IndexType>((k - 2) * (k - 1) / 2);  // n = k - 2
        for (IndexType l = 0; l < coefficients.size(); ++l) {
            value -= coefficients[l] * r_f[offset + l];
        }
    }

    return value;
}

} // namespace Kratos
