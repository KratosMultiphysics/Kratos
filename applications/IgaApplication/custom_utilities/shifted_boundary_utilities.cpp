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
#include <algorithm>
#include <cmath>
#include <initializer_list>

// External includes

// Project includes
#include "custom_utilities/shifted_boundary_utilities.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_surface_geometry.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace
{

using IndexType = ShiftedBoundaryUtilities::IndexType;
using SizeType = ShiftedBoundaryUtilities::SizeType;
using GeometryType = ShiftedBoundaryUtilities::GeometryType;
using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;

/**
 * @brief Taylor coefficient of a multi-index
 *
 *     c_\boldsymbol\alpha = \mathbf{d}^{\boldsymbol\alpha} / \boldsymbol\alpha!
 *
 * @param rDistance Distance vector d from the facet point to the skin.
 * @param MultiIndex The multi-index \boldsymbol\alpha, one entry per direction.
 * @return The coefficient c_\boldsymbol\alpha.
 */
double ComputeTaylorCoefficient(
    const array_1d<double, 3>& rDistance,
    const std::initializer_list<IndexType> MultiIndex)
{
    double power = 1.0;
    IndexType factorial = 1;

    // \prod_j{ d_j^{\alpha_j} / \alpha_j! }
    for (IndexType j = 0; const IndexType alpha_j : MultiIndex) {
        power *= std::pow(rDistance[j], alpha_j);
        factorial *= MathUtils<double>::Factorial(alpha_j);
        ++j;
    }
    return power / factorial;
}

/**
 * @brief Taylor terms of a point's shape functions along a vector, summed over
 *        a band of orders, from the 1D B-spline bases of its knot span
 *
 *     Factor \sum_{\boldsymbol\alpha} c_\boldsymbol\alpha(x)
 *     \partial^{\boldsymbol\alpha} N,
 *     MinimumOrder \le |\boldsymbol\alpha| \le MaximumOrder,
 *     \alpha_u \le p_u, \alpha_v \le p_v
 *
 * @param rPoint A point of a B-spline surface.
 * @param rKnotSpanIndices Knot span of the point, one index per direction.
 * @param rVector Vector x of the Taylor coefficients.
 * @param MinimumOrder Smallest order |\boldsymbol\alpha| of the sum.
 * @param MaximumOrder Largest order |\boldsymbol\alpha| of the sum.
 * @param Factor Factor of the sum.
 * @return The sum for each shape function.
 */
Vector ComputeTaylorSum(
    const GeometryType& rPoint,
    const std::vector<IndexType>& rKnotSpanIndices,
    const array_1d<double, 3>& rVector,
    const SizeType MinimumOrder,
    const SizeType MaximumOrder,
    const double Factor)
{
    const auto& r_surface = dynamic_cast<const NurbsSurfaceType&>(rPoint.GetGeometryParent(0));
    const SizeType p_u = r_surface.PolynomialDegree(0);
    const SizeType p_v = r_surface.PolynomialDegree(1);

    // 1D bases of the span at the point, derivatives up to the degree
    const auto& r_point = rPoint.IntegrationPoints()[0];
    NurbsCurveShapeFunction shape_functions_u(p_u, p_u);
    shape_functions_u.ComputeBSplineShapeFunctionValuesAtSpan(
        r_surface.KnotsU(), rKnotSpanIndices[0] + p_u - 1, r_point.X()
    );
    NurbsCurveShapeFunction shape_functions_v(p_v, p_v);
    shape_functions_v.ComputeBSplineShapeFunctionValuesAtSpan(
        r_surface.KnotsV(), rKnotSpanIndices[1] + p_v - 1, r_point.Y()
    );

    // c_{(i, j)} N_a^{(i)} N_b^{(j)} for N_{a + b (p_u + 1)}, uv = xy
    Vector taylor_sum = ZeroVector(rPoint.size());
    for (IndexType i = 0; i <= std::min(p_u, MaximumOrder); ++i) {
        for (IndexType j = 0; j <= std::min(p_v, MaximumOrder - i); ++j) {
            if (i + j < MinimumOrder) {
                continue;  // below the orders of the sum
            }
            const double factor = Factor * ComputeTaylorCoefficient(rVector, {i, j});
            for (IndexType b = 0; b <= p_v; ++b) {
                for (IndexType a = 0; a <= p_u; ++a) {
                    taylor_sum[a + b * (p_u + 1)] +=
                        factor * shape_functions_u(a, i) * shape_functions_v(b, j);
                }
            }
        }
    }

    return taylor_sum;
}

} // namespace

Vector ShiftedBoundaryUtilities::ComputeDirectionalDerivatives(
    const GeometryType& rPoint,
    const std::vector<IndexType>& rKnotSpanIndices,
    const array_1d<double, 3>& rDirection,
    const SizeType Order)
{
    // D_s^n N = n! \sum_{|\alpha| = n} c_\alpha(s) \partial^\alpha N
    const double factorial = static_cast<double>(MathUtils<double>::Factorial(Order));
    return ComputeTaylorSum(rPoint, rKnotSpanIndices, rDirection, Order, Order, factorial);
}

Vector ShiftedBoundaryUtilities::ComputeShiftedShapeFunctions(
    const GeometryType& rPoint,
    const std::vector<IndexType>& rKnotSpanIndices,
    const array_1d<double, 3>& rDistance,
    const SizeType Order)
{
    // S_d^m[N] = \sum_{|\alpha| \le m} c_\alpha(d) \partial^\alpha N
    return ComputeTaylorSum(rPoint, rKnotSpanIndices, rDistance, 0, Order, 1.0);
}

} // namespace Kratos
