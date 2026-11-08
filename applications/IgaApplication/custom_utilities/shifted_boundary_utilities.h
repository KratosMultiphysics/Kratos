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

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) ShiftedBoundaryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;
    using SizeType = std::size_t;
    using GeometryType = Geometry<Node>;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Derivatives of a point's shape functions along a unit vector,
     *        from the 1D B-spline bases of its knot span
     *
     *     D_s^n N = n! \sum_{|\boldsymbol\alpha| = n} c_\boldsymbol\alpha(s)
     *               \partial^{\boldsymbol\alpha} N,
     *               \alpha_u \le p_u, \alpha_v \le p_v
     *
     * with c_\boldsymbol\alpha(s) = s^{\boldsymbol\alpha} / \boldsymbol\alpha!,
     * the terms of order n of the Taylor series along s.
     *
     * @param rPoint A point of a B-spline surface.
     * @param rKnotSpanIndices Knot span of the point, one index per direction.
     * @param rDirection Unit vector s of the derivatives.
     * @param Order Order n of the derivatives.
     * @return The directional derivatives of the shape functions.
     */
    static Vector ComputeDirectionalDerivatives(
        const GeometryType& rPoint,
        const std::vector<IndexType>& rKnotSpanIndices,
        const array_1d<double, 3>& rDirection,
        const SizeType Order);

    /**
     * @brief Derivatives of a partial derivative of a point's shape functions
     *        along a unit vector, from the 1D B-spline bases of its knot span
     *
     *     D_s^n \partial_j N = n! \sum_{|\boldsymbol\alpha| = n}
     *                         c_\boldsymbol\alpha(s)
     *                         \partial^{\boldsymbol\alpha + e_j} N
     *
     * @param rPoint A point of a B-spline surface.
     * @param rKnotSpanIndices Knot span of the point, one index per direction.
     * @param rDirection Unit vector s of the derivatives.
     * @param Order Order n of the derivatives.
     * @param PartialDirection Direction j of the partial derivative.
     * @return The directional derivatives of the partial derivatives.
     */
    static Vector ComputeDirectionalDerivatives(
        const GeometryType& rPoint,
        const std::vector<IndexType>& rKnotSpanIndices,
        const array_1d<double, 3>& rDirection,
        const SizeType Order,
        const IndexType PartialDirection);

    /**
     * @brief Shape functions of a point shifted to the true boundary, from the
     *        1D B-spline bases of its knot span
     *
     *     S_d^m[N] = \sum_{|\boldsymbol\alpha| \le m} c_\boldsymbol\alpha(d)
     *                \partial^{\boldsymbol\alpha} N(\tilde{x}),
     *                \alpha_u \le p_u, \alpha_v \le p_v
     *
     * m = p_u + p_v is the full Taylor series, the span's polynomial at
     * \tilde{x} + d. The 1D derivatives never exceed the degree, so every order
     * m is available.
     *
     * @param rPoint A point of a B-spline surface, on a facet side.
     * @param rKnotSpanIndices Knot span of the point, one index per direction.
     * @param rDistance Distance vector d from the facet point to the skin.
     * @param Order Order m of the shift.
     * @return The shifted shape functions.
     */
    static Vector ComputeShiftedShapeFunctions(
        const GeometryType& rPoint,
        const std::vector<IndexType>& rKnotSpanIndices,
        const array_1d<double, 3>& rDistance,
        const SizeType Order);

    /**
     * @brief Gradients of a point's shape functions shifted to the true
     *        boundary, from the 1D B-spline bases of its knot span
     *
     *     S_d^m[\partial_j N] = \sum_{|\boldsymbol\alpha| \le m}
     *                           c_\boldsymbol\alpha(d)
     *                           \partial^{\boldsymbol\alpha + e_j} N(\tilde{x})
     *
     * @param rPoint A point of a B-spline surface, on a facet side.
     * @param rKnotSpanIndices Knot span of the point, one index per direction.
     * @param rDistance Distance vector d from the facet point to the skin.
     * @param Order Order m of the shift.
     * @return The shifted gradients, one row per shape function and one
     *         column per direction, as the shape function local gradients.
     */
    static Matrix ComputeShiftedGradients(
        const GeometryType& rPoint,
        const std::vector<IndexType>& rKnotSpanIndices,
        const array_1d<double, 3>& rDistance,
        const SizeType Order);

    ///@}
};

///@}

} // namespace Kratos
