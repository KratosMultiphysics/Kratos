//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Wataru Fukuda
//

#pragma once

// Project includes
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"

namespace Kratos {

/// Shape functions for a THB-spline surface (non-truncated: HB-spline variant).
///
/// At each evaluation point (u, v):
///   1. The finest active level is determined via ActiveLevelAtPoint.
///   2. Standard B-spline or NURBS shape functions are evaluated at that level.
///   3. The local CP indices are mapped to the packed Points() array of the
///      THBSurfaceGeometry via PackedControlPointIndex.
///
/// The interface is intentionally parallel to NurbsSurfaceShapeFunction so the
/// two classes are interchangeable inside CreateQuadraturePointGeometries.
///
/// @todo Implement ApplyTruncation to produce true THB shape functions.
///       Without truncation, partition of unity is not guaranteed at coarse/fine
///       boundaries (this is the HB variant, not THB).
class THBSurfaceShapeFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Static Helpers  (same layout as NurbsSurfaceShapeFunction)
    ///@{

    /// Number of shape function rows for a given derivative order.
    /// Row 0: N,  rows 1-2: dN/du dN/dv,  rows 3-5: second derivatives, ...
    static constexpr SizeType NumberOfShapeFunctionRows(SizeType DerivativeOrder) noexcept
    {
        return (1 + DerivativeOrder) * (2 + DerivativeOrder) / 2;
    }

    static constexpr IndexType IndexOfShapeFunctionRow(
        SizeType DerivativeOrderU, SizeType DerivativeOrderV) noexcept
    {
        return DerivativeOrderV +
               (DerivativeOrderU + DerivativeOrderV) *
               (1 + DerivativeOrderU + DerivativeOrderV) / 2;
    }

    ///@}
    ///@name Life Cycle
    ///@{

    THBSurfaceShapeFunction() = default;

    THBSurfaceShapeFunction(
        SizeType PolynomialDegreeU,
        SizeType PolynomialDegreeV,
        SizeType DerivativeOrder)
        : mDerivativeOrder(DerivativeOrder)
        , mNurbs(PolynomialDegreeU, PolynomialDegreeV, DerivativeOrder)
    {}

    ///@}
    ///@name Accessors
    ///@{

    SizeType DerivativeOrder() const { return mDerivativeOrder; }

    SizeType NumberOfNonzeroControlPoints() const
    {
        return static_cast<SizeType>(mControlPointIndices.size());
    }

    /// Shape function value: operator()(cp_index, derivative_row).
    /// cp_index runs over the nonzero CPs (0 .. NumberOfNonzeroControlPoints-1).
    /// derivative_row follows the NurbsSurfaceShapeFunction convention:
    ///   0 → N,  1 → dN/du,  2 → dN/dv,  3 → d²N/du², ...
    double operator()(IndexType ControlPointIndex, IndexType DerivativeRow) const
    {
        return mNurbs(ControlPointIndex, DerivativeRow);
    }

    /// Flat indices into the geometry's Points() array for each nonzero basis function.
    /// These already account for EliminateInactiveFunctions (packed layout).
    const std::vector<IndexType>& ControlPointIndices() const { return mControlPointIndices; }

    ///@}
    ///@name Evaluation
    ///@{

    /// Evaluate shape functions at parameter (u, v) for a given THBSurfaceGeometry.
    ///
    /// TTHBGeometry must expose:
    ///   - PolynomialDegree(direction)
    ///   - Levels()  →  container of THBLevel (with KnotsU, KnotsV, Weights)
    ///   - ActiveLevelAtPoint(u, v)
    ///   - PackedControlPointIndex(level, flat_index_within_level)
    template<class TTHBGeometry>
    void ComputeShapeFunctionValues(
        const TTHBGeometry& rGeometry,
        double u,
        double v)
    {
        const SizeType p = rGeometry.PolynomialDegree(0);
        const SizeType q = rGeometry.PolynomialDegree(1);
        const SizeType l = rGeometry.ActiveLevelAtPoint(u, v);
        const auto& lev = rGeometry.Levels()[l];

        mNurbs.ResizeDataContainers(p, q, mDerivativeOrder);

        if (lev.Weights.size() > 0) {
            mNurbs.ComputeNurbsShapeFunctionValues(lev.KnotsU, lev.KnotsV, lev.Weights, u, v);
        } else {
            mNurbs.ComputeBSplineShapeFunctionValues(lev.KnotsU, lev.KnotsV, u, v);
        }

        // Map local (within-level) flat CP indices to packed Points() indices.
        const SizeType nU_l = lev.KnotsU.size() - p + 1;
        const SizeType nV_l = lev.KnotsV.size() - q + 1;
        const auto local_flat = mNurbs.ControlPointIndices(nU_l, nV_l);

        mControlPointIndices.resize(local_flat.size());
        for (SizeType j = 0; j < local_flat.size(); ++j)
            mControlPointIndices[j] = rGeometry.PackedControlPointIndex(l, local_flat[j]);
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    SizeType                  mDerivativeOrder = 0;
    NurbsSurfaceShapeFunction mNurbs;
    std::vector<IndexType>    mControlPointIndices;

    ///@}
};

} // namespace Kratos
