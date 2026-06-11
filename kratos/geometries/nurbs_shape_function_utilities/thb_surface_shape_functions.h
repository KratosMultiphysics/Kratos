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

// System includes
#include <unordered_map>

// Project includes
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"

namespace Kratos {

/// Shape functions for a THB-spline surface.
///
/// At each evaluation point (u, v):
///   1. The finest active level l is determined via ActiveLevelAtPoint.
///   2. B-spline / NURBS shape functions at level l are evaluated.
///   3. Truncation is applied: for each nonzero coarse CP that has active
///      children at level l+1, the fine-level B-spline values at (u,v) are
///      evaluated and their refinement coefficients are subtracted.
///      This produces true THB shape functions that satisfy partition of unity.
///   4. CP indices are mapped to the packed Points() array of the geometry.
///
/// @note Truncation is single-level (level l → l+1 only).  For 3+ levels the
///       recursive case (subtracting trunc(B^{l+1}) rather than B^{l+1}) is
///       handled by iterating the correction from l up to L-2. @todo verify.
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
        return mValues[DerivativeRow * NumberOfNonzeroControlPoints() + ControlPointIndex];
    }

    /// Flat indices into the geometry's Points() array for each nonzero basis function.
    /// These already account for EliminateInactiveFunctions (packed layout).
    const std::vector<IndexType>& ControlPointIndices() const { return mControlPointIndices; }

    ///@}
    ///@name Evaluation
    ///@{

    /// Evaluate (truncated) shape functions at parameter (u, v).
    ///
    /// TTHBGeometry must expose:
    ///   - PolynomialDegree(direction)
    ///   - Levels()              → container of THBLevel (KnotsU, KnotsV, Weights)
    ///   - ActiveLevelAtPoint(u, v)
    ///   - PackedControlPointIndex(level, flat_index_within_level)
    ///   - GetTruncationData(level, flat_index)
    ///                           → const vector of TruncationEntry {FineFlatIndex, Coefficient}
    template<class TTHBGeometry>
    void ComputeShapeFunctionValues(
        const TTHBGeometry& rGeometry,
        double u,
        double v)
    {
        const SizeType p = rGeometry.PolynomialDegree(0);
        const SizeType q = rGeometry.PolynomialDegree(1);
        const SizeType l = rGeometry.ActiveLevelAtPoint(u, v);
        const auto& levels = rGeometry.Levels();
        const auto& lev = levels[l];

        mNurbs.ResizeDataContainers(p, q, mDerivativeOrder);
        if (lev.Weights.size() > 0)
            mNurbs.ComputeNurbsShapeFunctionValues(lev.KnotsU, lev.KnotsV, lev.Weights, u, v);
        else
            mNurbs.ComputeBSplineShapeFunctionValues(lev.KnotsU, lev.KnotsV, u, v);

        const SizeType num_nonzero = mNurbs.NumberOfNonzeroControlPoints();
        const SizeType nU_l = lev.KnotsU.size() - p + 1;
        const SizeType nV_l = lev.KnotsV.size() - q + 1;
        const auto local_flat = mNurbs.ControlPointIndices(nU_l, nV_l);

        // Copy B-spline values into mValues so they can be modified by truncation.
        const SizeType n_rows = NumberOfShapeFunctionRows(mDerivativeOrder);
        mValues.resize(n_rows * num_nonzero);
        for (SizeType row = 0; row < n_rows; ++row)
            for (SizeType j = 0; j < num_nonzero; ++j)
                mValues[row * num_nonzero + j] = mNurbs(j, row);

        // Map local flat indices to packed Points() indices.
        mControlPointIndices.resize(num_nonzero);
        for (SizeType j = 0; j < num_nonzero; ++j)
            mControlPointIndices[j] = rGeometry.PackedControlPointIndex(l, local_flat[j]);

        ApplyTruncation(rGeometry, u, v, l, local_flat, num_nonzero, p, q);
    }

    ///@}

private:
    ///@name Private Methods
    ///@{

    /// Subtracts fine-level B-spline contributions from coarse CPs that straddle
    /// the refinement boundary (those with non-empty truncation data).
    ///
    /// For each nonzero coarse CP i at level l:
    ///   trunc(N_i^l)(u,v) = N_i^l(u,v) - sum_{j active at l+1} c_{ij} * N_j^{l+1}(u,v)
    ///
    /// This is iterated from l up to L-2 to handle multi-level hierarchies.
    template<class TTHBGeometry, class TLocalFlat>
    void ApplyTruncation(
        const TTHBGeometry& rGeometry,
        double u, double v,
        SizeType l,
        const TLocalFlat& rLocalFlat,
        SizeType NumNonzero,
        SizeType p, SizeType q)
    {
        const auto& levels = rGeometry.Levels();
        if (l + 1 >= levels.size()) return;

        // Check whether any nonzero CP at level l has truncation entries.
        bool needs_truncation = false;
        for (SizeType j = 0; j < NumNonzero && !needs_truncation; ++j)
            if (!rGeometry.GetTruncationData(l, rLocalFlat[j]).empty())
                needs_truncation = true;
        if (!needs_truncation) return;

        // Evaluate fine (level l+1) B-splines at (u, v).
        // These can be nonzero at (u,v) even when (u,v) is outside Omega^{l+1}
        // because the fine functions near the boundary straddle it.
        const auto& lev_f = levels[l + 1];
        NurbsSurfaceShapeFunction fine_sf(p, q, mDerivativeOrder);
        if (lev_f.Weights.size() > 0)
            fine_sf.ComputeNurbsShapeFunctionValues(lev_f.KnotsU, lev_f.KnotsV, lev_f.Weights, u, v);
        else
            fine_sf.ComputeBSplineShapeFunctionValues(lev_f.KnotsU, lev_f.KnotsV, u, v);

        const SizeType nU_f = lev_f.KnotsU.size() - p + 1;
        const SizeType nV_f = lev_f.KnotsV.size() - q + 1;
        const auto fine_local = fine_sf.ControlPointIndices(nU_f, nV_f);

        // Build a flat-index → local-index lookup for the nonzero fine functions.
        std::unordered_map<SizeType, SizeType> fine_lookup;
        fine_lookup.reserve(fine_local.size());
        for (SizeType k = 0; k < fine_local.size(); ++k)
            fine_lookup[static_cast<SizeType>(fine_local[k])] = k;

        // Subtract: N_i^l -= sum_j c_{ij} * N_j^{l+1}
        const SizeType n_rows = NumberOfShapeFunctionRows(mDerivativeOrder);
        for (SizeType j = 0; j < NumNonzero; ++j) {
            const auto& entries = rGeometry.GetTruncationData(l, rLocalFlat[j]);
            for (const auto& e : entries) {
                auto it = fine_lookup.find(e.FineFlatIndex);
                if (it == fine_lookup.end()) continue;
                const SizeType k = it->second;
                for (SizeType row = 0; row < n_rows; ++row)
                    mValues[row * NumNonzero + j] -= e.Coefficient * fine_sf(k, row);
            }
        }
    }

    ///@}
    ///@name Member Variables
    ///@{

    SizeType                  mDerivativeOrder = 0;
    NurbsSurfaceShapeFunction mNurbs;
    std::vector<IndexType>    mControlPointIndices;
    /// Shape function values after truncation: layout [row * num_nonzero + cp_idx].
    std::vector<double>       mValues;

    ///@}
};

} // namespace Kratos
