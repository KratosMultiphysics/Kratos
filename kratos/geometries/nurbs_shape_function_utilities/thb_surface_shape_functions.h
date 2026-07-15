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
///   1. The finest active level L is determined via ActiveLevelAtPoint.
///   2. B-spline shape functions are evaluated at EVERY level l = 0 .. L.
///   3. For each level l, only ACTIVE basis functions are kept.
///      A function is active if its support is not fully covered by the next
///      refinement domain (ComputeActiveFunctions convention).
///   4. Each active function's truncated value is computed:
///        trunc(B_i^l)(u,v) = B_i^l(u,v)
///                           − sum_{j active at l+1} c_{ij} * B_j^{l+1}(u,v)
///      where c_{ij} are the refinement coefficients (GetTruncationData).
///   5. CP indices are mapped to the packed Points() array of the geometry.
///
/// This multi-level approach correctly handles evaluation points near refinement
/// domain boundaries, where some level-L B-splines may be inactive (their parent
/// straddles the boundary and is still represented by a coarser active function).
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
        , mTensorProductNurbs(PolynomialDegreeU, PolynomialDegreeV, DerivativeOrder)
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

    template<class TTHBGeometry>
    void ComputeShapeFunctionValues(
        const TTHBGeometry& rGeometry,
        double u,
        double v)
    {
        const SizeType ActiveLevel = rGeometry.ActiveLevelAtPoint(u, v);
        const auto& AllLevels = rGeometry.Levels();
        const SizeType number_of_shape_function_rows = NumberOfShapeFunctionRows(mDerivativeOrder);

        // evaluate B-splines at every level 0..ActiveLevel
        // Store raw values and build a flat-index to local-position lookup per level.
        struct LevelCache {
            SizeType number_of_nonzero_control_points = 0;
            std::vector<int> local_flat_indices;
            std::vector<double> values;
            std::unordered_map<SizeType, SizeType> flat_index_to_local_position;
        };

        std::vector<LevelCache> level_caches(ActiveLevel + 1);

        for (SizeType l = 0; l <= ActiveLevel; ++l) {
            const auto& CurrentLevel = AllLevels[l];
            const SizeType degree_u = CurrentLevel.DegreeU;
            const SizeType degree_v = CurrentLevel.DegreeV;
            mTensorProductNurbs.ResizeDataContainers(degree_u, degree_v, mDerivativeOrder);
            if (!CurrentLevel.Weights.empty())
                mTensorProductNurbs.ComputeNurbsShapeFunctionValues(
                    CurrentLevel.KnotsU, CurrentLevel.KnotsV, CurrentLevel.Weights, u, v);
            else
                mTensorProductNurbs.ComputeBSplineShapeFunctionValues(
                    CurrentLevel.KnotsU, CurrentLevel.KnotsV, u, v);

            auto& CurrentLevelCache = level_caches[l];
            CurrentLevelCache.number_of_nonzero_control_points = mTensorProductNurbs.NumberOfNonzeroControlPoints();
            const SizeType NumberOfControlPointsU = CurrentLevel.KnotsU.size() - degree_u + 1;
            const SizeType NumberOfControlPointsV = CurrentLevel.KnotsV.size() - degree_v + 1;
            CurrentLevelCache.local_flat_indices = mTensorProductNurbs.ControlPointIndices(
                                                        NumberOfControlPointsU, NumberOfControlPointsV);
            CurrentLevelCache.values.resize(
                number_of_shape_function_rows * CurrentLevelCache.number_of_nonzero_control_points);
            for (SizeType row = 0; row < number_of_shape_function_rows; ++row)
                for (SizeType k = 0; k < CurrentLevelCache.number_of_nonzero_control_points; ++k)
                    CurrentLevelCache.values[row * CurrentLevelCache.number_of_nonzero_control_points + k]
                        = mTensorProductNurbs(k, row);
            for (SizeType k = 0; k < CurrentLevelCache.number_of_nonzero_control_points; ++k)
                CurrentLevelCache.flat_index_to_local_position[static_cast<SizeType>(CurrentLevelCache.local_flat_indices[k])] = k;
        }

        // Collect active contributions with multi-level truncation.
        // τ(B_i^l) = B_i^l − Σ_{m>l} Σ_{j active at m} D_{ij}^m * B_j^m
        // where D coefficients are precomputed in ComputeTruncationData.
        mControlPointIndices.clear();
        std::vector<std::vector<double>> truncated_value_per_control_point;

        for (SizeType l = 0; l <= ActiveLevel; ++l) {
            const auto& CurrentLevelCache = level_caches[l];
            const auto& active_flags = rGeometry.GetActiveFunctions(l);

            for (SizeType j = 0; j < CurrentLevelCache.number_of_nonzero_control_points; ++j) {
                const SizeType flat_index = static_cast<SizeType>(CurrentLevelCache.local_flat_indices[j]);
                if (flat_index >= active_flags.size() || !active_flags[flat_index]) continue;

                std::vector<double> truncated_value(number_of_shape_function_rows);
                for (SizeType row = 0; row < number_of_shape_function_rows; ++row)
                    truncated_value[row] = CurrentLevelCache.values[row * CurrentLevelCache.number_of_nonzero_control_points + j];

                // Multi-level truncation: subtract active fine-level contributions.
                // Entries can target any FineLevel from l+1 to num_levels-1;
                // skip levels beyond the active range (those B-splines are zero here).
                for (const auto& TruncationEntry : rGeometry.GetTruncationData(l, flat_index)) {
                    if (TruncationEntry.FineLevel > ActiveLevel) continue;
                    const auto& FinerLevelCache = level_caches[TruncationEntry.FineLevel];
                    auto position_iterator = FinerLevelCache.flat_index_to_local_position.find(TruncationEntry.FineFlatIndex);
                    if (position_iterator == FinerLevelCache.flat_index_to_local_position.end()) continue;
                    const SizeType local_position = position_iterator->second;
                    for (SizeType row = 0; row < number_of_shape_function_rows; ++row)
                        truncated_value[row] -= TruncationEntry.Coefficient
                            * FinerLevelCache.values[row * FinerLevelCache.number_of_nonzero_control_points + local_position];
                }

                mControlPointIndices.push_back(rGeometry.PackedControlPointIndex(l, flat_index));
                truncated_value_per_control_point.push_back(std::move(truncated_value));
            }
        }

        // pack into [row * total + cp_idx] layout
        const SizeType total = mControlPointIndices.size();
        mValues.assign(number_of_shape_function_rows * total, 0.0);
        for (SizeType j = 0; j < total; ++j)
            for (SizeType row = 0; row < number_of_shape_function_rows; ++row)
                mValues[row * total + j] = truncated_value_per_control_point[j][row];
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    SizeType                  mDerivativeOrder = 0;
    NurbsSurfaceShapeFunction mTensorProductNurbs;
    std::vector<IndexType>    mControlPointIndices;
    /// Shape function values: layout [row * num_nonzero + cp_idx].
    std::vector<double>       mValues;

    ///@}
};

} // namespace Kratos
