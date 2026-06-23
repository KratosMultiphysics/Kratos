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
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"

namespace Kratos {

/// Shape functions for a THB-spline curve.
///
/// At each evaluation point t:
///   1. The finest active level L is determined via ActiveLevelAtPoint.
///   2. B-spline shape functions are evaluated at EVERY level l = 0 .. L.
///   3. For each level l, only ACTIVE basis functions are kept.
///   4. Each active function's truncated value is computed:
///        trunc(B_i^l)(t) = B_i^l(t) − sum_{j active at l+1} c_{ij} * B_j^{l+1}(t)
///      where c_{ij} are the 1-D refinement coefficients (GetTruncationData).
///   5. CP indices are mapped to the packed Points() array of the geometry.
///
/// mValues layout: [row * num_nonzero + cp_idx]
///   row 0 → N,   row 1 → dN/dt,   row 2 → d²N/dt², ...
class THBCurveShapeFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Static Helpers
    ///@{

    /// Number of shape function rows: N, dN/dt, d²N/dt², ...
    static constexpr SizeType NumberOfShapeFunctionRows(SizeType DerivativeOrder) noexcept
    {
        return DerivativeOrder + 1;
    }

    ///@}
    ///@name Life Cycle
    ///@{

    THBCurveShapeFunction() = default;

    THBCurveShapeFunction(SizeType PolynomialDegree, SizeType DerivativeOrder)
        : mDerivativeOrder(DerivativeOrder)
        , mTensorProductNurbs(PolynomialDegree, DerivativeOrder)
    {}

    ///@}
    ///@name Accessors
    ///@{

    SizeType DerivativeOrder() const { return mDerivativeOrder; }

    SizeType NumberOfNonzeroControlPoints() const
    {
        return static_cast<SizeType>(mControlPointIndices.size());
    }

    /// sf(cp_index, row): row 0 = N,  row 1 = dN/dt,  row 2 = d²N/dt², ...
    double operator()(IndexType ControlPointIndex, IndexType DerivativeRow) const
    {
        return mValues[DerivativeRow * NumberOfNonzeroControlPoints() + ControlPointIndex];
    }

    /// Packed indices into the geometry's Points() array for each nonzero basis function.
    const std::vector<IndexType>& ControlPointIndices() const { return mControlPointIndices; }

    ///@}
    ///@name Evaluation
    ///@{

    template<class TTHBGeometry>
    void ComputeShapeFunctionValues(const TTHBGeometry& rGeometry, double t)
    {
        const SizeType PolynomialDegreeT = rGeometry.PolynomialDegree(0);
        const SizeType ActiveLevel = rGeometry.ActiveLevelAtPoint(t);
        const auto& AllLevels = rGeometry.Levels();
        const SizeType number_of_shape_function_rows = NumberOfShapeFunctionRows(mDerivativeOrder);

        // evaluate B-splines at every level 0..ActiveLevel.
        // For a curve, the (PolynomialDegreeT+1) nonzero CPs at parameter t start at GetFirstNonzeroControlPoint().
        struct LevelCache {
            SizeType first_nonzero_cp = 0;
            SizeType number_of_nonzero_control_points = 0;
            std::vector<double> values;
            std::unordered_map<SizeType, SizeType> flat_index_to_local_position;
        };

        std::vector<LevelCache> level_caches(ActiveLevel + 1);

        for (SizeType l = 0; l <= ActiveLevel; ++l) {
            const auto& CurrentLevel = AllLevels[l];
            mTensorProductNurbs.ResizeDataContainers(PolynomialDegreeT, mDerivativeOrder);
            if (!CurrentLevel.Weights.empty())
                mTensorProductNurbs.ComputeNurbsShapeFunctionValues(CurrentLevel.Knots, CurrentLevel.Weights, t);
            else
                mTensorProductNurbs.ComputeBSplineShapeFunctionValues(CurrentLevel.Knots, t);

            auto& CurrentLevelCache = level_caches[l];
            CurrentLevelCache.first_nonzero_cp = mTensorProductNurbs.GetFirstNonzeroControlPoint();
            CurrentLevelCache.number_of_nonzero_control_points = mTensorProductNurbs.NumberOfNonzeroControlPoints();

            CurrentLevelCache.values.resize(
                number_of_shape_function_rows * CurrentLevelCache.number_of_nonzero_control_points);
            for (SizeType row = 0; row < number_of_shape_function_rows; ++row)
                for (SizeType k = 0; k < CurrentLevelCache.number_of_nonzero_control_points; ++k)
                    CurrentLevelCache.values[row * CurrentLevelCache.number_of_nonzero_control_points + k]
                        = mTensorProductNurbs(k, row);

            for (SizeType k = 0; k < CurrentLevelCache.number_of_nonzero_control_points; ++k)
                CurrentLevelCache.flat_index_to_local_position[CurrentLevelCache.first_nonzero_cp + k] = k;
        }

        // Collect active contributions with multi-level truncation.
        // τ(B_i^l) = B_i^l − Σ_{m>l} Σ_{j active at m} D_{ij}^m * B_j^m
        // where D coefficients are precomputed in ComputeTruncationData.
        mControlPointIndices.clear();
        std::vector<std::vector<double>> truncated_value_per_control_point;

        for (SizeType l = 0; l <= ActiveLevel; ++l) {
            const auto& CurrentLevelCache = level_caches[l];
            const auto& active_flags = rGeometry.GetActiveFunctions(l);

            for (SizeType k = 0; k < CurrentLevelCache.number_of_nonzero_control_points; ++k) {
                const SizeType flat_index = CurrentLevelCache.first_nonzero_cp + k;
                if (flat_index >= active_flags.size() || !active_flags[flat_index]) continue;

                std::vector<double> truncated_value(number_of_shape_function_rows);
                for (SizeType row = 0; row < number_of_shape_function_rows; ++row)
                    truncated_value[row] = CurrentLevelCache.values[row * CurrentLevelCache.number_of_nonzero_control_points + k];

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

        // pack into [row * total + cp_idx] layout.
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

    SizeType                mDerivativeOrder = 0;
    NurbsCurveShapeFunction mTensorProductNurbs;
    std::vector<IndexType>  mControlPointIndices;
    /// Shape function values: layout [row * num_nonzero + cp_idx].
    std::vector<double>     mValues;

    ///@}
};

} // namespace Kratos
