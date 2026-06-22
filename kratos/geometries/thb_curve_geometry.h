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
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "geometries/local_refined_curve_geometry.h"
#include "geometries/nurbs_shape_function_utilities/thb_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"
#include "utilities/quadrature_points_utility.h"
#include "integration/integration_point_utilities.h"

namespace Kratos {

/// Truncated Hierarchical B-spline (THB) curve geometry.
///
/// A single parametric direction t.  Levels are built by uniform knot bisection
/// (AddLevel) and refined regions are registered with AddRefinementDomain.
/// After setup, EliminateInactiveFunctions compacts the control-point array and
/// removes inactive nodes from the ModelPart.
///
/// @tparam TWorkingSpaceDimension  Embedding space dimension (2 or 3).
/// @tparam TContainerPointType     Container of control-point nodes.
template <int TWorkingSpaceDimension, class TContainerPointType>
class THBCurveGeometry
    : public LocalRefinedCurveGeometry<TWorkingSpaceDimension, TContainerPointType>
{
public:
    ///@name Type Definitions
    ///@{

    typedef LocalRefinedCurveGeometry<TWorkingSpaceDimension, TContainerPointType> BaseType;

    typedef typename BaseType::NodeType                   NodeType;
    typedef typename BaseType::IndexType                  IndexType;
    typedef typename BaseType::SizeType                   SizeType;
    typedef typename BaseType::CoordinatesArrayType       CoordinatesArrayType;
    typedef typename BaseType::PointsArrayType            PointsArrayType;
    typedef typename BaseType::GeometriesArrayType        GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    struct THBLevel {
        SizeType Degree;
        Vector   Knots;
        Vector   Weights;
    };

    struct THBRefinementDomain {
        SizeType Level;
        double MinT;
        double MaxT;
    };

    struct TruncationEntry {
        SizeType FineLevel;
        SizeType FineFlatIndex;
        double   Coefficient;
    };

    using THBLevelContainer         = std::vector<THBLevel>;
    using RefinementDomainContainer = std::vector<THBRefinementDomain>;
    using TruncationEntryContainer  = std::vector<TruncationEntry>;

    KRATOS_CLASS_POINTER_DEFINITION(THBCurveGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructs the initial level-0 B-Spline curve.
    THBCurveGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegree,
        const Vector& rKnots)
        : BaseType(rThisPoints)
    {
        mLevels.push_back(THBLevel{PolynomialDegree, rKnots, {}});
        const SizeType num_cps = rKnots.size() - PolynomialDegree + 1;
        mActiveFunctions.push_back(std::vector<bool>(num_cps, true));
    }

    /// Constructs the initial level-0 NURBS curve (with weights).
    THBCurveGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegree,
        const Vector& rKnots,
        const Vector& rWeights)
        : BaseType(rThisPoints)
    {
        KRATOS_ERROR_IF(rWeights.size() != rThisPoints.size())
            << "Number of control points and weights do not match!" << std::endl;

        mLevels.push_back(THBLevel{PolynomialDegree, rKnots, rWeights});
        const SizeType num_cps = rKnots.size() - PolynomialDegree + 1;
        mActiveFunctions.push_back(std::vector<bool>(num_cps, true));
    }

    ~THBCurveGeometry() override = default;

    ///@}
    ///@name Refinement
    ///@{

    /// Adds one level with an explicit knot vector. CPs computed via knot insertion;
    /// new nodes are orphan (not registered in any ModelPart).
    void AddLevel(const Vector& rKnots, const Vector& rWeights = Vector())
    {
        AddLevelImpl(rKnots, rWeights,
            [](IndexType id, double x, double y, double z) {
                return Kratos::make_intrusive<NodeType>(id, x, y, z);
            });
    }

    /// Adds one level with an explicit knot vector. CPs are registered in rModelPart.
    void AddLevel(const Vector& rKnots, ModelPart& rModelPart, const Vector& rWeights = Vector())
    {
        AddLevelImpl(rKnots, rWeights,
            [&rModelPart](IndexType id, double x, double y, double z)
                -> typename NodeType::Pointer {
                return rModelPart.CreateNewNode(id, x, y, z);
            });
    }

    /// Adds NLevels via uniform bisection; new nodes are orphan.
    void AddLevel(const SizeType NLevels)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBCurveGeometry::AddLevel: no base level defined." << std::endl;
        for (SizeType i = 0; i < NLevels; ++i) {
            AddLevel(BisectKnots(mLevels.back().Knots));
        }
    }

    /// Adds NLevels via uniform bisection; new nodes are registered in rModelPart.
    void AddLevel(const SizeType NLevels, ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBCurveGeometry::AddLevel: no base level defined." << std::endl;
        for (SizeType i = 0; i < NLevels; ++i) {
            AddLevel(BisectKnots(mLevels.back().Knots), rModelPart);
        }
    }

    /**
     * @brief Registers a refinement region on an existing level.
     * @param Level   Index of the level this region belongs to (≥ 1).
     * @param MinT    Lower bound of the refinement interval (inclusive).
     * @param MaxT    Upper bound of the refinement interval (exclusive).
     */
    void AddRefinementDomain(
        const SizeType Level,
        const double MinT, const double MaxT)
    {
        KRATOS_ERROR_IF(Level == 0 || Level >= mLevels.size())
            << "THBCurveGeometry::AddRefinementDomain: level " << Level
            << " does not exist (call AddLevel first)." << std::endl;

        mRefinementDomains.push_back(THBRefinementDomain{Level, MinT, MaxT});
        ComputeActiveFunctions();
        ComputeTruncationData();
    }

    /**
     * @brief Removes inactive control points from the ModelPart and compacts Points().
     *
     * After this call:
     *   - Inactive nodes are removed from rModelPart (if registered there).
     *   - this->Points() contains only active control points, ordered level-by-level.
     *   - CreateQuadraturePointGeometries uses the compacted index layout automatically.
     *
     * Precondition: AddRefinementDomain must have been called at least once.
     * May only be called once; calling again raises an error.
     */
    void EliminateInactiveFunctions(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(mIsEliminated)
            << "THBCurveGeometry::EliminateInactiveFunctions: already eliminated." << std::endl;
        KRATOS_ERROR_IF(mRefinementDomains.empty())
            << "THBCurveGeometry::EliminateInactiveFunctions: no refinement domains defined — "
               "call AddRefinementDomain first." << std::endl;

        const SizeType polynomial_degree = mLevels[0].Degree;
        const SizeType num_levels = mLevels.size();

        mActiveOffset.resize(num_levels);
        mActiveLocalIndex.resize(num_levels);

        PointsArrayType packed_points;
        SizeType running_offset = 0;

        for (SizeType l = 0; l < num_levels; ++l) {
            const SizeType num_cps_at_level =
                mLevels[l].Knots.size() - polynomial_degree + 1;
            const SizeType level_offset = ControlPointOffset(l);

            mActiveOffset[l] = running_offset;
            mActiveLocalIndex[l].assign(num_cps_at_level, -1);

            int packed_count = 0;
            for (SizeType i = 0; i < num_cps_at_level; ++i) {
                auto p_node = this->pGetPoint(level_offset + i);
                if (mActiveFunctions[l][i]) {
                    mActiveLocalIndex[l][i] = packed_count++;
                    packed_points.push_back(p_node);
                } else {
                    if (rModelPart.HasNode(p_node->Id()))
                        rModelPart.RemoveNode(p_node->Id());
                }
            }
            running_offset += static_cast<SizeType>(packed_count);
        }

        this->Points() = packed_points;
        mIsEliminated = true;
    }

    ///@}
    ///@name Getters
    ///@{

    const THBLevelContainer& Levels() const { return mLevels; }

    const RefinementDomainContainer& RefinementDomains() const { return mRefinementDomains; }

    SizeType NumberOfLevels() const { return mLevels.size(); }

    /// Active-function flags for a given level (indexed by flat CP index).
    const std::vector<bool>& GetActiveFunctions(SizeType Level) const
    {
        KRATOS_DEBUG_ERROR_IF(Level >= mActiveFunctions.size())
            << "THBCurveGeometry::GetActiveFunctions: level " << Level
            << " out of range." << std::endl;
        return mActiveFunctions[Level];
    }

    /// Number of active basis functions at the given level.
    SizeType NumberOfActiveFunctions(SizeType Level) const
    {
        KRATOS_DEBUG_ERROR_IF(Level >= mActiveFunctions.size())
            << "THBCurveGeometry::NumberOfActiveFunctions: level " << Level
            << " out of range." << std::endl;
        SizeType count = 0;
        for (bool active : mActiveFunctions[Level])
            if (active) ++count;
        return count;
    }

    /// Truncation coefficients for active coarse CP (Level, FlatIndex).
    const TruncationEntryContainer& GetTruncationData(SizeType Level, SizeType FlatIndex) const
    {
        static const TruncationEntryContainer empty;
        if (Level >= mTruncationData.size()) return empty;
        if (FlatIndex >= mTruncationData[Level].size()) return empty;
        return mTruncationData[Level][FlatIndex];
    }

    /// Evaluates THB shape functions at each parameter value in rTValues.
    /// Returns a (num_eval_pts × num_active_cps) matrix; result(i,j) is the value
    /// of the j-th active CP's basis function at rTValues[i], zero if not in support.
    Matrix EvaluateShapeFunctions(const Vector& rTValues) const
    {
        const SizeType num_eval_pts = rTValues.size();
        const SizeType num_active_cps = this->PointsNumber();
        Matrix result(num_eval_pts, num_active_cps, 0.0);

        THBCurveShapeFunction shape_function_container(PolynomialDegree(0), 0);
        for (SizeType i = 0; i < num_eval_pts; ++i) {
            shape_function_container.ComputeShapeFunctionValues(*this, rTValues[i]);
            const auto& cp_indices = shape_function_container.ControlPointIndices();
            for (SizeType j = 0; j < shape_function_container.NumberOfNonzeroControlPoints(); ++j)
                result(i, cp_indices[j]) = shape_function_container(j, 0);
        }
        return result;
    }

    /// Returns the finest level whose refinement domain contains t, or 0 if none.
    SizeType ActiveLevelAtPoint(double t) const
    {
        SizeType active_level = 0;
        for (SizeType l = 1; l < mLevels.size(); ++l) {
            if (IsInsideDomain(t, l))
                active_level = l;
        }
        return active_level;
    }

    /// Maps (level, flat-index-within-level) to the index in the packed Points() array.
    SizeType PackedControlPointIndex(SizeType Level, SizeType FlatIndex) const
    {
        if (mIsEliminated) {
            KRATOS_ERROR_IF(mActiveLocalIndex[Level][FlatIndex] < 0)
                << "THBCurveGeometry::PackedControlPointIndex: "
                   "inactive CP requested at level " << Level << " flat " << FlatIndex
                << " — the evaluation point used a level whose active CP set is incomplete. "
                   "Ensure refinement domain boundaries align with knot span boundaries." << std::endl;
            return mActiveOffset[Level] + static_cast<SizeType>(mActiveLocalIndex[Level][FlatIndex]);
        } else {
            return ControlPointOffset(Level) + FlatIndex;
        }
    }

    ///@}
    ///@name LocalRefinedCurveGeometry interface
    ///@{

    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        KRATOS_DEBUG_ERROR_IF(mLevels.empty())
            << "THBCurveGeometry::PolynomialDegree: no levels defined." << std::endl;
        return mLevels[0].Degree;
    }

    void SpansLocalSpace(
        std::vector<double>& rSpans,
        IndexType DirectionIndex) const override
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBCurveGeometry::SpansLocalSpace: no levels defined." << std::endl;

        const Vector& knots = mLevels.back().Knots;
        rSpans.clear();
        for (IndexType i = 0; i < static_cast<IndexType>(knots.size()); ++i) {
            if (rSpans.empty() || std::abs(knots[i] - rSpans.back()) > 1e-10)
                rSpans.push_back(knots[i]);
        }
    }

    /**
     * @brief Generates integration points for the hierarchical domain.
     *
     * Each knot span is assigned to the finest level that contains it:
     *   - Level 0: spans NOT covered by Omega^1
     *   - Level l: spans inside Omega^l but NOT inside Omega^{l+1}
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        const SizeType points_per_span = PolynomialDegree(0) + 1;

        rIntegrationPoints.clear();

        IntegrationPointsArrayType span_ips(points_per_span);

        for (SizeType l = 0; l < mLevels.size(); ++l) {
            const auto span_intervals = KnotSpanIntervals(mLevels[l].Knots);

            for (const auto& span : span_intervals) {
                const double span_mid = 0.5 * (span.GetT0() + span.GetT1());

                const bool in_this_level = (l == 0) || IsInsideDomain(span_mid, l);
                const bool in_next_level =
                    (l + 1 < mLevels.size()) && IsInsideDomain(span_mid, l + 1);

                if (!in_this_level || in_next_level) continue;

                typename IntegrationPointsArrayType::iterator integration_point_iterator =
                    span_ips.begin();
                IntegrationPointUtilities::IntegrationPoints1D(
                    integration_point_iterator, points_per_span, span.GetT0(), span.GetT1());

                for (const auto& ip : span_ips)
                    rIntegrationPoints.push_back(ip);
            }
        }
    }

    /**
     * @brief Builds quadrature-point geometries for all given integration points.
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        if (rResultGeometries.size() != rIntegrationPoints.size())
            rResultGeometries.resize(rIntegrationPoints.size());

        auto default_method = this->GetDefaultIntegrationMethod();

        THBCurveShapeFunction shape_function_container(
            PolynomialDegree(0), NumberOfShapeFunctionDerivatives);

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i) {
            const double t = rIntegrationPoints[i][0];

            std::vector<CoordinatesArrayType> global_space_derivatives(2);
            this->GlobalSpaceDerivatives(global_space_derivatives, rIntegrationPoints[i], 1);

            shape_function_container.ComputeShapeFunctionValues(*this, t);

            const SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();
            const auto& cp_indices = shape_function_container.ControlPointIndices();

            PointsArrayType nonzero_control_points(num_nonzero_cps);
            for (IndexType j = 0; j < num_nonzero_cps; ++j)
                nonzero_control_points(j) = this->pGetPoint(cp_indices[j]);

            Matrix N(1, num_nonzero_cps);
            for (IndexType j = 0; j < num_nonzero_cps; ++j)
                N(0, j) = shape_function_container(j, 0);

            DenseVector<Matrix> shape_function_derivatives(NumberOfShapeFunctionDerivatives - 1);
            for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n)
                shape_function_derivatives[n].resize(num_nonzero_cps, 1);

            if (NumberOfShapeFunctionDerivatives > 0) {
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n)
                    for (IndexType j = 0; j < num_nonzero_cps; ++j)
                        shape_function_derivatives[n](j, 0) = shape_function_container(j, n + 1);
            }

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[i], N, shape_function_derivatives);

            rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCurve(
                this->WorkingSpaceDimension(), 1, data_container, nonzero_control_points,
                global_space_derivatives[1][0], global_space_derivatives[1][1], this);
        }
    }

    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates) const override
    {
        THBCurveShapeFunction shape_function_container(PolynomialDegree(0), 0);
        shape_function_container.ComputeShapeFunctionValues(*this, rLocalCoordinates[0]);

        noalias(rResult) = ZeroVector(3);
        const auto& cp_indices = shape_function_container.ControlPointIndices();
        for (SizeType j = 0; j < shape_function_container.NumberOfNonzeroControlPoints(); ++j) {
            const CoordinatesArrayType& cp = this->GetPoint(cp_indices[j]);
            for (SizeType d = 0; d < 3; ++d)
                rResult[d] += shape_function_container(j, 0) * cp[d];
        }
        return rResult;
    }

    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override
    {
        THBCurveShapeFunction shape_function_container(PolynomialDegree(0), DerivativeOrder);
        shape_function_container.ComputeShapeFunctionValues(*this, rLocalCoordinates[0]);

        rGlobalSpaceDerivatives.resize(DerivativeOrder + 1);
        for (SizeType order = 0; order <= DerivativeOrder; ++order)
            noalias(rGlobalSpaceDerivatives[order]) = ZeroVector(3);

        const auto& cp_indices = shape_function_container.ControlPointIndices();
        for (SizeType j = 0; j < shape_function_container.NumberOfNonzeroControlPoints(); ++j) {
            const CoordinatesArrayType& cp = this->GetPoint(cp_indices[j]);
            for (SizeType order = 0; order <= DerivativeOrder; ++order)
                for (SizeType d = 0; d < 3; ++d)
                    rGlobalSpaceDerivatives[order][d] += shape_function_container(j, order) * cp[d];
        }
    }

    ///@}
    ///@name Integration Info
    ///@{

    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        return IntegrationInfo(1,
            PolynomialDegree(0) + 1,
            IntegrationInfo::QuadratureMethod::GAUSS);
    }

    ///@}
    ///@name Information
    ///@{

    std::string Info() const override
    {
        return std::to_string(TWorkingSpaceDimension)
            + "D THBCurveGeometry (" + std::to_string(mLevels.size()) + " level(s))";
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    void PrintData(std::ostream& rOStream) const override
    {
        for (SizeType l = 0; l < mLevels.size(); ++l) {
            rOStream << "  Level " << l
                     << ": p=" << mLevels[l].Degree
                     << " nKnots=" << mLevels[l].Knots.size() << "\n";
        }
        for (const auto& dom : mRefinementDomains) {
            rOStream << "  Omega (level " << dom.Level << ") = ["
                     << dom.MinT << ", " << dom.MaxT << ")\n";
        }
    }

    ///@}

private:
    ///@name Private Member Variables
    ///@{

    THBLevelContainer         mLevels;
    RefinementDomainContainer mRefinementDomains;
    std::vector<std::vector<bool>> mActiveFunctions;

    bool                          mIsEliminated = false;
    std::vector<SizeType>         mActiveOffset;
    std::vector<std::vector<int>> mActiveLocalIndex;

    std::vector<std::vector<TruncationEntryContainer>> mTruncationData;

    ///@}
    ///@name Private Helpers
    ///@{

    /// Shared implementation for all AddLevel overloads.
    template<typename TNodeFactory>
    void AddLevelImpl(
        const Vector& rKnots,
        const Vector& rWeights,
        TNodeFactory CreateNode)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBCurveGeometry::AddLevel: no base level defined." << std::endl;
        KRATOS_ERROR_IF(mIsEliminated)
            << "THBCurveGeometry::AddLevel: cannot add levels after EliminateInactiveFunctions." << std::endl;

        const SizeType last_level_index = mLevels.size() - 1;
        const THBLevel& previous_level  = mLevels[last_level_index];
        const SizeType polynomial_degree = previous_level.Degree;
        const SizeType num_cps_coarse    = previous_level.Knots.size() - polynomial_degree + 1;

        const Matrix M = ComputeRefinementMatrix1D(previous_level.Knots, rKnots, polynomial_degree);
        const SizeType num_cps_fine = M.size1();

        const SizeType coarse_cp_offset = ControlPointOffset(last_level_index);

        IndexType next_node_id = 1;
        for (SizeType i = 0; i < this->PointsNumber(); ++i)
            next_node_id = std::max(next_node_id, this->GetPoint(i).Id() + IndexType(1));

        // P'_j = sum_i M[j, i] * P_i
        for (SizeType j = 0; j < num_cps_fine; ++j) {
            double x = 0.0, y = 0.0, z = 0.0;
            for (SizeType i = 0; i < num_cps_coarse; ++i) {
                const double coeff = M(j, i);
                if (std::abs(coeff) < 1e-15) continue;
                const auto& pt = this->GetPoint(coarse_cp_offset + i);
                x += coeff * pt.X();
                y += coeff * pt.Y();
                z += coeff * pt.Z();
            }
            this->Points().push_back(CreateNode(next_node_id++, x, y, z));
        }

        mLevels.push_back(THBLevel{previous_level.Degree, rKnots, rWeights});
        mActiveFunctions.push_back(std::vector<bool>(num_cps_fine, true));
    }

    /// Inserts the midpoint of every unique span (uniform h-refinement by factor 2).
    static Vector BisectKnots(const Vector& rKnots)
    {
        std::vector<double> result;
        result.reserve(rKnots.size() * 2);
        for (SizeType i = 0; i + 1 < rKnots.size(); ++i) {
            result.push_back(rKnots[i]);
            if (rKnots[i + 1] - rKnots[i] > 1e-10)
                result.push_back(0.5 * (rKnots[i] + rKnots[i + 1]));
        }
        result.push_back(rKnots[rKnots.size() - 1]);
        Vector v(result.size());
        for (SizeType i = 0; i < result.size(); ++i) v[i] = result[i];
        return v;
    }

    /**
     * @brief Computes the 1-D B-spline refinement matrix M (n_fine × n_coarse).
     *
     * P_fine[j] = sum_i M[j, i] * P_coarse[i]
     *
     * Algorithm: iterative single-knot insertion (Boehm's algorithm).
     */
    static Matrix ComputeRefinementMatrix1D(
        const Vector& rCoarseKnots,
        const Vector& rFineKnots,
        const SizeType polynomial_degree)
    {
        const SizeType num_cps_coarse = rCoarseKnots.size() - polynomial_degree + 1;

        // Convert internal → external format: prepend/append first/last knot
        auto to_ext = [](const Vector& k) {
            std::vector<double> ext;
            ext.reserve(k.size() + 2);
            ext.push_back(k[0]);
            for (SizeType i = 0; i < k.size(); ++i) ext.push_back(k[i]);
            ext.push_back(k[k.size() - 1]);
            return ext;
        };

        const std::vector<double> ext_coarse = to_ext(rCoarseKnots);
        const std::vector<double> ext_fine   = to_ext(rFineKnots);

        // new_knots = ext_fine multiset minus ext_coarse multiset (both sorted)
        std::vector<double> new_knots;
        {
            SizeType ci = 0;
            for (SizeType fi = 0; fi < ext_fine.size(); ++fi) {
                if (ci < ext_coarse.size() &&
                    std::abs(ext_fine[fi] - ext_coarse[ci]) < 1e-10)
                    ++ci;
                else
                    new_knots.push_back(ext_fine[fi]);
            }
        }

        // Start with the num_cps_coarse × num_cps_coarse identity; grow by inserting knots.
        Matrix M(num_cps_coarse, num_cps_coarse, 0.0);
        for (SizeType i = 0; i < num_cps_coarse; ++i) M(i, i) = 1.0;

        std::vector<double> current_knots = ext_coarse;
        SizeType num_current_cps = num_cps_coarse;

        for (const double t : new_knots) {
            // knot_span_index: last index with current_knots[k] <= t
            SizeType knot_span_index = 0;
            for (SizeType i = 0; i + 1 < current_knots.size(); ++i)
                if (current_knots[i] <= t + 1e-10) knot_span_index = i;

            // insertion_matrix: (num_current_cps+1) × num_current_cps
            Matrix insertion_matrix(num_current_cps + 1, num_current_cps, 0.0);
            const int span_index = static_cast<int>(knot_span_index);
            const int degree     = static_cast<int>(polynomial_degree);
            for (int i = 0; i <= static_cast<int>(num_current_cps); ++i) {
                if (i <= span_index - degree) {
                    insertion_matrix(i, i) = 1.0;
                } else if (i > span_index) {
                    insertion_matrix(i, i - 1) = 1.0;
                } else {
                    const double knot_span_length =
                        current_knots[i + polynomial_degree] - current_knots[i];
                    const double alpha =
                        (std::abs(knot_span_length) < 1e-15)
                            ? 1.0
                            : (t - current_knots[i]) / knot_span_length;
                    if (i < static_cast<int>(num_current_cps)) insertion_matrix(i, i) = alpha;
                    if (i > 0) insertion_matrix(i, i - 1) = 1.0 - alpha;
                }
            }

            M = prod(insertion_matrix, M);
            current_knots.insert(current_knots.begin() + knot_span_index + 1, t);
            ++num_current_cps;
        }

        return M;
    }

    /// Unique non-zero knot spans from an internal-format knot vector.
    static std::vector<NurbsInterval> KnotSpanIntervals(const Vector& rKnots)
    {
        std::vector<NurbsInterval> result;
        for (SizeType i = 0; i + 1 < rKnots.size(); ++i) {
            if (rKnots[i + 1] - rKnots[i] > 1e-10)
                result.emplace_back(rKnots[i], rKnots[i + 1]);
        }
        return result;
    }

    /// Returns true if t is inside any refinement domain at level >= min_level.
    bool IsInsideRefinedRegion(double t, SizeType min_level) const
    {
        for (const auto& dom : mRefinementDomains)
            if (dom.Level >= min_level && t >= dom.MinT && t <= dom.MaxT)
                return true;
        return false;
    }

    /**
     * @brief Two-phase THB active-function selection (supports the coexistence scenario).
     *
     * For each level pair (l, l+1):
     *
     *   Phase 1 (T^B_{l+1}): activate fine children
     *      Activate fine CP j if supp(B_j^{l+1}) ⊆ Ω^{≥l+1}.
     *      This is independent of parent status → coexistence is possible.
     *
     *   Phase 2 (T^A_{l+1}): deactivate coarse parents
     *      Deactivate active coarse CP i if supp(B_i^l) ⊆ Ω^{≥l+1}.
     *      CPs whose support straddles the boundary are retained (truncated at eval time).
     */
    void ComputeActiveFunctions()
    {
        const SizeType polynomial_degree = mLevels[0].Degree;

        std::fill(mActiveFunctions[0].begin(), mActiveFunctions[0].end(), true);
        for (SizeType l = 1; l < mLevels.size(); ++l)
            std::fill(mActiveFunctions[l].begin(), mActiveFunctions[l].end(), false);

        for (SizeType l = 0; l + 1 < mLevels.size(); ++l) {
            const THBLevel& coarse_level = mLevels[l];
            const THBLevel& fine_level   = mLevels[l + 1];
            const SizeType num_cps_coarse = coarse_level.Knots.size() - polynomial_degree + 1;
            const SizeType num_cps_fine   = fine_level.Knots.size()   - polynomial_degree + 1;

            // Phase 1 (T^B_{l+1}): activate fine children
            for (SizeType j = 0; j < num_cps_fine; ++j) {
                const double fine_supp_min =
                    (j == 0) ? fine_level.Knots[0] : fine_level.Knots[j - 1];
                const double fine_supp_max =
                    fine_level.Knots[std::min(j + polynomial_degree,
                                              fine_level.Knots.size() - 1)];

                bool support_in_omega = true;
                for (SizeType k = 0; k + 1 < fine_level.Knots.size() && support_in_omega; ++k) {
                    if (fine_level.Knots[k + 1] - fine_level.Knots[k] < 1e-10) continue;
                    if (fine_level.Knots[k] < fine_supp_min - 1e-10) continue;
                    if (fine_level.Knots[k + 1] > fine_supp_max + 1e-10) continue;
                    const double cell_mid = 0.5 * (fine_level.Knots[k] + fine_level.Knots[k + 1]);
                    if (!IsInsideRefinedRegion(cell_mid, l + 1))
                        support_in_omega = false;
                }
                if (support_in_omega)
                    mActiveFunctions[l + 1][j] = true;
            }

            // Phase 2 (T^A_{l+1}): deactivate coarse parents
            for (SizeType i = 0; i < num_cps_coarse; ++i) {
                if (!mActiveFunctions[l][i]) continue;

                const double support_min =
                    (i == 0) ? coarse_level.Knots[0] : coarse_level.Knots[i - 1];
                const double support_max =
                    coarse_level.Knots[std::min(i + polynomial_degree,
                                                coarse_level.Knots.size() - 1)];

                bool all_cells_covered = true;
                for (SizeType k = 0; k + 1 < coarse_level.Knots.size() && all_cells_covered; ++k) {
                    if (coarse_level.Knots[k + 1] - coarse_level.Knots[k] < 1e-10) continue;
                    if (coarse_level.Knots[k]     < support_min - 1e-10) continue;
                    if (coarse_level.Knots[k + 1] > support_max + 1e-10) continue;
                    const double cell_mid = 0.5 * (coarse_level.Knots[k] + coarse_level.Knots[k + 1]);
                    if (!IsInsideRefinedRegion(cell_mid, l + 1))
                        all_cells_covered = false;
                }
                if (all_cells_covered)
                    mActiveFunctions[l][i] = false;
            }
        }
    }

    /**
     * @brief Builds multi-level truncation coefficient lists for active coarse CPs.
     *
     * For each active CP i at level l, propagates refinement coefficients through the
     * hierarchy: pending_coefficients coefficients for inactive level-m CPs are forwarded to level
     * m+1 via the refinement relation; pending_coefficientss landing on active CPs are stored as
     * TruncationEntry {m, j, accumulated_coefficient}.
     *
     * This correctly handles 3+ level hierarchies where a level-l CP may need to subtract
     * active functions at levels l+2, l+3, … (not just l+1).
     */
    void ComputeTruncationData()
    {
        const SizeType polynomial_degree = mLevels[0].Degree;
        const SizeType num_levels = mLevels.size();

        // Pre-compute per-level CP counts and consecutive refinement matrices.
        std::vector<SizeType> num_cps(num_levels);
        for (SizeType l = 0; l < num_levels; ++l)
            num_cps[l] = mLevels[l].Knots.size() - polynomial_degree + 1;
        std::vector<Matrix> ref_mat(num_levels - 1);
        for (SizeType l = 0; l + 1 < num_levels; ++l)
            ref_mat[l] = ComputeRefinementMatrix1D(
                mLevels[l].Knots, mLevels[l + 1].Knots, polynomial_degree);

        mTruncationData.resize(num_levels);
        for (SizeType l = 0; l < num_levels; ++l)
            mTruncationData[l].assign(num_cps[l], TruncationEntryContainer{});

        for (SizeType l = 0; l + 1 < num_levels; ++l) {
            for (SizeType coarse_cp = 0; coarse_cp < num_cps[l]; ++coarse_cp) {
                if (!mActiveFunctions[l][coarse_cp]) continue;

                auto& entries = mTruncationData[l][coarse_cp];

                // Initialize pending_coefficients at level l+1 from the l->l+1 refinement relation.
                std::unordered_map<SizeType, double> pending_coefficients;
                for (SizeType fine_cp = 0; fine_cp < num_cps[l + 1]; ++fine_cp) {
                    const double coeff = ref_mat[l](fine_cp, coarse_cp);
                    if (std::abs(coeff) >= 1e-15)
                        pending_coefficients[fine_cp] = coeff;
                }

                // Propagate through levels l+1, l+2, ... stopping at active CPs.
                for (SizeType m = l + 1; m < num_levels; ++m) {
                    std::unordered_map<SizeType, double> next_pending_coefficients;
                    for (const auto& [cp, coeff] : pending_coefficients) {
                        if (std::abs(coeff) < 1e-15) continue;
                        if (mActiveFunctions[m][cp]) {
                            // Active: store as truncation entry.
                            entries.push_back({m, cp, coeff});
                        } else if (m + 1 < num_levels) {
                            // Inactive: propagate through the m->m+1 refinement relation.
                            for (SizeType next_cp = 0; next_cp < num_cps[m + 1]; ++next_cp) {
                                const double ck = ref_mat[m](next_cp, cp);
                                if (std::abs(ck) >= 1e-15)
                                    next_pending_coefficients[next_cp] += coeff * ck;
                            }
                        }
                    }
                    pending_coefficients = std::move(next_pending_coefficients);
                }
            }
        }
    }

    /// Returns true if t is inside the refinement domain registered for the given level.
    /// Half-open interval [MinT, MaxT) to avoid double-counting shared boundaries.
    bool IsInsideDomain(double t, SizeType Level) const
    {
        for (const auto& dom : mRefinementDomains) {
            if (dom.Level == Level && t >= dom.MinT && t < dom.MaxT)
                return true;
        }
        return false;
    }

    SizeType NumberOfControlPoints(SizeType Level) const
    {
        return mLevels[Level].Knots.size() - mLevels[0].Degree + 1;
    }

    SizeType ControlPointOffset(SizeType Level) const
    {
        SizeType offset = 0;
        for (SizeType l = 0; l < Level; ++l)
            offset += NumberOfControlPoints(l);
        return offset;
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    THBCurveGeometry() : BaseType(PointsArrayType()) {}

    ///@}

}; // class THBCurveGeometry

} // namespace Kratos
