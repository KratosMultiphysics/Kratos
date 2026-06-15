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
        const SizeType n = rKnots.size() - PolynomialDegree + 1;
        mActiveFunctions.push_back(std::vector<bool>(n, true));
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
        const SizeType n = rKnots.size() - PolynomialDegree + 1;
        mActiveFunctions.push_back(std::vector<bool>(n, true));
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
        for (SizeType n = 0; n < NLevels; ++n) {
            AddLevel(BisectKnots(mLevels.back().Knots));
        }
    }

    /// Adds NLevels via uniform bisection; new nodes are registered in rModelPart.
    void AddLevel(const SizeType NLevels, ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBCurveGeometry::AddLevel: no base level defined." << std::endl;
        for (SizeType n = 0; n < NLevels; ++n) {
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

        const SizeType p = mLevels[0].Degree;
        const SizeType nLevels = mLevels.size();

        mActiveOffset.resize(nLevels);
        mActiveLocalIndex.resize(nLevels);

        PointsArrayType packed_points;
        SizeType running_offset = 0;

        for (SizeType l = 0; l < nLevels; ++l) {
            const SizeType n          = mLevels[l].Knots.size() - p + 1;
            const SizeType old_offset = ControlPointOffset(l);

            mActiveOffset[l] = running_offset;
            mActiveLocalIndex[l].assign(n, -1);

            int packed = 0;
            for (SizeType i = 0; i < n; ++i) {
                auto p_node = this->pGetPoint(old_offset + i);
                if (mActiveFunctions[l][i]) {
                    mActiveLocalIndex[l][i] = packed++;
                    packed_points.push_back(p_node);
                } else {
                    if (rModelPart.HasNode(p_node->Id()))
                        rModelPart.RemoveNode(p_node->Id());
                }
            }
            running_offset += static_cast<SizeType>(packed);
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

    /// Returns the finest level whose refinement domain contains t, or 0 if none.
    SizeType ActiveLevelAtPoint(double t) const
    {
        SizeType active = 0;
        for (SizeType l = 1; l < mLevels.size(); ++l) {
            if (IsInsideDomain(t, l))
                active = l;
        }
        return active;
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
        const SizeType p      = PolynomialDegree(0);
        const SizeType points = p + 1;

        rIntegrationPoints.clear();

        IntegrationPointsArrayType span_ips(points);

        for (SizeType l = 0; l < mLevels.size(); ++l) {
            const auto spans = KnotSpanIntervals(mLevels[l].Knots);

            for (const auto& s : spans) {
                const double mid = 0.5 * (s.GetT0() + s.GetT1());

                const bool in_this_level = (l == 0) || IsInsideDomain(mid, l);
                const bool in_next_level = (l + 1 < mLevels.size()) && IsInsideDomain(mid, l + 1);

                if (!in_this_level || in_next_level) continue;

                auto it = span_ips.begin();
                IntegrationPointUtilities::IntegrationPoints1D(
                    it, points, s.GetT0(), s.GetT1());

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
        const SizeType p = PolynomialDegree(0);

        THBCurveShapeFunction sf(p, NumberOfShapeFunctionDerivatives);

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i) {
            const double t = rIntegrationPoints[i][0];

            std::vector<CoordinatesArrayType> global_space_derivatives(2);
            this->GlobalSpaceDerivatives(global_space_derivatives, rIntegrationPoints[i], 1);

            sf.ComputeShapeFunctionValues(*this, t);

            const SizeType num_nonzero = sf.NumberOfNonzeroControlPoints();
            const auto& cp_indices     = sf.ControlPointIndices();

            PointsArrayType nonzero_control_points(num_nonzero);
            for (IndexType j = 0; j < num_nonzero; ++j)
                nonzero_control_points(j) = this->pGetPoint(cp_indices[j]);

            Matrix N(1, num_nonzero);
            for (IndexType j = 0; j < num_nonzero; ++j)
                N(0, j) = sf(j, 0);

            DenseVector<Matrix> sf_derivatives(NumberOfShapeFunctionDerivatives - 1);
            for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n)
                sf_derivatives[n].resize(num_nonzero, 1);

            if (NumberOfShapeFunctionDerivatives > 0) {
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n)
                    for (IndexType j = 0; j < num_nonzero; ++j)
                        sf_derivatives[n](j, 0) = sf(j, n + 1);
            }

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[i], N, sf_derivatives);

            rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCurve(
                this->WorkingSpaceDimension(), 1, data_container, nonzero_control_points,
                global_space_derivatives[1][0], global_space_derivatives[1][1], this);
        }
    }

    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates) const override
    {
        THBCurveShapeFunction sf(PolynomialDegree(0), 0);
        sf.ComputeShapeFunctionValues(*this, rLocalCoordinates[0]);

        noalias(rResult) = ZeroVector(3);
        const auto& cp_indices = sf.ControlPointIndices();
        for (SizeType j = 0; j < sf.NumberOfNonzeroControlPoints(); ++j) {
            const CoordinatesArrayType& cp = this->GetPoint(cp_indices[j]);
            for (SizeType d = 0; d < 3; ++d)
                rResult[d] += sf(j, 0) * cp[d];
        }
        return rResult;
    }

    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override
    {
        THBCurveShapeFunction sf(PolynomialDegree(0), DerivativeOrder);
        sf.ComputeShapeFunctionValues(*this, rLocalCoordinates[0]);

        rGlobalSpaceDerivatives.resize(DerivativeOrder + 1);
        for (SizeType order = 0; order <= DerivativeOrder; ++order)
            noalias(rGlobalSpaceDerivatives[order]) = ZeroVector(3);

        const auto& cp_indices = sf.ControlPointIndices();
        for (SizeType j = 0; j < sf.NumberOfNonzeroControlPoints(); ++j) {
            const CoordinatesArrayType& cp = this->GetPoint(cp_indices[j]);
            for (SizeType order = 0; order <= DerivativeOrder; ++order)
                for (SizeType d = 0; d < 3; ++d)
                    rGlobalSpaceDerivatives[order][d] += sf(j, order) * cp[d];
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

        const SizeType l = mLevels.size() - 1;
        const THBLevel& prev = mLevels[l];
        const SizeType p     = prev.Degree;
        const SizeType n_old = prev.Knots.size() - p + 1;

        const Matrix M = ComputeRefinementMatrix1D(prev.Knots, rKnots, p);
        const SizeType n_new = M.size1();

        const SizeType offset = ControlPointOffset(l);

        IndexType new_id = 1;
        for (SizeType i = 0; i < this->PointsNumber(); ++i)
            new_id = std::max(new_id, this->GetPoint(i).Id() + IndexType(1));

        // P'_j = sum_i M[j, i] * P_i
        for (SizeType j = 0; j < n_new; ++j) {
            double x = 0.0, y = 0.0, z = 0.0;
            for (SizeType i = 0; i < n_old; ++i) {
                const double coeff = M(j, i);
                if (std::abs(coeff) < 1e-15) continue;
                const auto& pt = this->GetPoint(offset + i);
                x += coeff * pt.X();
                y += coeff * pt.Y();
                z += coeff * pt.Z();
            }
            this->Points().push_back(CreateNode(new_id++, x, y, z));
        }

        mLevels.push_back(THBLevel{prev.Degree, rKnots, rWeights});
        mActiveFunctions.push_back(std::vector<bool>(n_new, true));
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
        const SizeType p)
    {
        const SizeType n_c = rCoarseKnots.size() - p + 1;

        auto to_ext = [](const Vector& k) {
            std::vector<double> ext;
            ext.reserve(k.size() + 2);
            ext.push_back(k[0]);
            for (SizeType i = 0; i < k.size(); ++i) ext.push_back(k[i]);
            ext.push_back(k[k.size() - 1]);
            return ext;
        };

        const std::vector<double> ext_c = to_ext(rCoarseKnots);
        const std::vector<double> ext_f = to_ext(rFineKnots);

        std::vector<double> new_knots;
        {
            SizeType ci = 0;
            for (SizeType fi = 0; fi < ext_f.size(); ++fi) {
                if (ci < ext_c.size() && std::abs(ext_f[fi] - ext_c[ci]) < 1e-10)
                    ++ci;
                else
                    new_knots.push_back(ext_f[fi]);
            }
        }

        Matrix M(n_c, n_c, 0.0);
        for (SizeType i = 0; i < n_c; ++i) M(i, i) = 1.0;

        std::vector<double> cur = ext_c;
        SizeType n_cur = n_c;

        for (const double t : new_knots) {
            SizeType k = 0;
            for (SizeType i = 0; i + 1 < cur.size(); ++i)
                if (cur[i] <= t + 1e-10) k = i;

            Matrix T(n_cur + 1, n_cur, 0.0);
            const int ik = static_cast<int>(k);
            const int ip = static_cast<int>(p);
            for (int i = 0; i <= static_cast<int>(n_cur); ++i) {
                if (i <= ik - ip) {
                    T(i, i) = 1.0;
                } else if (i > ik) {
                    T(i, i - 1) = 1.0;
                } else {
                    const double denom = cur[i + p] - cur[i];
                    const double alpha = (std::abs(denom) < 1e-15) ? 1.0 : (t - cur[i]) / denom;
                    if (i < static_cast<int>(n_cur)) T(i, i)     = alpha;
                    if (i > 0)                       T(i, i - 1) = 1.0 - alpha;
                }
            }

            M = prod(T, M);
            cur.insert(cur.begin() + k + 1, t);
            ++n_cur;
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
     * @brief Top-down propagation that assigns active/inactive flags.
     *
     * A coarse CP at level l is deactivated only if ALL knot cells of its
     * support lie inside a refinement domain at level ≥ l+1.
     * Its true children (those whose support is contained within the coarse support)
     * are then activated at level l+1.
     */
    void ComputeActiveFunctions()
    {
        const SizeType p = mLevels[0].Degree;

        std::fill(mActiveFunctions[0].begin(), mActiveFunctions[0].end(), true);
        for (SizeType l = 1; l < mLevels.size(); ++l)
            std::fill(mActiveFunctions[l].begin(), mActiveFunctions[l].end(), false);

        for (SizeType l = 0; l + 1 < mLevels.size(); ++l) {
            const THBLevel& lev_c = mLevels[l];
            const THBLevel& lev_f = mLevels[l + 1];
            const SizeType  n_c   = lev_c.Knots.size() - p + 1;
            const SizeType  n_f   = lev_f.Knots.size() - p + 1;

            for (SizeType i = 0; i < n_c; ++i) {
                if (!mActiveFunctions[l][i]) continue;

                // Support of B_i^l in internal knot format
                const double t_min = (i == 0) ? lev_c.Knots[0] : lev_c.Knots[i - 1];
                const double t_max = lev_c.Knots[std::min(i + p, lev_c.Knots.size() - 1)];

                // Is every cell of this support inside Ω^{≥l+1}?
                bool all_cells_covered = true;
                for (SizeType k = 0; k + 1 < lev_c.Knots.size() && all_cells_covered; ++k) {
                    if (lev_c.Knots[k + 1] - lev_c.Knots[k] < 1e-10) continue;
                    if (lev_c.Knots[k] < t_min - 1e-10 || lev_c.Knots[k + 1] > t_max + 1e-10) continue;
                    const double mid = 0.5 * (lev_c.Knots[k] + lev_c.Knots[k + 1]);
                    if (!IsInsideRefinedRegion(mid, l + 1))
                        all_cells_covered = false;
                }
                if (!all_cells_covered) continue;

                // Deactivate coarse CP and activate its true children.
                mActiveFunctions[l][i] = false;

                for (SizeType j = 0; j < n_f; ++j) {
                    const double ft_min = (j == 0) ? lev_f.Knots[0] : lev_f.Knots[j - 1];
                    const double ft_max = lev_f.Knots[std::min(j + p, lev_f.Knots.size() - 1)];
                    if (ft_min < t_min - 1e-10 || ft_max > t_max + 1e-10) continue;
                    mActiveFunctions[l + 1][j] = true;
                }
            }
        }
    }

    /**
     * @brief Builds truncation coefficient lists for active coarse CPs.
     *
     * For active CP i at level l, stores entries {j, M[j,i]} for every
     * active fine CP j at level l+1 with nonzero refinement coefficient.
     */
    void ComputeTruncationData()
    {
        const SizeType p  = mLevels[0].Degree;
        const SizeType nL = mLevels.size();

        mTruncationData.resize(nL);

        for (SizeType l = 0; l + 1 < nL; ++l) {
            const THBLevel& lev_c = mLevels[l];
            const THBLevel& lev_f = mLevels[l + 1];
            const SizeType  n_c   = lev_c.Knots.size() - p + 1;
            const SizeType  n_f   = lev_f.Knots.size() - p + 1;

            const Matrix M = ComputeRefinementMatrix1D(lev_c.Knots, lev_f.Knots, p);

            mTruncationData[l].assign(n_c, TruncationEntryContainer{});

            for (SizeType i = 0; i < n_c; ++i) {
                if (!mActiveFunctions[l][i]) continue;

                auto& entries = mTruncationData[l][i];
                for (SizeType j = 0; j < n_f; ++j) {
                    const double coeff = M(j, i);
                    if (std::abs(coeff) < 1e-15) continue;
                    if (!mActiveFunctions[l + 1][j]) continue;
                    entries.push_back({j, coeff});
                }
            }
        }

        // Finest level: no truncation.
        const SizeType n_last = mLevels[nL - 1].Knots.size() - p + 1;
        mTruncationData[nL - 1].assign(n_last, TruncationEntryContainer{});
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
        const SizeType p = mLevels[0].Degree;
        return mLevels[Level].Knots.size() - p + 1;
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
