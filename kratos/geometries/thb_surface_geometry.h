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
#include "geometries/local_refined_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/thb_surface_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_utilities.h"
#include "utilities/quadrature_points_utility.h"
#include "integration/integration_point_utilities.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TContainerPointType>
class THBSurfaceGeometry
    : public LocalRefinedSurfaceGeometry<TWorkingSpaceDimension, TContainerPointType>
{
public:
    ///@name Type Definitions
    ///@{

    typedef LocalRefinedSurfaceGeometry<TWorkingSpaceDimension, TContainerPointType> BaseType;

    typedef typename BaseType::NodeType                   NodeType;
    typedef typename BaseType::IndexType                  IndexType;
    typedef typename BaseType::SizeType                   SizeType;
    typedef typename BaseType::CoordinatesArrayType       CoordinatesArrayType;
    typedef typename BaseType::PointsArrayType            PointsArrayType;
    typedef typename BaseType::GeometriesArrayType        GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    struct THBLevel {
        SizeType DegreeU;
        SizeType DegreeV;
        Vector   KnotsU;
        Vector   KnotsV;
        Vector   Weights;
    };

    struct THBRefinementDomain {
        SizeType Level;
        double MinU;
        double MaxU;
        double MinV;
        double MaxV;
    };

    struct TruncationEntry {
        SizeType FineFlatIndex;
        double   Coefficient;
    };

    using THBLevelContainer         = std::vector<THBLevel>;
    using RefinementDomainContainer = std::vector<THBRefinementDomain>;
    using TruncationEntryContainer  = std::vector<TruncationEntry>;

    KRATOS_CLASS_POINTER_DEFINITION(THBSurfaceGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructs the initial level-0 B-Spline surface.
    THBSurfaceGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const Vector& rKnotsU,
        const Vector& rKnotsV)
        : BaseType(rThisPoints)
    {
        mLevels.push_back(THBLevel{PolynomialDegreeU, PolynomialDegreeV, rKnotsU, rKnotsV, {}});
        const SizeType num_cps_u = rKnotsU.size() - PolynomialDegreeU + 1;
        const SizeType num_cps_v = rKnotsV.size() - PolynomialDegreeV + 1;
        mActiveFunctions.push_back(std::vector<bool>(num_cps_u * num_cps_v, true));
    }

    /// Constructs the initial level-0 NURBS surface (with weights).
    THBSurfaceGeometry(
        const PointsArrayType& rThisPoints,
        const SizeType PolynomialDegreeU,
        const SizeType PolynomialDegreeV,
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rWeights)
        : BaseType(rThisPoints)
    {
        KRATOS_ERROR_IF(rWeights.size() != rThisPoints.size())
            << "Number of control points and weights do not match!" << std::endl;

        mLevels.push_back(THBLevel{PolynomialDegreeU, PolynomialDegreeV, rKnotsU, rKnotsV, rWeights});
        const SizeType num_cps_u = rKnotsU.size() - PolynomialDegreeU + 1;
        const SizeType num_cps_v = rKnotsV.size() - PolynomialDegreeV + 1;
        mActiveFunctions.push_back(std::vector<bool>(num_cps_u * num_cps_v, true));
    }

    ~THBSurfaceGeometry() override = default;

    ///@}
    ///@name Refinement
    ///@{

    /// Adds one level with explicit knot vectors. CPs computed via knot insertion;
    /// new nodes are orphan (not registered in any ModelPart).
    void AddLevel(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rWeights = Vector())
    {
        AddLevelImpl(rKnotsU, rKnotsV, rWeights,
            [](IndexType id, double x, double y, double z) {
                return Kratos::make_intrusive<NodeType>(id, x, y, z);
            });
    }

    /// Adds one level with explicit knot vectors. CPs computed via knot insertion;
    /// new nodes are created via rModelPart.CreateNewNode and registered there.
    void AddLevel(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        ModelPart& rModelPart,
        const Vector& rWeights = Vector())
    {
        AddLevelImpl(rKnotsU, rKnotsV, rWeights,
            [&rModelPart](IndexType id, double x, double y, double z)
                -> typename NodeType::Pointer {
                return rModelPart.CreateNewNode(id, x, y, z);
            });
    }

    /// Adds NLevels via uniform bisection; new nodes are orphan.
    void AddLevel(const SizeType NLevels)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::AddLevel: no base level defined." << std::endl;
        for (SizeType i = 0; i < NLevels; ++i) {
            const THBLevel& previous_level = mLevels.back();
            AddLevel(BisectKnots(previous_level.KnotsU), BisectKnots(previous_level.KnotsV));
        }
    }

    /// Adds NLevels via uniform bisection; new nodes are registered in rModelPart.
    void AddLevel(const SizeType NLevels, ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::AddLevel: no base level defined." << std::endl;
        for (SizeType i = 0; i < NLevels; ++i) {
            const THBLevel& previous_level = mLevels.back();
            AddLevel(BisectKnots(previous_level.KnotsU), BisectKnots(previous_level.KnotsV), rModelPart);
        }
    }

    /**
     * @brief Adds a refinement region on an existing level.
     *
     * Can be called multiple times for the same level to register
     * multiple disjoint refinement regions.
     *
     * @param Level  Index of the level this region belongs to (≥ 1).
     */
    void AddRefinementDomain(
        const SizeType Level,
        const double OmegaMinU, const double OmegaMaxU,
        const double OmegaMinV, const double OmegaMaxV)
    {
        KRATOS_ERROR_IF(Level == 0 || Level >= mLevels.size())
            << "THBSurfaceGeometry::AddRefinementDomain: level " << Level
            << " does not exist (call AddLevel first)." << std::endl;

        mRefinementDomains.push_back(THBRefinementDomain{Level, OmegaMinU, OmegaMaxU, OmegaMinV, OmegaMaxV});
        ComputeActiveFunctions();
        ComputeTruncationData();
    }

    /**
     * @brief Removes inactive control points from the ModelPart and compacts Points().
     *
     * After this call:
     *   - Inactive nodes are removed from rModelPart (if they are registered there).
     *   - this->Points() contains only active control points, ordered level-by-level.
     *   - CreateQuadraturePointGeometries uses the compacted index layout automatically.
     *
     * Preconditions: AddRefinementDomain must have been called at least once.
     * This method may only be called once; calling it again raises an error.
     * AddLevel must not be called after this method.
     */
    void EliminateInactiveFunctions(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(mIsEliminated)
            << "THBSurfaceGeometry::EliminateInactiveFunctions: already eliminated." << std::endl;
        KRATOS_ERROR_IF(mRefinementDomains.empty())
            << "THBSurfaceGeometry::EliminateInactiveFunctions: no refinement domains defined — "
               "call AddRefinementDomain first." << std::endl;

        const SizeType polynomial_degree_u = mLevels[0].DegreeU;
        const SizeType polynomial_degree_v = mLevels[0].DegreeV;
        const SizeType num_levels = mLevels.size();

        mActiveOffset.resize(num_levels);
        mActiveLocalIndex.resize(num_levels);

        PointsArrayType packed_points;
        SizeType running_offset = 0;

        for (SizeType l = 0; l < num_levels; ++l) {
            const SizeType num_cps_u_at_level =
                mLevels[l].KnotsU.size() - polynomial_degree_u + 1;
            const SizeType num_cps_v_at_level =
                mLevels[l].KnotsV.size() - polynomial_degree_v + 1;
            const SizeType num_cps_at_level = num_cps_u_at_level * num_cps_v_at_level;
            const SizeType level_offset = ControlPointOffset(l);

            mActiveOffset[l] = running_offset;
            mActiveLocalIndex[l].assign(num_cps_at_level, -1);

            int packed_count = 0;
            for (SizeType flat = 0; flat < num_cps_at_level; ++flat) {
                auto p_node = this->pGetPoint(level_offset + flat);
                if (mActiveFunctions[l][flat]) {
                    mActiveLocalIndex[l][flat] = packed_count++;
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

    /**
     * @brief Returns the active-function flags for a given level.
     *
     * mActiveFunctions[l][i_v * nU + i_u] is true when basis function (i_u, i_v)
     * at level l is active (its support is NOT fully covered by finer refinement
     * domains) and false when it has been eliminated.
     *
     * Recomputed automatically every time AddRefinementDomain is called.
     */
    const std::vector<bool>& GetActiveFunctions(SizeType Level) const
    {
        KRATOS_DEBUG_ERROR_IF(Level >= mActiveFunctions.size())
            << "THBSurfaceGeometry::GetActiveFunctions: level " << Level
            << " out of range." << std::endl;
        return mActiveFunctions[Level];
    }

    /// Returns the number of active basis functions at the given level.
    SizeType NumberOfActiveFunctions(SizeType Level) const
    {
        KRATOS_DEBUG_ERROR_IF(Level >= mActiveFunctions.size())
            << "THBSurfaceGeometry::NumberOfActiveFunctions: level " << Level
            << " out of range." << std::endl;
        SizeType count = 0;
        for (bool active : mActiveFunctions[Level])
            if (active) ++count;
        return count;
    }

    /// Returns the truncation entries for active coarse CP at (Level, FlatIndex).
    /// Empty when the CP has no active children — no truncation needed.
    const TruncationEntryContainer& GetTruncationData(SizeType Level, SizeType FlatIndex) const
    {
        static const TruncationEntryContainer empty;
        if (Level >= mTruncationData.size()) return empty;
        if (FlatIndex >= mTruncationData[Level].size()) return empty;
        return mTruncationData[Level][FlatIndex];
    }

    /// Returns the finest level whose refinement domain contains (u, v), or 0 if none.
    SizeType ActiveLevelAtPoint(double u, double v) const
    {
        SizeType active_level = 0;
        for (SizeType l = 1; l < mLevels.size(); ++l) {
            if (IsInsideDomain(u, v, l))
                active_level = l;
        }
        return active_level;
    }

    /// Maps (level, flat-index-within-level) to the index in the packed Points() array.
    /// Works both before and after EliminateInactiveFunctions.
    SizeType PackedControlPointIndex(SizeType Level, SizeType FlatIndex) const
    {
        if (mIsEliminated) {
            KRATOS_ERROR_IF(mActiveLocalIndex[Level][FlatIndex] < 0)
                << "THBSurfaceGeometry::PackedControlPointIndex: "
                   "inactive CP requested at level " << Level << " flat " << FlatIndex
                << " -- the evaluation point used a level whose active CP set is incomplete. "
                   "Ensure refinement domain boundaries align with knot span boundaries." << std::endl;
            return mActiveOffset[Level] + static_cast<SizeType>(mActiveLocalIndex[Level][FlatIndex]);
        } else {
            return ControlPointOffset(Level) + FlatIndex;
        }
    }

    ///@}
    ///@name LocalRefinedSurfaceGeometry interface
    ///@{

    /// Polynomial degree is the same on all THB levels; read from level 0.
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        KRATOS_DEBUG_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::PolynomialDegree: no levels defined." << std::endl;

        return (LocalDirectionIndex == 0) ? mLevels[0].DegreeU : mLevels[0].DegreeV;
    }

    void SpansLocalSpace(
        std::vector<double>& rSpans,
        IndexType DirectionIndex) const override
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::SpansLocalSpace: no levels defined." << std::endl;

        // The finest level's knot vector covers the full domain and is a superset of
        // all coarser levels' breakpoints (each level is a uniform bisection of the previous).
        const Vector& knots = (DirectionIndex == 0)
            ? mLevels.back().KnotsU
            : mLevels.back().KnotsV;

        rSpans.clear();
        for (IndexType i = 0; i < static_cast<IndexType>(knots.size()); ++i) {
            if (rSpans.empty() || std::abs(knots[i] - rSpans.back()) > 1e-10)
                rSpans.push_back(knots[i]);
        }
    }

    /**
     * @brief Generates integration points for the hierarchical domain.
     *
     * Integration cells are assigned level by level:
     *   - Level 0: spans NOT covered by Omega^1
     *   - Level l: spans inside Omega^l but NOT inside Omega^{l+1}
     */
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        const SizeType points_in_u = PolynomialDegree(0) + 1;
        const SizeType points_in_v = PolynomialDegree(1) + 1;

        rIntegrationPoints.clear();

        IntegrationPointsArrayType span_ips(points_in_u * points_in_v);

        for (SizeType l = 0; l < mLevels.size(); ++l) {
            const auto knot_span_intervals_u = KnotSpanIntervals(mLevels[l].KnotsU);
            const auto knot_span_intervals_v = KnotSpanIntervals(mLevels[l].KnotsV);

            for (const auto& span_u : knot_span_intervals_u) {
                for (const auto& span_v : knot_span_intervals_v) {
                    const double span_mid_u = 0.5 * (span_u.GetT0() + span_u.GetT1());
                    const double span_mid_v = 0.5 * (span_v.GetT0() + span_v.GetT1());

                    // Active at level l: inside Omega^l (or level 0), AND not refined further
                    const bool in_this_level =
                        (l == 0) || IsInsideDomain(span_mid_u, span_mid_v, l);
                    const bool in_next_level =
                        (l + 1 < mLevels.size()) &&
                        IsInsideDomain(span_mid_u, span_mid_v, l + 1);

                    if (!in_this_level || in_next_level)
                        continue;

                    typename IntegrationPointsArrayType::iterator integration_point_iterator =
                        span_ips.begin();
                    IntegrationPointUtilities::IntegrationPoints2D(
                        integration_point_iterator,
                        points_in_u, points_in_v,
                        span_u.GetT0(), span_u.GetT1(),
                        span_v.GetT0(), span_v.GetT1());

                    for (const auto& ip : span_ips)
                        rIntegrationPoints.push_back(ip);
                }
            }
        }
    }

    /**
     * @brief Builds quadrature-point geometries for all given integration points.
     *
     * For each point (u, v), the finest active level is determined and
     * THBSurfaceShapeFunction evaluates the B-spline shape functions at that
     * level.  The resulting nonzero control points are looked up from the packed
     * Points() array (works both before and after EliminateInactiveFunctions).
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

        THBSurfaceShapeFunction shape_function_container(
            PolynomialDegree(0), PolynomialDegree(1), NumberOfShapeFunctionDerivatives);

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i) {
            const double u = rIntegrationPoints[i][0];
            const double v = rIntegrationPoints[i][1];

            shape_function_container.ComputeShapeFunctionValues(*this, u, v);

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
                shape_function_derivatives[n].resize(num_nonzero_cps, n + 2);

            if (NumberOfShapeFunctionDerivatives > 0) {
                IndexType shape_derivative_index = 1;
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n) {
                    for (IndexType k = 0; k < n + 2; ++k)
                        for (IndexType j = 0; j < num_nonzero_cps; ++j)
                            shape_function_derivatives[n](j, k) =
                                shape_function_container(j, shape_derivative_index + k);
                    shape_derivative_index += n + 2;
                }
            }

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[i], N, shape_function_derivatives);

            rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePoint(
                this->WorkingSpaceDimension(), 2, data_container, nonzero_control_points, this);
        }
    }

    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates) const override
    {
        THBSurfaceShapeFunction shape_function_container(
            PolynomialDegree(0), PolynomialDegree(1), 0);
        shape_function_container.ComputeShapeFunctionValues(
            *this, rLocalCoordinates[0], rLocalCoordinates[1]);

        noalias(rResult) = ZeroVector(3);
        const auto& cp_indices = shape_function_container.ControlPointIndices();
        for (SizeType j = 0; j < shape_function_container.NumberOfNonzeroControlPoints(); ++j) {
            const CoordinatesArrayType& cp = this->GetPoint(cp_indices[j]);
            for (SizeType d = 0; d < 3; ++d)
                rResult[d] += shape_function_container(j, 0) * cp[d];
        }
        return rResult;
    }

    /// @todo Implement using THBSurfaceShapeFunction.
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const override
    {
        KRATOS_ERROR << "THBSurfaceGeometry::GlobalSpaceDerivatives: not yet implemented." << std::endl;
    }

    ///@}
    ///@name Information
    ///@{

    std::string Info() const override
    {
        return std::to_string(TWorkingSpaceDimension)
            + "D THBSurfaceGeometry (" + std::to_string(mLevels.size()) + " level(s))";
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    void PrintData(std::ostream& rOStream) const override
    {
        for (SizeType l = 0; l < mLevels.size(); ++l) {
            rOStream << "  Level " << l
                     << ": p=" << mLevels[l].DegreeU << " q=" << mLevels[l].DegreeV
                     << " nKU=" << mLevels[l].KnotsU.size()
                     << " nKV=" << mLevels[l].KnotsV.size() << "\n";
        }
        for (const auto& dom : mRefinementDomains) {
            rOStream << "  Omega (level " << dom.Level << ") = ["
                     << dom.MinU << "," << dom.MaxU << "] x ["
                     << dom.MinV << "," << dom.MaxV << "]\n";
        }
    }

    ///@}

private:
    ///@name Private Member Variables
    ///@{

    THBLevelContainer         mLevels;
    RefinementDomainContainer mRefinementDomains;
    std::vector<std::vector<bool>> mActiveFunctions;

    bool                           mIsEliminated = false;
    std::vector<SizeType>          mActiveOffset;
    std::vector<std::vector<int>>  mActiveLocalIndex;

    std::vector<std::vector<TruncationEntryContainer>> mTruncationData;

    ///@}
    ///@name Private Helpers
    ///@{

    /// Shared implementation for all AddLevel overloads.
    /// TNodeFactory: callable (IndexType id, double x, double y, double z) -> NodeType::Pointer
    template<typename TNodeFactory>
    void AddLevelImpl(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rWeights,
        TNodeFactory CreateNode)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::AddLevel: no base level defined." << std::endl;
        KRATOS_ERROR_IF(mIsEliminated)
            << "THBSurfaceGeometry::AddLevel: cannot add levels after EliminateInactiveFunctions." << std::endl;

        const SizeType last_level_index   = mLevels.size() - 1;
        const THBLevel& previous_level    = mLevels[last_level_index];
        const SizeType polynomial_degree_u = previous_level.DegreeU;
        const SizeType polynomial_degree_v = previous_level.DegreeV;
        const SizeType num_cps_u_coarse    = previous_level.KnotsU.size() - polynomial_degree_u + 1;
        const SizeType num_cps_v_coarse    = previous_level.KnotsV.size() - polynomial_degree_v + 1;

        const Matrix M_U = ComputeRefinementMatrix1D(previous_level.KnotsU, rKnotsU, polynomial_degree_u);
        const Matrix M_V = ComputeRefinementMatrix1D(previous_level.KnotsV, rKnotsV, polynomial_degree_v);
        const SizeType num_cps_u_fine = M_U.size1();
        const SizeType num_cps_v_fine = M_V.size1();

        const SizeType coarse_cp_offset = ControlPointOffset(last_level_index);

        IndexType next_node_id = 1;
        for (SizeType i = 0; i < this->PointsNumber(); ++i)
            next_node_id = std::max(next_node_id, this->GetPoint(i).Id() + IndexType(1));

        // P'_{j_v, j_u} = sum_{i_v, i_u} M_V[j_v, i_v] * M_U[j_u, i_u] * P_{i_v, i_u}
        for (SizeType j_v = 0; j_v < num_cps_v_fine; ++j_v) {
            for (SizeType j_u = 0; j_u < num_cps_u_fine; ++j_u) {
                double x = 0.0, y = 0.0, z = 0.0;
                for (SizeType i_v = 0; i_v < num_cps_v_coarse; ++i_v) {
                    const double coeff_v = M_V(j_v, i_v);
                    if (std::abs(coeff_v) < 1e-15) continue;
                    for (SizeType i_u = 0; i_u < num_cps_u_coarse; ++i_u) {
                        const double coeff = coeff_v * M_U(j_u, i_u);
                        if (std::abs(coeff) < 1e-15) continue;
                        const auto& pt = this->GetPoint(
                            coarse_cp_offset + i_v * num_cps_u_coarse + i_u);
                        x += coeff * pt.X();
                        y += coeff * pt.Y();
                        z += coeff * pt.Z();
                    }
                }
                this->Points().push_back(CreateNode(next_node_id++, x, y, z));
            }
        }

        mLevels.push_back(THBLevel{
            previous_level.DegreeU, previous_level.DegreeV, rKnotsU, rKnotsV, rWeights});
        mActiveFunctions.push_back(
            std::vector<bool>(num_cps_u_fine * num_cps_v_fine, true));
    }

    /// Returns a new internal-format knot vector with the midpoint of every unique
    /// span inserted (uniform bisection / h-refinement by factor 2).
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
     * Given internal-format knot vectors for a coarse and a fine level (fine must be
     * a refinement of coarse), returns the matrix M such that:
     *   P_fine[j] = sum_i  M[j, i] * P_coarse[i]
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
                    std::abs(ext_fine[fi] - ext_coarse[ci]) < 1e-10) {
                    ++ci;
                } else {
                    new_knots.push_back(ext_fine[fi]);
                }
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
                    // Blend: alpha * P_i + (1-alpha) * P_{i-1}
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

        return M;  // size n_fine × n_coarse
    }

    /// Returns unique non-zero knot span intervals from an internal-format knot vector.
    static std::vector<NurbsInterval> KnotSpanIntervals(const Vector& rKnots)
    {
        std::vector<NurbsInterval> result;
        for (SizeType i = 0; i + 1 < rKnots.size(); ++i) {
            if (rKnots[i + 1] - rKnots[i] > 1e-10)
                result.emplace_back(rKnots[i], rKnots[i + 1]);
        }
        return result;
    }

    /// Returns true if (u, v) is inside any refinement domain at level >= min_level.
    bool IsInsideRefinedRegion(double u, double v, SizeType min_level) const
    {
        for (const auto& dom : mRefinementDomains)
            if (dom.Level >= min_level &&
                u >= dom.MinU && u <= dom.MaxU &&
                v >= dom.MinV && v <= dom.MaxV)
                return true;
        return false;
    }

    /**
     * @brief Builds truncation coefficient lists for active coarse CPs.
     *
     * For each active CP i at level l, stores entries {j_flat, M_U[j_u,i_u]*M_V[j_v,i_v]}
     * for every active fine CP j at level l+1 with nonzero refinement coefficient.
     */
    void ComputeTruncationData()
    {
        const SizeType polynomial_degree_u = mLevels[0].DegreeU;
        const SizeType polynomial_degree_v = mLevels[0].DegreeV;
        const SizeType num_levels = mLevels.size();

        mTruncationData.resize(num_levels);

        for (SizeType l = 0; l + 1 < num_levels; ++l) {
            const THBLevel& coarse_level = mLevels[l];
            const THBLevel& fine_level   = mLevels[l + 1];
            const SizeType num_cps_u_coarse = coarse_level.KnotsU.size() - polynomial_degree_u + 1;
            const SizeType num_cps_v_coarse = coarse_level.KnotsV.size() - polynomial_degree_v + 1;
            const SizeType num_cps_u_fine   = fine_level.KnotsU.size()   - polynomial_degree_u + 1;
            const SizeType num_cps_v_fine   = fine_level.KnotsV.size()   - polynomial_degree_v + 1;

            const Matrix M_U = ComputeRefinementMatrix1D(
                coarse_level.KnotsU, fine_level.KnotsU, polynomial_degree_u);
            const Matrix M_V = ComputeRefinementMatrix1D(
                coarse_level.KnotsV, fine_level.KnotsV, polynomial_degree_v);

            mTruncationData[l].assign(
                num_cps_u_coarse * num_cps_v_coarse, TruncationEntryContainer{});

            for (SizeType i_v = 0; i_v < num_cps_v_coarse; ++i_v) {
                for (SizeType i_u = 0; i_u < num_cps_u_coarse; ++i_u) {
                    const SizeType i_flat = i_v * num_cps_u_coarse + i_u;
                    if (!mActiveFunctions[l][i_flat]) continue;

                    auto& entries = mTruncationData[l][i_flat];
                    for (SizeType j_v = 0; j_v < num_cps_v_fine; ++j_v) {
                        const double coeff_v = M_V(j_v, i_v);
                        if (std::abs(coeff_v) < 1e-15) continue;
                        for (SizeType j_u = 0; j_u < num_cps_u_fine; ++j_u) {
                            const double coeff_u = M_U(j_u, i_u);
                            if (std::abs(coeff_u) < 1e-15) continue;
                            const SizeType j_flat = j_v * num_cps_u_fine + j_u;
                            if (!mActiveFunctions[l + 1][j_flat]) continue;
                            entries.push_back({j_flat, coeff_u * coeff_v});
                        }
                    }
                }
            }
        }
        // Finest level: no truncation (nothing finer to subtract).
        const SizeType num_cps_u_finest =
            mLevels[num_levels - 1].KnotsU.size() - polynomial_degree_u + 1;
        const SizeType num_cps_v_finest =
            mLevels[num_levels - 1].KnotsV.size() - polynomial_degree_v + 1;
        mTruncationData[num_levels - 1].assign(
            num_cps_u_finest * num_cps_v_finest, TruncationEntryContainer{});
    }

    /**
     * @brief Top-down propagation that assigns active/inactive flags to every CP.
     *
     * Initialization:
     *   mActiveFunctions[0]   = all true   (level 0 covers the full domain)
     *   mActiveFunctions[l>0] = all false  (activated only by coarser propagation)
     *
     * Main loop  l = 0 .. L-2:
     *   For each active CP (i_u, i_v) at level l:
     *     1. Compute its support [support_min_u, support_max_u] x [support_min_v, support_max_v].
     *     2. Check whether EVERY non-zero cell in the support has its midpoint inside
     *        a refinement domain at level >= l+1.  If not → coarse CP stays active.
     *     3. If all cells covered → mark coarse CP INACTIVE and activate true children.
     *        A fine CP j is a child iff supp(B_j^{l+1}) ⊆ supp(B_i^l).
     *
     * Called automatically by AddRefinementDomain.
     */
    void ComputeActiveFunctions()
    {
        const SizeType polynomial_degree_u = mLevels[0].DegreeU;
        const SizeType polynomial_degree_v = mLevels[0].DegreeV;

        // Level 0: fully active. Finer levels: start inactive, activated by propagation.
        std::fill(mActiveFunctions[0].begin(), mActiveFunctions[0].end(), true);
        for (SizeType l = 1; l < mLevels.size(); ++l)
            std::fill(mActiveFunctions[l].begin(), mActiveFunctions[l].end(), false);

        for (SizeType l = 0; l + 1 < mLevels.size(); ++l) {
            const THBLevel& coarse_level = mLevels[l];
            const THBLevel& fine_level   = mLevels[l + 1];
            const SizeType num_cps_u_coarse = coarse_level.KnotsU.size() - polynomial_degree_u + 1;
            const SizeType num_cps_v_coarse = coarse_level.KnotsV.size() - polynomial_degree_v + 1;
            const SizeType num_cps_u_fine   = fine_level.KnotsU.size()   - polynomial_degree_u + 1;
            const SizeType num_cps_v_fine   = fine_level.KnotsV.size()   - polynomial_degree_v + 1;

            for (SizeType i_v = 0; i_v < num_cps_v_coarse; ++i_v) {
                for (SizeType i_u = 0; i_u < num_cps_u_coarse; ++i_u) {
                    if (!mActiveFunctions[l][i_v * num_cps_u_coarse + i_u]) continue;

                    // Support of coarse B_{i_u, i_v}^l
                    const double support_min_u =
                        (i_u == 0) ? coarse_level.KnotsU[0] : coarse_level.KnotsU[i_u - 1];
                    const double support_max_u =
                        coarse_level.KnotsU[std::min(i_u + polynomial_degree_u,
                                                     coarse_level.KnotsU.size() - 1)];
                    const double support_min_v =
                        (i_v == 0) ? coarse_level.KnotsV[0] : coarse_level.KnotsV[i_v - 1];
                    const double support_max_v =
                        coarse_level.KnotsV[std::min(i_v + polynomial_degree_v,
                                                     coarse_level.KnotsV.size() - 1)];

                    // Is every cell of this support inside Ω^{≥l+1}?
                    bool all_cells_covered = true;
                    for (SizeType k_u = 0;
                         k_u + 1 < coarse_level.KnotsU.size() && all_cells_covered; ++k_u) {
                        if (coarse_level.KnotsU[k_u + 1] - coarse_level.KnotsU[k_u] < 1e-10) continue;
                        if (coarse_level.KnotsU[k_u]     < support_min_u - 1e-10) continue;
                        if (coarse_level.KnotsU[k_u + 1] > support_max_u + 1e-10) continue;
                        for (SizeType k_v = 0;
                             k_v + 1 < coarse_level.KnotsV.size() && all_cells_covered; ++k_v) {
                            if (coarse_level.KnotsV[k_v + 1] - coarse_level.KnotsV[k_v] < 1e-10) continue;
                            if (coarse_level.KnotsV[k_v]     < support_min_v - 1e-10) continue;
                            if (coarse_level.KnotsV[k_v + 1] > support_max_v + 1e-10) continue;
                            const double cell_mid_u =
                                0.5 * (coarse_level.KnotsU[k_u] + coarse_level.KnotsU[k_u + 1]);
                            const double cell_mid_v =
                                0.5 * (coarse_level.KnotsV[k_v] + coarse_level.KnotsV[k_v + 1]);
                            if (!IsInsideRefinedRegion(cell_mid_u, cell_mid_v, l + 1))
                                all_cells_covered = false;
                        }
                    }
                    if (!all_cells_covered) continue;

                    // Mark coarse CP INACTIVE and activate its true children.
                    mActiveFunctions[l][i_v * num_cps_u_coarse + i_u] = false;

                    for (SizeType j_u = 0; j_u < num_cps_u_fine; ++j_u) {
                        const double fine_supp_min_u =
                            (j_u == 0) ? fine_level.KnotsU[0] : fine_level.KnotsU[j_u - 1];
                        const double fine_supp_max_u =
                            fine_level.KnotsU[std::min(j_u + polynomial_degree_u,
                                                       fine_level.KnotsU.size() - 1)];
                        if (fine_supp_min_u < support_min_u - 1e-10) continue;
                        if (fine_supp_max_u > support_max_u + 1e-10) continue;

                        for (SizeType j_v = 0; j_v < num_cps_v_fine; ++j_v) {
                            const double fine_supp_min_v =
                                (j_v == 0) ? fine_level.KnotsV[0] : fine_level.KnotsV[j_v - 1];
                            const double fine_supp_max_v =
                                fine_level.KnotsV[std::min(j_v + polynomial_degree_v,
                                                           fine_level.KnotsV.size() - 1)];
                            if (fine_supp_min_v < support_min_v - 1e-10) continue;
                            if (fine_supp_max_v > support_max_v + 1e-10) continue;

                            mActiveFunctions[l + 1][j_v * num_cps_u_fine + j_u] = true;
                        }
                    }
                }
            }
        }
    }

    /// Returns true if (u, v) is inside any refinement domain registered for the given level.
    bool IsInsideDomain(double u, double v, SizeType Level) const
    {
        for (const auto& dom : mRefinementDomains) {
            if (dom.Level == Level &&
                u >= dom.MinU && u <= dom.MaxU &&
                v >= dom.MinV && v <= dom.MaxV)
                return true;
        }
        return false;
    }

    /// Returns the number of control points in u for a given level.
    SizeType NumberOfControlPointsU(SizeType Level) const
    {
        return mLevels[Level].KnotsU.size() - mLevels[0].DegreeU + 1;
    }

    /// Returns the number of control points in v for a given level.
    SizeType NumberOfControlPointsV(SizeType Level) const
    {
        return mLevels[Level].KnotsV.size() - mLevels[0].DegreeV + 1;
    }

    /// Returns the index offset into this->Points() where level Level's control points begin.
    SizeType ControlPointOffset(SizeType Level) const
    {
        SizeType offset = 0;
        for (SizeType l = 0; l < Level; ++l)
            offset += NumberOfControlPointsU(l) * NumberOfControlPointsV(l);
        return offset;
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    THBSurfaceGeometry() : BaseType(PointsArrayType()) {}

    ///@}

}; // class THBSurfaceGeometry

} // namespace Kratos
