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

    using THBLevelContainer         = std::vector<THBLevel>;
    using RefinementDomainContainer = std::vector<THBRefinementDomain>;

    /// Placeholder — truncation metadata structure is not yet defined.
    using TruncationRuleContainer   = std::vector<int>;

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
        const SizeType nU = rKnotsU.size() - PolynomialDegreeU + 1;
        const SizeType nV = rKnotsV.size() - PolynomialDegreeV + 1;
        mActiveFunctions.push_back(std::vector<bool>(nU * nV, true));
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
        const SizeType nU = rKnotsU.size() - PolynomialDegreeU + 1;
        const SizeType nV = rKnotsV.size() - PolynomialDegreeV + 1;
        mActiveFunctions.push_back(std::vector<bool>(nU * nV, true));
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
        for (SizeType n = 0; n < NLevels; ++n) {
            const THBLevel& prev = mLevels.back();
            AddLevel(BisectKnots(prev.KnotsU), BisectKnots(prev.KnotsV));
        }
    }

    /// Adds NLevels via uniform bisection; new nodes are registered in rModelPart.
    void AddLevel(const SizeType NLevels, ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::AddLevel: no base level defined." << std::endl;
        for (SizeType n = 0; n < NLevels; ++n) {
            const THBLevel& prev = mLevels.back();
            AddLevel(BisectKnots(prev.KnotsU), BisectKnots(prev.KnotsV), rModelPart);
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

        const SizeType p = mLevels[0].DegreeU;
        const SizeType q = mLevels[0].DegreeV;
        const SizeType nLevels = mLevels.size();

        mActiveOffset.resize(nLevels);
        mActiveLocalIndex.resize(nLevels);

        PointsArrayType packed_points;
        SizeType running_offset = 0;

        for (SizeType l = 0; l < nLevels; ++l) {
            const SizeType nU    = mLevels[l].KnotsU.size() - p + 1;
            const SizeType nV    = mLevels[l].KnotsV.size() - q + 1;
            const SizeType total = nU * nV;
            const SizeType old_offset = ControlPointOffset(l);  // cumulative before compaction

            mActiveOffset[l] = running_offset;
            mActiveLocalIndex[l].assign(total, -1);

            int packed = 0;
            for (SizeType flat = 0; flat < total; ++flat) {
                auto p_node = this->pGetPoint(old_offset + flat);
                if (mActiveFunctions[l][flat]) {
                    mActiveLocalIndex[l][flat] = packed++;
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

    /// Returns the finest level whose refinement domain contains (u, v), or 0 if none.
    SizeType ActiveLevelAtPoint(double u, double v) const
    {
        SizeType active = 0;
        for (SizeType l = 1; l < mLevels.size(); ++l) {
            if (IsInsideDomain(u, v, l))
                active = l;
        }
        return active;
    }

    /// Maps (level, flat-index-within-level) to the index in the packed Points() array.
    /// Works both before and after EliminateInactiveFunctions.
    SizeType PackedControlPointIndex(SizeType Level, SizeType FlatIndex) const
    {
        if (mIsEliminated) {
            KRATOS_DEBUG_ERROR_IF(mActiveLocalIndex[Level][FlatIndex] < 0)
                << "THBSurfaceGeometry::PackedControlPointIndex: "
                   "inactive function requested — logic error." << std::endl;
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

    /// @todo Implement region-aware span decomposition.
    ///       Blocked on CreateIntegrationPoints; a flat global union is not correct for THB.
    void SpansLocalSpace(
        std::vector<double>& rSpans,
        IndexType DirectionIndex) const override
    {
        KRATOS_ERROR << "THBSurfaceGeometry::SpansLocalSpace: "
            "not yet implemented — requires THB integration-cell algorithm." << std::endl;
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
        const SizeType p = PolynomialDegree(0);
        const SizeType q = PolynomialDegree(1);
        const SizeType points_u = p + 1;
        const SizeType points_v = q + 1;

        rIntegrationPoints.clear();

        IntegrationPointsArrayType span_ips(points_u * points_v);

        for (SizeType l = 0; l < mLevels.size(); ++l) {
            auto spans_u = KnotSpanIntervals(mLevels[l].KnotsU);
            auto spans_v = KnotSpanIntervals(mLevels[l].KnotsV);

            for (const auto& su : spans_u) {
                for (const auto& sv : spans_v) {
                    const double mid_u = 0.5 * (su.GetT0() + su.GetT1());
                    const double mid_v = 0.5 * (sv.GetT0() + sv.GetT1());

                    // Active at level l: inside Omega^l (or level 0), AND not refined further
                    const bool in_this_level = (l == 0) || IsInsideDomain(mid_u, mid_v, l);
                    const bool in_next_level = (l + 1 < mLevels.size()) && IsInsideDomain(mid_u, mid_v, l + 1);

                    if (!in_this_level || in_next_level)
                        continue;

                    auto it = span_ips.begin();
                    IntegrationPointUtilities::IntegrationPoints2D(
                        it, points_u, points_v,
                        su.GetT0(), su.GetT1(),
                        sv.GetT0(), sv.GetT1());

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
     *
     * @note Truncation is not yet applied — this computes HB shape functions.
     *       Partition of unity is not guaranteed at coarse/fine boundaries until
     *       ApplyTruncation is implemented in THBSurfaceShapeFunction.
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
        const SizeType q = PolynomialDegree(1);

        THBSurfaceShapeFunction sf(p, q, NumberOfShapeFunctionDerivatives);

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i) {
            const double u = rIntegrationPoints[i][0];
            const double v = rIntegrationPoints[i][1];

            sf.ComputeShapeFunctionValues(*this, u, v);

            const SizeType num_nonzero = sf.NumberOfNonzeroControlPoints();
            const auto& cp_indices = sf.ControlPointIndices();

            PointsArrayType nonzero_control_points(num_nonzero);
            for (IndexType j = 0; j < num_nonzero; ++j)
                nonzero_control_points(j) = this->pGetPoint(cp_indices[j]);

            Matrix N(1, num_nonzero);
            for (IndexType j = 0; j < num_nonzero; ++j)
                N(0, j) = sf(j, 0);

            DenseVector<Matrix> sf_derivatives(NumberOfShapeFunctionDerivatives - 1);
            for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n)
                sf_derivatives[n].resize(num_nonzero, n + 2);

            if (NumberOfShapeFunctionDerivatives > 0) {
                IndexType deriv_index = 1;
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n) {
                    for (IndexType k = 0; k < n + 2; ++k)
                        for (IndexType j = 0; j < num_nonzero; ++j)
                            sf_derivatives[n](j, k) = sf(j, deriv_index + k);
                    deriv_index += n + 2;
                }
            }

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[i], N, sf_derivatives);

            rResultGeometries(i) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePoint(
                this->WorkingSpaceDimension(), 2, data_container, nonzero_control_points, this);
        }
    }

    /// @todo Implement using THBSurfaceShapeFunction.
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates) const override
    {
        KRATOS_ERROR << "THBSurfaceGeometry::GlobalCoordinates: not yet implemented." << std::endl;
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
    TruncationRuleContainer   mTruncationRules; ///< @todo define once truncation metadata is specified
    std::vector<std::vector<bool>> mActiveFunctions;

    bool                           mIsEliminated = false;
    std::vector<SizeType>          mActiveOffset;      // packed Points() start index for level l
    std::vector<std::vector<int>>  mActiveLocalIndex;  // [l][flat_idx] → packed pos within l, -1 if inactive

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

        const SizeType l = mLevels.size() - 1;
        const THBLevel& prev = mLevels[l];
        const SizeType p = prev.DegreeU;
        const SizeType q = prev.DegreeV;
        const SizeType n_u_old = prev.KnotsU.size() - p + 1;
        const SizeType n_v_old = prev.KnotsV.size() - q + 1;

        const Matrix M_U = ComputeRefinementMatrix1D(prev.KnotsU, rKnotsU, p);
        const Matrix M_V = ComputeRefinementMatrix1D(prev.KnotsV, rKnotsV, q);
        const SizeType n_u_new = M_U.size1();
        const SizeType n_v_new = M_V.size1();

        const SizeType offset = ControlPointOffset(l);

        IndexType new_id = 1;
        for (SizeType i = 0; i < this->PointsNumber(); ++i)
            new_id = std::max(new_id, this->GetPoint(i).Id() + IndexType(1));

        // P'_{j_v, j_u} = sum_{i_v, i_u} M_V[j_v, i_v] * M_U[j_u, i_u] * P_{i_v, i_u}
        for (SizeType j_v = 0; j_v < n_v_new; ++j_v) {
            for (SizeType j_u = 0; j_u < n_u_new; ++j_u) {
                double x = 0.0, y = 0.0, z = 0.0;
                for (SizeType i_v = 0; i_v < n_v_old; ++i_v) {
                    const double cv = M_V(j_v, i_v);
                    if (std::abs(cv) < 1e-15) continue;
                    for (SizeType i_u = 0; i_u < n_u_old; ++i_u) {
                        const double coeff = cv * M_U(j_u, i_u);
                        if (std::abs(coeff) < 1e-15) continue;
                        const auto& pt = this->GetPoint(offset + i_v * n_u_old + i_u);
                        x += coeff * pt.X();
                        y += coeff * pt.Y();
                        z += coeff * pt.Z();
                    }
                }
                this->Points().push_back(CreateNode(new_id++, x, y, z));
            }
        }

        mLevels.push_back(THBLevel{prev.DegreeU, prev.DegreeV, rKnotsU, rKnotsV, rWeights});
        mActiveFunctions.push_back(std::vector<bool>(n_u_new * n_v_new, true));
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
     * @brief Computes the 1-D B-spline refinement matrix M (size n_fine × n_coarse).
     *
     * Given internal-format knot vectors for a coarse and a fine level (fine must be
     * a refinement of coarse), returns the matrix M such that:
     *   P_fine[j] = sum_i  M[j, i] * P_coarse[i]
     *
     * Algorithm: iterative single-knot insertion.
     * Works in external (clamped) format internally; inputs/outputs use Kratos
     * internal format (n + p - 1 knots).
     */
    static Matrix ComputeRefinementMatrix1D(
        const Vector& rCoarseKnots,
        const Vector& rFineKnots,
        const SizeType p)
    {
        const SizeType n_c = rCoarseKnots.size() - p + 1;

        // Convert internal → external format: prepend/append first/last knot
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

        // New knots = ext_f multiset minus ext_c multiset (both sorted)
        std::vector<double> new_knots;
        {
            SizeType ci = 0;
            for (SizeType fi = 0; fi < ext_f.size(); ++fi) {
                if (ci < ext_c.size() && std::abs(ext_f[fi] - ext_c[ci]) < 1e-10) {
                    ++ci;
                } else {
                    new_knots.push_back(ext_f[fi]);
                }
            }
        }

        // Start with the n_c × n_c identity; grow to n_f × n_c by inserting knots one by one.
        Matrix M(n_c, n_c, 0.0);
        for (SizeType i = 0; i < n_c; ++i) M(i, i) = 1.0;

        std::vector<double> cur = ext_c;
        SizeType n_cur = n_c;

        for (const double t : new_knots) {
            // Span k: last index with cur[k] <= t
            SizeType k = 0;
            for (SizeType i = 0; i + 1 < cur.size(); ++i)
                if (cur[i] <= t + 1e-10) k = i;

            // Single-insertion matrix T: (n_cur+1) × n_cur
            Matrix T(n_cur + 1, n_cur, 0.0);
            const int ik = static_cast<int>(k);
            const int ip = static_cast<int>(p);
            for (int i = 0; i <= static_cast<int>(n_cur); ++i) {
                if (i <= ik - ip) {
                    T(i, i) = 1.0;
                } else if (i > ik) {
                    T(i, i - 1) = 1.0;
                } else {
                    // Blend: alpha * P_i + (1-alpha) * P_{i-1}
                    const double denom = cur[i + p] - cur[i];
                    const double alpha = (std::abs(denom) < 1e-15) ? 1.0 : (t - cur[i]) / denom;
                    if (i < static_cast<int>(n_cur)) T(i, i) = alpha;
                    if (i > 0) T(i, i - 1) = 1.0 - alpha;
                }
            }

            M = prod(T, M);
            cur.insert(cur.begin() + k + 1, t);
            ++n_cur;
        }

        return M;  // size n_f × n_c
    }

    /// Returns unique non-zero knot span intervals from an internal-format knot vector.
    /// Kratos convention: n + p - 1 knots (one repeated knot dropped from each end
    /// of the full clamped vector). Iterates the full stored vector for transitions.
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
    /// Used to check whether a knot cell is covered by a finer hierarchy level.
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
     * @brief Top-down propagation that assigns active/inactive flags to every CP.
     *
     * Initialization:
     *   mActiveFunctions[0]   = all true   (level 0 covers the full domain)
     *   mActiveFunctions[l>0] = all false  (activated only by coarser propagation)
     *
     * Main loop  l = 0 .. L-2:
     *   For each active CP (i_u, i_v) at level l:
     *     1. Compute its support [cu_min, cu_max] x [cv_min, cv_max].
     *     2. Check whether EVERY non-zero cell in the support has its midpoint
     *        inside a refinement domain at level >= l+1.
     *        If not → the coarse CP stays active; skip to the next CP.
     *     3. If all cells are covered → mark the coarse CP INACTIVE and
     *        activate its children at level l+1.
     *        A fine CP j is a child of coarse CP i iff
     *          supp(B_j^{l+1}) ⊆ supp(B_i^l)
     *        i.e. fu_min >= cu_min AND fu_max <= cu_max (same in v).
     *        Fine CPs whose support straddles the boundary of the coarse support
     *        are NOT children and are left inactive; the still-active coarse CPs
     *        at the boundary already represent that region.
     *
     * Cascade: if a newly activated child's support is fully inside Ω^{l+2}, the
     * next loop iteration (l+1) finds that child's support fully covered, deactivates it,
     * and activates its own children — no special handling needed.
     *
     * Called automatically by AddRefinementDomain.
     */
    void ComputeActiveFunctions()
    {
        const SizeType p = mLevels[0].DegreeU;
        const SizeType q = mLevels[0].DegreeV;

        // Level 0: fully active. Finer levels: start inactive, activated by propagation.
        std::fill(mActiveFunctions[0].begin(), mActiveFunctions[0].end(), true);
        for (SizeType l = 1; l < mLevels.size(); ++l)
            std::fill(mActiveFunctions[l].begin(), mActiveFunctions[l].end(), false);

        for (SizeType l = 0; l + 1 < mLevels.size(); ++l) {
            const THBLevel& lev_c = mLevels[l];
            const THBLevel& lev_f = mLevels[l + 1];
            const SizeType  nU_c  = lev_c.KnotsU.size() - p + 1;
            const SizeType  nV_c  = lev_c.KnotsV.size() - q + 1;
            const SizeType  nU_f  = lev_f.KnotsU.size() - p + 1;
            const SizeType  nV_f  = lev_f.KnotsV.size() - q + 1;

            for (SizeType i_v = 0; i_v < nV_c; ++i_v) {
                for (SizeType i_u = 0; i_u < nU_c; ++i_u) {
                    if (!mActiveFunctions[l][i_v * nU_c + i_u]) continue;

                    // Support of coarse B_{i_u, i_v}^l
                    const double cu_min = (i_u == 0) ? lev_c.KnotsU[0] : lev_c.KnotsU[i_u - 1];
                    const double cu_max = lev_c.KnotsU[std::min(i_u + p, lev_c.KnotsU.size() - 1)];
                    const double cv_min = (i_v == 0) ? lev_c.KnotsV[0] : lev_c.KnotsV[i_v - 1];
                    const double cv_max = lev_c.KnotsV[std::min(i_v + q, lev_c.KnotsV.size() - 1)];

                    // Is every cell of this support inside Ω^{≥l+1}?
                    bool all_cells_covered = true;
                    for (SizeType k_u = 0; k_u + 1 < lev_c.KnotsU.size() && all_cells_covered; ++k_u) {
                        if (lev_c.KnotsU[k_u + 1] - lev_c.KnotsU[k_u] < 1e-10) continue;
                        if (lev_c.KnotsU[k_u] < cu_min - 1e-10 || lev_c.KnotsU[k_u + 1] > cu_max + 1e-10) continue;
                        for (SizeType k_v = 0; k_v + 1 < lev_c.KnotsV.size() && all_cells_covered; ++k_v) {
                            if (lev_c.KnotsV[k_v + 1] - lev_c.KnotsV[k_v] < 1e-10) continue;
                            if (lev_c.KnotsV[k_v] < cv_min - 1e-10 || lev_c.KnotsV[k_v + 1] > cv_max + 1e-10) continue;
                            const double mid_u = 0.5 * (lev_c.KnotsU[k_u] + lev_c.KnotsU[k_u + 1]);
                            const double mid_v = 0.5 * (lev_c.KnotsV[k_v] + lev_c.KnotsV[k_v + 1]);
                            if (!IsInsideRefinedRegion(mid_u, mid_v, l + 1))
                                all_cells_covered = false;
                        }
                    }
                    if (!all_cells_covered) continue; // coarse CP stays active, no propagation

                    // Mark coarse CP INACTIVE and activate its true children.
                    mActiveFunctions[l][i_v * nU_c + i_u] = false;

                    for (SizeType j_u = 0; j_u < nU_f; ++j_u) {
                        const double fu_min = (j_u == 0) ? lev_f.KnotsU[0] : lev_f.KnotsU[j_u - 1];
                        const double fu_max = lev_f.KnotsU[std::min(j_u + p, lev_f.KnotsU.size() - 1)];
                        if (fu_min < cu_min - 1e-10 || fu_max > cu_max + 1e-10) continue;

                        for (SizeType j_v = 0; j_v < nV_f; ++j_v) {
                            const double fv_min = (j_v == 0) ? lev_f.KnotsV[0] : lev_f.KnotsV[j_v - 1];
                            const double fv_max = lev_f.KnotsV[std::min(j_v + q, lev_f.KnotsV.size() - 1)];
                            if (fv_min < cv_min - 1e-10 || fv_max > cv_max + 1e-10) continue;

                            mActiveFunctions[l + 1][j_v * nU_f + j_u] = true;
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
    /// Internal knot format: size = n + p - 1  →  n = size - p + 1
    SizeType NumberOfControlPointsU(SizeType Level) const
    {
        const SizeType p = mLevels[0].DegreeU;
        return mLevels[Level].KnotsU.size() - p + 1;
    }

    /// Returns the number of control points in v for a given level.
    SizeType NumberOfControlPointsV(SizeType Level) const
    {
        const SizeType q = mLevels[0].DegreeV;
        return mLevels[Level].KnotsV.size() - q + 1;
    }

    /// Returns the index offset into this->Points() where level Level's control points begin.
    /// Assumes control points are stored level by level (level 0 first, then 1, ...).
    SizeType ControlPointOffset(SizeType Level) const
    {
        SizeType offset = 0;
        for (SizeType l = 0; l < Level; ++l)
            offset += NumberOfControlPointsU(l) * NumberOfControlPointsV(l);
        return offset;
    }

    ///@}

}; // class THBSurfaceGeometry

} // namespace Kratos
