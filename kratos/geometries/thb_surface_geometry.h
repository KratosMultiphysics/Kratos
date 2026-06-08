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
#include "geometries/local_refined_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"
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
    }

    ~THBSurfaceGeometry() override = default;

    ///@}
    ///@name Refinement
    ///@{

    /**
     * @brief Adds a refinement level and computes its control points via knot insertion.
     *
     * The new level's control points are derived from the previous level's control points
     * using the tensor-product B-spline refinement operator (Oslo algorithm / iterative
     * single-knot insertion). The computed nodes are appended to this->Points().
     *
     * Call once per level, in order (level 1, then 2, ...).
     * Polynomial degree is inherited from level 0.
     * rKnotsU/rKnotsV must be a refinement of the previous level's knot vectors
     * (i.e., the previous knots must be a subset of the new ones).
     */
    void AddLevel(
        const Vector& rKnotsU,
        const Vector& rKnotsV,
        const Vector& rWeights = Vector())
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::AddLevel: no base level defined." << std::endl;

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

        // Assign IDs beyond any existing node in this geometry
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
                this->Points().push_back(
                    Kratos::make_intrusive<NodeType>(new_id++, x, y, z));
            }
        }

        mLevels.push_back(THBLevel{prev.DegreeU, prev.DegreeV, rKnotsU, rKnotsV, rWeights});
    }

    /**
     * @brief Adds NLevels refinement levels via uniform bisection.
     *
     * Each new level is obtained by inserting the midpoint of every unique knot
     * span of the previous level (halving all spans). Equivalent to calling
     * AddLevel(knots_u, knots_v) NLevels times with auto-generated knot vectors.
     */
    void AddLevel(const SizeType NLevels)
    {
        KRATOS_ERROR_IF(mLevels.empty())
            << "THBSurfaceGeometry::AddLevel: no base level defined." << std::endl;

        for (SizeType n = 0; n < NLevels; ++n) {
            const THBLevel& prev = mLevels.back();
            AddLevel(BisectKnots(prev.KnotsU), BisectKnots(prev.KnotsV));
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
    }

    ///@}
    ///@name Getters
    ///@{

    const THBLevelContainer& Levels() const { return mLevels; }

    const RefinementDomainContainer& RefinementDomains() const { return mRefinementDomains; }

    SizeType NumberOfLevels() const { return mLevels.size(); }

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
     * @brief Builds quadrature-point geometries using B-spline shape functions
     *        at the locally active level for each integration point.
     *
     * NOTE: Truncation is ignored — shape functions from coarser levels
     * are fully deactivated inside refinement domains (not truncated).
     * This violates partition of unity and is a temporary simplification.
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

        for (IndexType i = 0; i < rIntegrationPoints.size(); ++i) {
            const double u = rIntegrationPoints[i][0];
            const double v = rIntegrationPoints[i][1];

            const SizeType l = ActiveLevelAtPoint(u, v);
            const THBLevel& level = mLevels[l];
            // Internal knot format: size = n + p - 1  →  n = size - p + 1
            const SizeType nU_l = level.KnotsU.size() - p + 1;
            const SizeType nV_l = level.KnotsV.size() - q + 1;
            const SizeType offset = ControlPointOffset(l);

            NurbsSurfaceShapeFunction sf(p, q, NumberOfShapeFunctionDerivatives);
            /// THBSurfaceShapeFunction sf(p, q, NumberOfShapeFunctionDerivatives);  ///todo
            if (level.Weights.size() > 0) {
                sf.ComputeNurbsShapeFunctionValues(level.KnotsU, level.KnotsV, level.Weights, u, v);
            } else {
                sf.ComputeBSplineShapeFunctionValues(level.KnotsU, level.KnotsV, u, v);
            }

            const SizeType num_nonzero = sf.NumberOfNonzeroControlPoints();
            auto local_cp_indices = sf.ControlPointIndices(nU_l, nV_l);

            PointsArrayType nonzero_control_points(num_nonzero);
            for (IndexType j = 0; j < num_nonzero; ++j)
                nonzero_control_points(j) = this->pGetPoint(offset + local_cp_indices[j]);

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

    ///@}
    ///@name Private Helpers
    ///@{

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
