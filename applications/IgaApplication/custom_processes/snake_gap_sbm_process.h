//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//

#pragma once

// System includes
#include <utility>
#include <unordered_map>
#include <vector>

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"
#include "processes/process.h"
#include "geometries/nurbs_curve_geometry.h"
#include "snake_sbm_process.h"
#include "custom_utilities/create_breps_sbm_utilities.h"
#include "includes/global_pointer_variables.h"

#include "spaces/ublas_space.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class SnakeGapSbmProcess
 * @brief Process class for implementing Snake-based Surrogate Boundary Method (SBM).
 * This class provides various functions to create and manipulate surrogate boundaries
 * for Iga models in Kratos.
 */
class KRATOS_API(IGA_APPLICATION) SnakeGapSbmProcess
    : public SnakeSbmProcess
{

public:

    using NodeType = Node;
    using ContainerNodeType = PointerVector<Node>;
    using ContainerEmbeddedNodeType = PointerVector<Point>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;
    using ElementsContainerType = ModelPart::ElementsContainerType;

    using GeometryType = Geometry<NodeType>;
    using GeometriesArrayType = GeometryType::GeometriesArrayType;
    using IntegrationPointsArrayType = GeometryType::IntegrationPointsArrayType;
    using CoordinatesArrayType = GeometryType::CoordinatesArrayType;
    using PropertiesPointerType = Properties::Pointer;
    using BrepCurveType = BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType>;
    
    using NurbsCurveGeometryType = NurbsCurveGeometry<3, PointerVector<Node>>;
    using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<NodeType>>;
    using NodePointerVector = GlobalPointersVector<NodeType>;

    using SparseSpaceType  = UblasSpace<double, CompressedMatrix, Vector>;
    using SparseMatrixType = SparseSpaceType::MatrixType;

    struct KnotSpanIdsCSR
    {
        // Metadata (as in your previous struct)
        std::size_t NumberOfSpansX = 0;
        std::size_t NumberOfSpansY = 0;
        double MinU = 0.0, MaxU = 0.0;
        double MinV = 0.0, MaxV = 0.0;
        double SpanSizeX = 0.0, SpanSizeY = 0.0;

        // Sparse topology and counts per nonzero
        SparseMatrixType Occupancy;          // CSR, (i,j) value = count of node ids in that cell

        // Compact payload: for each nonzero k, node ids are stored in [pool[nnz_off[k]], ..., length nnz_len[k]]
        std::vector<std::size_t> nnz_off;       // offsets into pool, size = nnz
        std::vector<std::size_t> nnz_len;       // lengths per nnz,  size = nnz
        std::vector<IndexType> pool;         // flat contiguous storage of node ids
    };

    struct IdsView {
        const IndexType* data = nullptr;
        std::size_t size = 0;
    };

    // Return the CSR "k" index for entry (i,j), or static_cast<std::size_t>(-1) if not found.
    static std::size_t FindNnzIndex(const SparseMatrixType& A, std::size_t i, std::size_t j)
    {
        if (i >= A.size1() || j >= A.size2()) return static_cast<std::size_t>(-1);

        const auto& rp = A.index1_data();   // row_ptr, size = size1+1
        const auto& ci = A.index2_data();   // col_ind, size = nnz

        if (rp.size() < A.size1() + 1) return static_cast<std::size_t>(-1);
        std::size_t b = rp[i];
        std::size_t e = rp[i + 1];

        // guard rails against corrupted/unfinalized structure
        const std::size_t nnz = static_cast<std::size_t>(ci.size());
        if (b > e) return static_cast<std::size_t>(-1);
        if (b > nnz) b = nnz;
        if (e > nnz) e = nnz;

        const auto first = ci.begin() + b;
        const auto last  = ci.begin() + e;
        const auto it = std::lower_bound(first, last, j);
        if (it != last && *it == j) return static_cast<std::size_t>(it - ci.begin());
        return static_cast<std::size_t>(-1);
    }


    // Pack temporary per-cell vectors into the compact pool and set offsets/lengths.
    // Also sorts and unique-erases each list.
    static void CommitPayload(std::vector<std::vector<IndexType>>& tmp,
                                    KnotSpanIdsCSR& out)
    {
        out.nnz_off.resize(tmp.size());
        out.nnz_len.resize(tmp.size());

        for (std::size_t k = 0; k < tmp.size(); ++k) {
            auto& v = tmp[k];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());

            out.nnz_off[k] = static_cast<std::size_t>(out.pool.size());
            out.nnz_len[k] = static_cast<std::size_t>(v.size());
            out.pool.insert(out.pool.end(), v.begin(), v.end());
        }
    }

    // Access node ids for a given nonzero k
    static IdsView CellIdsByK(const KnotSpanIdsCSR& S, const std::size_t k)
    {
        return IdsView{ S.pool.data() + S.nnz_off[k], S.nnz_len[k] };
    }

    // Access node ids for a given cell (i,j). Returns empty view if cell empty.
    static IdsView CellIds(const KnotSpanIdsCSR& S, const std::size_t i, const std::size_t j)
    {
        const std::size_t k = FindNnzIndex(S.Occupancy, i, j);
        if (k == static_cast<std::size_t>(-1)) return {};
        return CellIdsByK(S, k);
    }

    ///@name Type Definitions
    ///@{

    /// Pointer definition of SnakeGapSbmProcess
    KRATOS_CLASS_POINTER_DEFINITION(SnakeGapSbmProcess);

    ///@name Life Cycle
    ///@{

    /// Constructor
    SnakeGapSbmProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~SnakeGapSbmProcess() = default;

    // Get the default parameters
    const Parameters GetDefaultParameters() const override;
    const Parameters GetValidParameters() const;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {   
        CreateSbmExtendedGeometries();
    };
    
    void ExecuteInitialize() override
    {
        SnakeSbmProcess::CreateTheSnakeCoordinates();
    };

    void ExecuteInitializeSolutionStep() override
    {};

    void ExecuteBeforeSolutionLoop() override
    {};

    KnotSpanIdsCSR CreateSkinNodesPerKnotSpanMatrix(
        const ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart) const;

    KnotSpanIdsCSR CreateSkinConditionsPerKnotSpanMatrix(
        const ModelPart& rSkinSubModelPart,
        const SnakeGapSbmProcess::KnotSpanIdsCSR& rReferenceMatrix) const;


struct BinSearchParameters
{
    DynamicBins&                                   rTestBins;
    std::size_t                                       NumberOfResults;
    ModelPart::NodesContainerType::ContainerType   Results;
    std::vector<double>                            ListOfDistances;
    double                                         SearchRadius;

    BinSearchParameters(
        DynamicBins&                               rBins,
        std::size_t                                   NumberOfResultsInput,
        ModelPart::NodesContainerType::ContainerType ResultsInput,
        std::vector<double>                        DistancesInput,
        double                                     SearchRadiusInput)
        : rTestBins(rBins)
        , NumberOfResults(NumberOfResultsInput)
        , Results(std::move(ResultsInput))
        , ListOfDistances(std::move(DistancesInput))
        , SearchRadius(SearchRadiusInput)
    {}

    void reset()
    {
        Results.clear();
        ListOfDistances.clear();
    }
};


struct IntegrationParameters
{
    std::size_t  NumberOfShapeFunctionsDerivatives;
    IntegrationInfo CurveIntegrationInfo;
    const Vector&           rKnotSpanSizes;
    const KnotSpanIdsCSR* pSkinNodesPerSpan = nullptr;
    const KnotSpanIdsCSR* pSkinConditionsPerSpan = nullptr;
    IntegrationParameters(
        std::size_t                    NumberOfDerivatives,
        const IntegrationInfo&      rInfo,
        const Vector&               rKnotSpans)
        : NumberOfShapeFunctionsDerivatives(NumberOfDerivatives)
        , CurveIntegrationInfo(rInfo)
        , rKnotSpanSizes(rKnotSpans)
    {}
};
    
private:

    ModelPart* mpGapElementsSubModelPart = nullptr; 
    ModelPart* mpGapInterfaceSubModelPart = nullptr; 
    std::string mGapElementName;
    std::string mGapInterfaceConditionName;
    std::size_t mInternalDivisions;
    std::size_t mGapApproximationOrder;
    std::string mGapSbmType;

    /**
     * @brief Creates the gap SBM geometries for the configured model part.
     */
    void CreateSbmExtendedGeometries();

    /**
     * @brief Builds the gap SBM geometries for a given skin/surrogate pair.
     * @tparam TIsInnerLoop True when processing the inner loop orientation.
     * @param rSkinSubModelPart Reference to the target skin sub model part.
     * @param rSurrogateSubModelPart Reference to the surrogate sub model part.
     */
    template <bool TIsInnerLoop>
    void CreateSbmExtendedGeometries(
        const ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart);

    /**
     * @brief Creates gap and skin quadrature point data for the supplied entities.
     * @tparam TIsInnerLoop True when processing inner loop quadrature information.
     * @param rIntegrationParameters Reference to the integration parameters container.
     * @param rBinSearchParameters Reference to the bin search configuration.
     * @param pNurbsSurface Pointer to the NURBS surface used for projections.
     * @param pSurrogateNode1 Pointer to the first surrogate node of the pair.
     * @param pSurrogateNode2 Pointer to the second surrogate node of the pair.
     * @param rSurrogateBrepMiddleGeometry Reference to the intermediate surrogate geometry.
     * @param rIgaModelPart Reference to the destination IGA model part.
     * @param rSkinSubModelPart Reference to the skin model part.
     */
    template <bool TIsInnerLoop>
    void CreateGapAndSkinQuadraturePoints(
        IntegrationParameters& rIntegrationParameters,
        NurbsSurfaceType::Pointer& pNurbsSurface,
        const Node::Pointer& pSurrogateNode1, 
        const Node::Pointer& pSurrogateNode2, 
        const GeometryType::Pointer& rSurrogateBrepMiddleGeometry,
        ModelPart& rIgaModelPart,
        const ModelPart& rSkinSubModelPart);

    /**
     * @brief Populates the surrogate-to-skin projection mapping.
     * @tparam TIsInnerLoop True when working with the inner loop orientation.
     * @param rSurrogateSubModelPart Reference to the surrogate sub model part.
     * @param rSkinSubModelPart Reference to the skin sub model part.
     * @param rSkinNodesPerSpan Reference to the skin node CSR structure.
     */
    template <bool TIsInnerLoop>
    void SetSurrogateToSkinProjections(
        const ModelPart& rSurrogateSubModelPart,
        const ModelPart& rSkinSubModelPart,
        const KnotSpanIdsCSR& rSkinNodesPerSpan);

    /**
     * @brief Checks that surrogate boundary vertices lie on the same skin layer.
     * @param rSkinSubModelPart Reference to the skin sub model part containing the layer information.
     * @param pSurrogateNode1 Pointer to the first surrogate vertex to be validated.
     * @param pSurrogateNode2 Pointer to the second surrogate vertex to be validated.
     */
    void AssestProjectionsFeasibility(
        const ModelPart& rSkinSubModelPart,
        Node::Pointer pSurrogateNode1, 
        Node::Pointer pSurrogateNode2);

    void FindClosestTruePointToSurrogateVertexByNurbs(); // FIXME: unused; remove or implement

    bool ProjectToSkinBoundary( // FIXME: unused; remove or implement
        const ModelPart* pSkinModelPart,
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rProjectedPoint,
        CoordinatesArrayType& rProjectedPointLocal,
        std::vector<array_1d<double, 3>>& rCurveDerivatives,
        int nInitialGuesses);
    
    /**
     * @brief Creates the gap conditions associated with the provided geometries.
     * @param rGeometriesBegin Iterator to the first geometry to be converted.
     * @param rGeometriesEnd Iterator past the last geometry to be converted.
     * @param rModelPart Reference to the model part where the conditions are created.
     * @param rConditionName Reference to the registered condition name.
     * @param rIdCounter Reference to the identifier counter used for new entities.
     * @param pProperties Pointer to the properties assigned to the created conditions.
     * @param KnotSpanSizes Knot span sizes driving the quadrature creation.
     * @param pSurrogateReferenceGeometries Reference geometries used for the surrogate association.
     */
    void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        const std::string& rConditionName,
        std::size_t& rIdCounter,
        PropertiesPointerType pProperties,
        const Vector KnotSpanSizes,
        const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const;

    /// Creates elements from geometries
    /**
     * @brief Creates gap elements from the generated surrogate geometries.
     * @param rGeometriesBegin Iterator to the first source geometry.
     * @param rGeometriesEnd Iterator past the last source geometry.
     * @param rDestinationModelPart Reference to the model part receiving the elements.
     * @param rElementName Reference to the registered element name.
     * @param rIdCounter Reference to the identifier counter used for the elements.
     * @param pProperties Pointer to the material properties applied to the new elements.
     * @param pSurrogateReferenceGeometries Surrogate reference geometries paired with each element.
     */
    void CreateElements(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        const std::string& rElementName,
        std::size_t& rIdCounter,
        PropertiesPointerType pProperties,
        const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const;
    
    /**
     * @brief Creates a Brep curve for the active knot range.
     *
     * @param pFirstBrepPoint Pointer to the first Brep point.
     * @param pSecondBrepPoint Pointer to the second Brep point.
     * @param rActiveRangeKnotVector Reference to the active range knot vector.
     * @return NurbsCurveGeometry<3, PointerVector<Node>>::Pointer Pointer to the generated curve geometry.
     */
    static typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer CreateBrepCurve(
        const Node::Pointer pFirstBrepPoint,
        const Node::Pointer pSecondBrepPoint, 
        const Vector& rActiveRangeKnotVector)
        {
            // Create the data for the trimming curves
            PointerVector<Node> control_points;
            control_points.push_back(pFirstBrepPoint);
            control_points.push_back(pSecondBrepPoint);
            const int polynomial_degree = 1;
            Vector knot_vector = ZeroVector(4) ;
            knot_vector[0] = rActiveRangeKnotVector[0] ;
            knot_vector[1] = rActiveRangeKnotVector[0] ;
            knot_vector[2] = rActiveRangeKnotVector[1] ;
            knot_vector[3] = rActiveRangeKnotVector[1] ;
            // Create the trimming curves
            typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer p_trimming_curve(
                new NurbsCurveGeometry<3, PointerVector<Node>>(
                    control_points,
                    polynomial_degree,
                    knot_vector));   
            return p_trimming_curve;
        }
    
    /// 
    /**  
     * @brief Return n×n Gauss points mapped on the curved quadrilateral
     * 
     * @param Order
     * @param rB0
     * @param rL0
     * @param rL1
     * @param rB1
     * @param rP00
     * @param rP01
     * @param rP10
     * @param rP11
     * @return IntegrationPointsArrayType
     */
    IntegrationPointsArrayType CreateCoonsPatchGaussPoints(
        std::size_t                  Order,
        const BrepCurveType&         rB0,
        const BrepCurveType&         rL0,
        const BrepCurveType&         rL1,
        const BrepCurveType&         rB1,
        const array_1d<double,3>&    rP00,
        const array_1d<double,3>&    rP01,
        const array_1d<double,3>&    rP10,
        const array_1d<double,3>&    rP11) const;

    // --- static utility helpers ----------------------------------

    /// Global coordinates of a point on a Brep-curve at parameter t
    static array_1d<double,3> GlobalPoint(
        const GeometryType& rCurve,
        double               T);

    /// Coons patch mapping X(ξ,η)
    static array_1d<double,3> CoonsPoint(
        double                     Xi,
        double                     Eta,
        const BrepCurveType&       rB0,
        const BrepCurveType&       rL0,
        const BrepCurveType&       rL1,
        const BrepCurveType&       rB1,
        const array_1d<double,3>&  rP00,
        const array_1d<double,3>&  rP01,
        const array_1d<double,3>&  rP10,
        const array_1d<double,3>&  rP11);

    static array_1d<double,3> CoonsDerivative(
        const double                 Xi,
        const double                 Eta,
        const bool                   WithRespectToXi,
        const BrepCurveType&          rB0,   // bordo in eta=0   -> B0(ξ)
        const BrepCurveType&          rL0,   // bordo in xi=0    -> L0(η)
        const BrepCurveType&          rL1,   // bordo in xi=1    -> L1(η)
        const BrepCurveType&          rB1,   // bordo in eta=1   -> B1(ξ)
        const array_1d<double,3>&    rP00,
        const array_1d<double,3>&    rP01,
        const array_1d<double,3>&    rP10,
        const array_1d<double,3>&    rP11);

    /// Finite-difference derivative of CoonsPoint
    static array_1d<double,3> CoonsDerivativeFD(
        const double                     Xi,
        const double                     Eta,
        const bool                       WithRespectToXi,
        const BrepCurveType&       rB0,
        const BrepCurveType&       rL0,
        const BrepCurveType&       rL1,
        const BrepCurveType&       rB1,
        const array_1d<double,3>&  rP00,
        const array_1d<double,3>&  rP01,
        const array_1d<double,3>&  rP10,
        const array_1d<double,3>&  rP11,
        double                     Step = 1.0e-8);


/** Collect skin points between two condition IDs, project them to UV (inclusive). */
template <bool TIsInnerLoop>
std::vector<array_1d<double,3>> CollectSkinUVBetween(
    const ModelPart& rSkinSubModelPart,
    const NurbsSurfaceType& rSurface,
    IndexType id_node_1,
    IndexType id_node_2) const
{
    std::vector<array_1d<double,3>> uv_pts;
    uv_pts.reserve(rSkinSubModelPart.NumberOfConditions());

    const IndexType first_id = rSkinSubModelPart.NodesBegin()->Id();
    const IndexType last_id  = first_id + rSkinSubModelPart.NumberOfNodes() - 1;
    auto next_id = [&](IndexType id){ return (id < last_id) ? (id + 1) : first_id; };
    auto previous_id = [&](IndexType id){ return (id > first_id) ? (id - 1) : last_id; };

    IndexType id = id_node_1;

    int iter = 0;
    while (true) {
        const auto& p_node = rSkinSubModelPart.pGetNode(id);
        array_1d<double,3> uv = p_node->Coordinates();
        uv_pts.push_back(uv);
        if (id == id_node_2) break;

        if constexpr (TIsInnerLoop)  
        { // inner loop: go against the condition orientation
            id = previous_id(id);
        } else { // outer loop: follow the condition orientation
            id = next_id(id);
        }
        iter++;
        KRATOS_ERROR_IF(iter > rSkinSubModelPart.NumberOfConditions())
            << "CollectSkinUVBetween: infinite loop detected between IDs "
            << id_node_1 << " and " << id_node_2 << ".\n";
    }
    return uv_pts;
}



/** Chord-length parametrization on [0,1] for UV points. */
std::vector<double> ChordLengthParams01_UV(const std::vector<array_1d<double,3>>& UV) const
{
    const std::size_t n = UV.size();
    std::vector<double> t(n,0.0);
    if (n <= 1) return t;

    double L = 0.0;
    std::vector<double> d(n,0.0);
    for (std::size_t i=1; i<n; ++i) {
        const double du = UV[i][0]-UV[i-1][0];
        const double dv = UV[i][1]-UV[i-1][1];
        d[i] = std::sqrt(du*du+dv*dv);
        L += d[i];
    }
    if (L <= 1e-16) return t;

    double acc = 0.0;
    for (std::size_t i=1; i<n; ++i) {
        acc += d[i];
        t[i] = acc / L;
    }
    t.front() = 0.0;
    t.back()  = 1.0;
    return t;
}



// ================================================================
// Generic p-degree Bézier/NURBS fitting in UV with fixed endpoints
// ================================================================
//
// Given ordered UV samples {Q_i} between two skin conditions (inclusive),
// fit a p-degree Bézier curve with endpoints fixed to Q_0 and Q_N:
//     C(t) = sum_{k=0}^p B_k^p(t) P_k,  t in [0,1],  weights = 1
// unknowns are the (p-1) interior control points P_1..P_{p-1}.
//
// We solve LS: minimize sum_i ||Q_i - B_0 P_0 - B_p P_p - sum_{k=1}^{p-1} B_k P_k||^2
// Assemble normal equations for U and V independently:
//     (A^T A) X = A^T R,  where A_{i,k} = B_k^p(t_i),  k=1..p-1
//     R_i = Q_i - B_0 P_0 - B_p P_p
// Add a tiny ridge on the diagonal for robustness.
//
// Notes:
// - Parametrization: chord-length in UV (centripetal is also possible).
// - All weights = 1 (polynomial). To switch to rational, extend with per-node weights.
// - Works for p >= 2 (p=1 is trivial polyline).
//
// ================================================================

/** Evaluate Bernstein basis of degree p at t in [0,1]. B has size p+1. */
void BernsteinBasis(const int p, const double t, std::vector<double>& B) const
{
    const double s = std::min(1.0, std::max(0.0, t));
    B.assign(p+1, 0.0);

    // De Casteljau-like stable evaluation for Bernstein basis
    // Start with B[0]=1 and build iteratively.
    B[0] = 1.0;
    const double om = 1.0 - s;
    for (int j = 1; j <= p; ++j) {
        double saved = 0.0;
        for (int k = 0; k < j; ++k) {
            const double tmp = B[k];
            B[k] = saved + om * tmp;
            saved = s * tmp;
        }
        B[j] = saved;
    }
}

/** Build Bézier knot vector [0...0, 1...1] with multiplicity p+1 each side. */
void BuildBezierKnots(const int p, Vector& rKnots) const
{
    const int nkn = 2*(p+1);
    if ((int)rKnots.size() != nkn) rKnots.resize(nkn);
    for (int i=0; i<=p;   ++i) rKnots[i]      = 0.0;
    for (int i=p+1;i<nkn;++i) rKnots[i]      = 1.0;
}

/** Tiny dense solver for symmetric positive definite (normal eq.) with ridge. */
bool SolveSPD(std::vector<std::vector<double>>& A, std::vector<double>& b) const
{
    // Cholesky (naive) with simple pivot guard
    const int n = (int)A.size();
    for (int i=0; i<n; ++i) {
        // Diagonal
        double sum = A[i][i];
        for (int k=0; k<i; ++k) sum -= A[i][k]*A[i][k];
        if (sum <= 1e-20) return false;
        const double Lii = std::sqrt(sum);
        A[i][i] = Lii;

        // Off-diagonals
        for (int j=i+1; j<n; ++j) {
            double s = A[j][i];
            for (int k=0; k<i; ++k) s -= A[j][k]*A[i][k];
            A[j][i] = s / Lii;
        }
    }

    // Solve L y = b
    for (int i=0; i<n; ++i) {
        double s = b[i];
        for (int k=0; k<i; ++k) s -= A[i][k]*b[k];
        b[i] = s / A[i][i];
    }
    // Solve L^T x = y
    for (int i=n-1; i>=0; --i) {
        double s = b[i];
        for (int k=i+1; k<n; ++k) s -= A[k][i]*b[k];
        b[i] = s / A[i][i];
    }
    return true;
}

/** Build a UV Bézier/NURBS of degree p from control points (weights=1). */
typename NurbsCurveGeometryType::Pointer MakeBezierUV_FromControls(
    const std::vector<array_1d<double,3>>& CtrlsUV, const int p) const
{
    const int ncp = p + 1;
    KRATOS_ERROR_IF((int)CtrlsUV.size() != ncp)
        << "MakeBezierUV_FromControls: expected " << ncp << " control points.\n";

    PointerVector<Node> ctrl;
    ctrl.reserve(ncp);
    for (int i=0; i<ncp; ++i) {
        array_1d<double,3> P = CtrlsUV[i];
        P[2] = 0.0; // UV is 2D
        // Virtual node (not stored in a ModelPart)
        ctrl.push_back( Node::Pointer(new Node(1, P)) );
    }

    Vector knots; BuildBezierKnots(p, knots);
    Vector w(ncp); for (int i=0;i<ncp;++i) w[i]=1.0;

    return typename NurbsCurveGeometryType::Pointer(
        new NurbsCurveGeometryType(ctrl, p, knots, w));
}

/** Generic LS fit of degree p for UV samples with fixed endpoints. */
typename NurbsCurveGeometryType::Pointer FitBezierUV_LS_Generic(
    const std::vector<array_1d<double,3>>& UVsamples, // ordered, inclusive
    const int p,
    const double ridge = 1e-12) const
{
    KRATOS_ERROR_IF(p < 2) << "FitBezierUV_LS_Generic: degree p must be >= 2.\n";
    const std::size_t N = UVsamples.size();
    KRATOS_ERROR_IF(N < (std::size_t)(p+1))
        << "FitBezierUV_LS_Generic: need at least p+1 samples.\n";

    // Endpoints fixed
    const array_1d<double,3> P0 = UVsamples.front();
    const array_1d<double,3> Pp = UVsamples.back();

    // Parameterization
    auto t = ChordLengthParams01_UV(UVsamples);

    // Assemble normal equations for interior controls P1..P_{p-1}
    const int m  = p - 1;              // number of unknown control points
    std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.0)); // A^T A
    std::vector<double> bu(m, 0.0), bv(m, 0.0);                         // A^T R (for U and V)

    std::vector<double> B; B.reserve(p+1);

    for (std::size_t i=0; i<N; ++i) {
        BernsteinBasis(p, t[i], B);

        // Residual after removing endpoint contributions
        const double Ru = UVsamples[i][0] - (B[0]*P0[0] + B[p]*Pp[0]);
        const double Rv = UVsamples[i][1] - (B[0]*P0[1] + B[p]*Pp[1]);

        // Fill (A^T A) and (A^T R)
        for (int a=1; a<=p-1; ++a) {
            const double Ba = B[a];
            bu[a-1] += Ba * Ru;
            bv[a-1] += Ba * Rv;
            for (int b=1; b<=p-1; ++b) {
                M[a-1][b-1] += Ba * B[b];
            }
        }
    }

    // Ridge regularization
    for (int d=0; d<m; ++d) M[d][d] += ridge;

    // Solve for U and V components
    auto Mu = M; auto Mv = M; // same matrix for both
    bool ok_u = SolveSPD(Mu, bu);
    bool ok_v = SolveSPD(Mv, bv);
    KRATOS_ERROR_IF(!(ok_u && ok_v)) << "FitBezierUV_LS_Generic: SPD solve failed.\n";

    // Collect control points: P0, P1..P_{p-1}, Pp
    std::vector<array_1d<double,3>> Ctrls(p+1, ZeroVector(3));
    Ctrls[0]     = P0;
    Ctrls[p]     = Pp;
    for (int k=1; k<=p-1; ++k) {
        Ctrls[k][0] = bu[k-1];
        Ctrls[k][1] = bv[k-1];
        Ctrls[k][2] = 0.0;
    }

    return MakeBezierUV_FromControls(Ctrls, p);
}

/** Reverse orientation of a UV Bézier/NURBS (any degree): returns the same geometry, reversed. */
typename NurbsCurveGeometryType::Pointer ReverseBezierUV_Generic(
    const typename NurbsCurveGeometryType::Pointer& p_forward) const
{
    KRATOS_ERROR_IF(!p_forward) << "ReverseBezierUV_Generic: null input.\n";

    const int p = (int)p_forward->PolynomialDegree(0);
    const int n = (int)p_forward->size();
    KRATOS_ERROR_IF(n != p+1) << "ReverseBezierUV_Generic: unexpected control count.\n";

    PointerVector<Node> ctrl_rev; ctrl_rev.reserve(n);
    for (int i=n-1; i>=0; --i) ctrl_rev.push_back(p_forward->pGetPoint(i));

    Vector knots; BuildBezierKnots(p, knots);
    Vector w(n); for (int i=0;i<n;++i) w[i]=1.0;

    return typename NurbsCurveGeometryType::Pointer(
        new NurbsCurveGeometryType(ctrl_rev, p, knots, w));
}

/** Fit between two skin conditions (inclusive) with general degree p in UV. */
template <bool TIsInnerLoop>
typename NurbsCurveGeometryType::Pointer FitUV_BetweenSkinNodes_Generic(
    const ModelPart& rSkinSubModelPart,
    const NurbsSurfaceType& rSurface,
    IndexType id_node_1,
    IndexType id_node_2,
    const int p,
    const double ridge = 1e-12) const
{
    // Reuse your robust UV collector (ordered with wrap-around).
    auto UV = CollectSkinUVBetween<TIsInnerLoop>(rSkinSubModelPart, rSurface, id_node_1, id_node_2);
    KRATOS_ERROR_IF(UV.size() < (std::size_t)(p+1))
        << "FitUV_BetweenSkinConditions_Generic: not enough samples for degree p=" << p << ".\n";

    // Ensure endpoints match exactly the first and last
    UV.front() = UV.front();
    UV.back()  = UV.back();

    return FitBezierUV_LS_Generic(UV, p, ridge);
}

double Orientation(
    const Node::Pointer& p, const Node::Pointer& q, const Node::Pointer& r)
{
    return (q->X() - p->X()) * (r->Y() - p->Y()) -
           (q->Y() - p->Y()) * (r->X() - p->X());
}

bool OnSegment(
    const Node::Pointer& p, const Node::Pointer& q, const Node::Pointer& r)
{
    return (std::min(p->X(), r->X()) <= q->X() && q->X() <= std::max(p->X(), r->X()) &&
            std::min(p->Y(), r->Y()) <= q->Y() && q->Y() <= std::max(p->Y(), r->Y()));
}

bool SegmentsIntersect(
    const Node::Pointer& A, const Node::Pointer& B, const Node::Pointer& C, const Node::Pointer& D)
{
    double o1 = Orientation(A, B, C);
    double o2 = Orientation(A, B, D);
    double o3 = Orientation(C, D, A);
    double o4 = Orientation(C, D, B);

    // Caso generale
    if (o1 * o2 < 0.0 && o3 * o4 < 0.0)
        return true;

    // Collinearità (con bounding box)
    if (std::abs(o1) < 1e-14 && OnSegment(A, C, B)) return true;
    if (std::abs(o2) < 1e-14 && OnSegment(A, D, B)) return true;
    if (std::abs(o3) < 1e-14 && OnSegment(C, A, D)) return true;
    if (std::abs(o4) < 1e-14 && OnSegment(C, B, D)) return true;

    return false;
}

// IndexType FindClosestNodeInLayer( //FIXME: remove
//     const DynamicBinsPointerType& rStartPoint,
//     BinSearchParameters& rSearchParameters,
//     const std::string& rLayer,
//     const ModelPart& rSkinSubModelPart);

template <bool TIsInnerLoop>
IndexType FindClosestNodeInLayerWithDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanIdsCSR& rSkinConditionsPerSpan,
    const Vector& rDirection);

        
}; // Class SnakeGapSbmProcess

}  // namespace Kratos.
