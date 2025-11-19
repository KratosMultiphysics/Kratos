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
#include <algorithm>

// Project includes
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"
#include "snake_sbm_process.h"
#include "custom_utilities/create_breps_sbm_utilities.h"
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

    /**
     * @brief Compact CSR container used to associate knot spans with node or condition identifiers.
     */
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

    /**
     * @brief Lightweight view exposing a contiguous block of identifiers.
     */
    struct IdsView {
        const IndexType* data = nullptr;
        std::size_t size = 0;
    };

    /**
     * @brief Returns the CSR non-zero index for a given matrix entry.
     * @param A Constant reference to the sparse matrix inspected.
     * @param i Row index requested for the lookup.
     * @param j Column index requested for the lookup.
     * @return The flattened CSR index when the entry exists, or static_cast<std::size_t>(-1) otherwise.
     */
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


    /**
     * @brief Packs per-cell node identifiers into the CSR payload.
     * @param tmp Temporary storage containing unsorted node identifiers per cell.
     * @param out Reference to the CSR container receiving the compacted data.
     * @details The entries are sorted and deduplicated prior to being appended to the pool.
     */
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

    /**
     * @brief Provides the node identifiers associated with a CSR entry index.
     * @param S Reference to the CSR structure storing the identifiers.
     * @param k Flat CSR index representing the queried cell.
     * @return Lightweight view on the node identifiers stored in the entry.
     */
    static IdsView CellIdsByK(const KnotSpanIdsCSR& S, const std::size_t k)
    {
        if (k >= S.nnz_off.size() || k >= S.nnz_len.size()) {
            KRATOS_WARNING_ONCE("SnakeGapSbmProcess")
                << "CellIdsByK received index " << k << " but nnz_off has size "
                << S.nnz_off.size() << " and nnz_len has size " << S.nnz_len.size() << ". "
                << "Returning empty view to avoid accessing invalid data." << std::endl;
            return {};
        }

        const std::size_t offset = S.nnz_off[k];
        const std::size_t length = S.nnz_len[k];
        const std::size_t pool_size = S.pool.size();

        if (offset >= pool_size) {
            KRATOS_WARNING_ONCE("SnakeGapSbmProcess")
                << "CellIdsByK computed offset " << offset << " beyond pool size "
                << pool_size << ". Returning empty view to avoid invalid access." << std::endl;
            return {};
        }

        const std::size_t max_length = pool_size - offset;
        const std::size_t safe_length = std::min(length, max_length);
        if (safe_length != length) {
            KRATOS_WARNING_ONCE("SnakeGapSbmProcess")
                << "CellIdsByK clipped entry length from " << length << " to "
                << safe_length << " to stay within pool bounds (" << pool_size << ")." << std::endl;
        }

        return IdsView{ S.pool.data() + offset, safe_length };
    }

    /**
     * @brief Retrieves the node identifiers stored for a specific cell location.
     * @param S Reference to the CSR structure storing the identifiers.
     * @param i Row index of the inspected cell.
     * @param j Column index of the inspected cell.
     * @return Lightweight view on the node identifiers, or an empty view when the cell is empty.
     */
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

    /**
     * @brief Constructs the SnakeGapSbmProcess with the provided model and parameters.
     * @param rModel Reference to the parent model instance.
     * @param ThisParameters Parameter container describing the SBM configuration.
     */
    SnakeGapSbmProcess(
        Model& rModel,
        Parameters ThisParameters);

    /**
     * @brief Defaulted destructor.
     */
    ~SnakeGapSbmProcess() = default;

    /**
     * @brief Retrieves the default parameter set for the process.
     * @return Parameter container holding the default values.
     */
    const Parameters GetDefaultParameters() const override;

    /**
     * @brief Retrieves the validated parameter set for the process.
     * @return Parameter container with validated entries.
     */
    const Parameters GetValidParameters() const;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Executes the process and builds the extended gap SBM geometries.
     */
    void Execute() override
    {   
        CreateSbmExtendedGeometries();
    };
    
    /**
     * @brief Prepares the process before the simulation starts.
     */
    void ExecuteInitialize() override
    {
        SnakeSbmProcess::CreateTheSnakeCoordinates();
    };

    /**
     * @brief Initializes data required at the beginning of the solution step.
     */
    void ExecuteInitializeSolutionStep() override
    {};

    /**
     * @brief Executes tasks required before entering the solution loop.
     */
    void ExecuteBeforeSolutionLoop() override
    {};

    /**
     * @brief Builds a CSR matrix linking skin nodes to knot span bins.
     * @param rSkinSubModelPart Reference to the skin model part.
     * @param rSurrogateSubModelPart Reference to the surrogate model part used for binning.
     * @return CSR structure collecting the node identifiers per knot span.
     */
    KnotSpanIdsCSR CreateSkinNodesPerKnotSpanMatrix(
        const ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart) const;

    /**
     * @brief Builds a CSR matrix linking skin conditions to knot span bins.
     * @param rSkinSubModelPart Reference to the skin model part containing the conditions.
     * @param rReferenceMatrix Reference to the node-based CSR used as topology template.
     * @return CSR structure collecting the condition identifiers per knot span.
     */
    KnotSpanIdsCSR CreateSkinConditionsPerKnotSpanMatrix(
        const ModelPart& rSkinSubModelPart,
        const SnakeGapSbmProcess::KnotSpanIdsCSR& rReferenceMatrix) const;


/**
 * @brief Aggregates integration settings shared across gap SBM computations.
 */
struct IntegrationParameters
{
    std::size_t  NumberOfShapeFunctionsDerivatives;
    IntegrationInfo CurveIntegrationInfo;
    const Vector&           rKnotSpanSizes;
    const KnotSpanIdsCSR* pSkinNodesPerSpan = nullptr;
    const KnotSpanIdsCSR* pSkinConditionsPerSpan = nullptr;
    /**
     * @brief Initializes the integration parameters container.
     * @param NumberOfDerivatives Number of shape function derivatives required.
     * @param rInfo Reference to the IntegrationInfo used to build quadrature.
     * @param rKnotSpans Reference to the vector storing knot span lengths.
     */
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
    std::size_t mGapInterpolationOrder;
    std::string mGapSbmType;
    bool mUseForMultipatch = false;

    /**
     * @brief Attaches the surrogate middle geometry to skin nodes.
     * @param rSkinSubModelPart Skin model part containing the nodes.
     * @param rSegmentProjectionIds Projection node ids for each surrogate segment node.
     * @param pSurrogateMiddleGeometry Surrogate middle brep geometry.
     * @param pProjectionNode1 First projection node on the skin.
     * @param pProjectionNode2 Second projection node on the skin.
     */
    void AttachSurrogateMiddleGeometryToSkinNodes(
        const ModelPart& rSkinSubModelPart,
        const std::vector<IndexType>& rSegmentProjectionIds,
        GeometryType::Pointer pSurrogateMiddleGeometry,
        const NodeType::Pointer& pProjectionNode1,
        const NodeType::Pointer& pProjectionNode2);

    /**
     * @brief Synchronizes neighbour geometries of the first and last skin nodes.
     * @param rSkinSubModelPart Skin model part whose end nodes will be updated.
     */
    void SynchronizeEndSkinNodeNeighbourGeometries(const ModelPart& rSkinSubModelPart);

    /**
     * @brief Creates multipatch coupling conditions on the inner skin.
     * @param rSkinSubModelPart Skin model part containing the inner-loop conditions.
     * @param rKnotSpanSizes Knot-span sizes used to scale integration data.
     * @param pNurbsSurface Pointer to the NURBS surface hosting the SBM geometry.
     */
    void CreateInnerSkinMultipatchCouplingConditions(
        const ModelPart& rSkinSubModelPart,
        const Vector& rKnotSpanSizes,
        NurbsSurfaceType::Pointer& pNurbsSurface);

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
    
    /**
     * @brief Creates the gap conditions associated with the provided geometries.
     * @param rGeometriesBegin Iterator to the first geometry to be converted.
     * @param rGeometriesEnd Iterator past the last geometry to be converted.
     * @param rModelPart Reference to the model part where the conditions are created.
     * @param rConditionName Reference to the registered condition name.
     * @param rIdCounter Reference to the identifier counter used for new entities.
     * @param pProperties Pointer to the properties assigned to the created conditions.
     * @param KnotSpanSizes Vector storing the knot span sizes guiding the quadrature creation.
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
    
    /**
     * @brief Returns Gauss points mapped on the curved quadrilateral defined by the Coons patch.
     * @param Order Quadrature order used in each parametric direction.
     * @param rB0 Brep curve describing the edge at eta = 0.
     * @param rL0 Brep curve describing the edge at xi = 0.
     * @param rL1 Brep curve describing the edge at xi = 1.
     * @param rB1 Brep curve describing the edge at eta = 1.
     * @param rP00 Corner point at (xi, eta) = (0,0).
     * @param rP01 Corner point at (xi, eta) = (0,1).
     * @param rP10 Corner point at (xi, eta) = (1,0).
     * @param rP11 Corner point at (xi, eta) = (1,1).
     * @return Array of integration points mapped onto the Coons surface.
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

    /**
     * @brief Computes the global coordinates of a point on a Brep curve at the requested parameter.
     * @param rCurve Reference to the Brep curve geometry.
     * @param T Parametric coordinate where the curve is evaluated.
     * @return Global coordinates corresponding to the evaluated parameter.
     */
    static array_1d<double,3> GlobalPoint(
        const GeometryType& rCurve,
        double               T);

    /**
     * @brief Evaluates the Coons patch mapping at the given parametric location.
     * @param Xi Parametric coordinate along the first direction.
     * @param Eta Parametric coordinate along the second direction.
     * @param rB0 Brep curve representing the edge at eta = 0.
     * @param rL0 Brep curve representing the edge at xi = 0.
     * @param rL1 Brep curve representing the edge at xi = 1.
     * @param rB1 Brep curve representing the edge at eta = 1.
     * @param rP00 Corner point at (xi, eta) = (0,0).
     * @param rP01 Corner point at (xi, eta) = (0,1).
     * @param rP10 Corner point at (xi, eta) = (1,0).
     * @param rP11 Corner point at (xi, eta) = (1,1).
     * @return Global coordinates of the Coons patch at (Xi, Eta).
     */
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

    /**
     * @brief Computes the derivative of the Coons patch with respect to xi or eta.
     * @param Xi Parametric coordinate along the first direction.
     * @param Eta Parametric coordinate along the second direction.
     * @param WithRespectToXi True to differentiate with respect to xi, false for eta.
     * @param rB0 Brep curve representing the edge at eta = 0.
     * @param rL0 Brep curve representing the edge at xi = 0.
     * @param rL1 Brep curve representing the edge at xi = 1.
     * @param rB1 Brep curve representing the edge at eta = 1.
     * @param rP00 Corner point at (xi, eta) = (0,0).
     * @param rP01 Corner point at (xi, eta) = (0,1).
     * @param rP10 Corner point at (xi, eta) = (1,0).
     * @param rP11 Corner point at (xi, eta) = (1,1).
     * @return Derivative vector evaluated at (Xi, Eta).
     */
    static array_1d<double,3> CoonsDerivative(
        const double                 Xi,
        const double                 Eta,
        const bool                   WithRespectToXi,
        const BrepCurveType&          rB0,
        const BrepCurveType&          rL0,
        const BrepCurveType&          rL1,
        const BrepCurveType&          rB1,
        const array_1d<double,3>&    rP00,
        const array_1d<double,3>&    rP01,
        const array_1d<double,3>&    rP10,
        const array_1d<double,3>&    rP11);

    /**
     * @brief Approximates the Coons patch derivative by finite differences.
     * @param Xi Parametric coordinate along the first direction.
     * @param Eta Parametric coordinate along the second direction.
     * @param WithRespectToXi True to differentiate with respect to xi, false for eta.
     * @param rB0 Brep curve representing the edge at eta = 0.
     * @param rL0 Brep curve representing the edge at xi = 0.
     * @param rL1 Brep curve representing the edge at xi = 1.
     * @param rB1 Brep curve representing the edge at eta = 1.
     * @param rP00 Corner point at (xi, eta) = (0,0).
     * @param rP01 Corner point at (xi, eta) = (0,1).
     * @param rP10 Corner point at (xi, eta) = (1,0).
     * @param rP11 Corner point at (xi, eta) = (1,1).
     * @param Step Increment used for the finite difference approximation.
     * @return Approximate derivative vector at (Xi, Eta).
     */
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


/**
 * @brief Collects skin nodes between two indices and projects them to the UV space.
 * @tparam TIsInnerLoop True when traversing the inner loop orientation.
 * @param rSkinSubModelPart Reference to the skin sub model part considered.
 * @param rSurface Reference to the surface used for the projection.
 * @param id_node_1 Identifier of the starting skin node.
 * @param id_node_2 Identifier of the ending skin node.
 * @return Ordered list of projected UV coordinates.
 */
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

    ModelPart::SizeType iter = 0;
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



/**
 * @brief Computes the chord-length parametrization on [0,1] for UV samples.
 * @param UV Collection of UV points describing the sampled curve.
 * @return Normalized cumulative arc-length parameters.
 */
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

/**
 * @brief Evaluates the Bernstein basis of degree p at a parameter in [0,1].
 * @param p Polynomial degree of the basis.
 * @param t Parametric value where the basis is evaluated.
 * @param B Vector storing the evaluated basis functions.
 */
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

/**
 * @brief Builds a Bézier knot vector with multiplicity p+1 at the boundaries.
 * @param p Polynomial degree of the Bézier curve.
 * @param rKnots Vector receiving the generated knot values.
 */
void BuildBezierKnots(const int p, Vector& rKnots) const
{
    const int nkn = 2*(p+1);
    if ((int)rKnots.size() != nkn) rKnots.resize(nkn);
    for (int i=0; i<=p;   ++i) rKnots[i]      = 0.0;
    for (int i=p+1;i<nkn;++i) rKnots[i]      = 1.0;
}

/**
 * @brief Solves a symmetric positive definite linear system using a simple Cholesky factorization.
 * @param A Reference to the system matrix, overwritten with its factorization.
 * @param b Right-hand side vector, overwritten with the solution.
 * @return True when the system is successfully solved.
 */
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

/**
 * @brief Creates a UV Bézier curve of degree p from a set of control points.
 * @param CtrlsUV Control points expressed in UV coordinates.
 * @param p Polynomial degree of the resulting curve.
 * @return Pointer to the constructed Bézier curve geometry.
 */
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

/**
 * @brief Performs a least-squares Bézier fit of degree p for UV samples with fixed endpoints.
 * @param UVsamples Ordered UV samples, inclusive of both endpoints.
 * @param p Polynomial degree of the fitted curve.
 * @param ridge Ridge regularization magnitude applied to the normal equations.
 * @return Pointer to the fitted Bézier curve geometry.
 */
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

/**
 * @brief Reverses the orientation of a UV Bézier curve.
 * @param p_forward Pointer to the curve whose orientation is reversed.
 * @return Pointer to a curve geometry with inverted control point order.
 */
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

/**
 * @brief Fits a UV Bézier curve between two skin nodes with the requested degree.
 * @tparam TIsInnerLoop True when traversing the inner loop orientation.
 * @param rSkinSubModelPart Reference to the skin sub model part containing the nodes.
 * @param rSurface Reference to the surface used for UV projection.
 * @param id_node_1 Identifier of the first node of the fitting interval.
 * @param id_node_2 Identifier of the last node of the fitting interval.
 * @param p Polynomial degree of the Bézier curve.
 * @param ridge Ridge regularization applied to the fitting system.
 * @return Pointer to the fitted Bézier curve geometry.
 */
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

/**
 * @brief Computes the orientation sign of the triangle (p, q, r).
 * @param p Pointer to the first point.
 * @param q Pointer to the second point.
 * @param r Pointer to the third point.
 * @return Signed area representing the orientation.
 */
double Orientation(
    const Node::Pointer& p, const Node::Pointer& q, const Node::Pointer& r)
{
    return (q->X() - p->X()) * (r->Y() - p->Y()) -
           (q->Y() - p->Y()) * (r->X() - p->X());
}

/**
 * @brief Checks whether a point lies on the segment joining two other points.
 * @param p Pointer to the first endpoint.
 * @param q Pointer to the point to be tested.
 * @param r Pointer to the second endpoint.
 * @return True when q lies within the bounding box of segment pr.
 */
bool OnSegment(
    const Node::Pointer& p, const Node::Pointer& q, const Node::Pointer& r)
{
    return (std::min(p->X(), r->X()) <= q->X() && q->X() <= std::max(p->X(), r->X()) &&
            std::min(p->Y(), r->Y()) <= q->Y() && q->Y() <= std::max(p->Y(), r->Y()));
}

/**
 * @brief Determines if two segments intersect in the plane.
 * @param A Pointer to the first endpoint of the first segment.
 * @param B Pointer to the second endpoint of the first segment.
 * @param C Pointer to the first endpoint of the second segment.
 * @param D Pointer to the second endpoint of the second segment.
 * @return True when the segments intersect.
 */
bool SegmentsIntersect(
    const Node::Pointer& A, const Node::Pointer& B, const Node::Pointer& C, const Node::Pointer& D)
{
    double o1 = Orientation(A, B, C);
    double o2 = Orientation(A, B, D);
    double o3 = Orientation(C, D, A);
    double o4 = Orientation(C, D, B);

    // General case
    if (o1 * o2 < 0.0 && o3 * o4 < 0.0)
        return true;

    // Collinearity (with bounding box check)
    if (std::abs(o1) < 1e-14 && OnSegment(A, C, B)) return true;
    if (std::abs(o2) < 1e-14 && OnSegment(A, D, B)) return true;
    if (std::abs(o3) < 1e-14 && OnSegment(C, A, D)) return true;
    if (std::abs(o4) < 1e-14 && OnSegment(C, B, D)) return true;

    return false;
}

/**
 * @brief Determines if two segments intersect in the plane and returns the intersection point.
 * @param A Pointer to the first endpoint of the first segment.
 * @param B Pointer to the second endpoint of the first segment.
 * @param C Pointer to the first endpoint of the second segment.
 * @param D Pointer to the second endpoint of the second segment.
 * @param rIntersection Output: intersection coordinates (valid only if function returns true).
 * @return True when the segments intersect.
 */
bool SegmentsIntersect(
    const Node::Pointer& A, const Node::Pointer& B,
    const Node::Pointer& C, const Node::Pointer& D,
    array_1d<double,3>& rIntersection)
{
    // Vector representation in 2D
    const double x1 = A->X(), y1 = A->Y();
    const double x2 = B->X(), y2 = B->Y();
    const double x3 = C->X(), y3 = C->Y();
    const double x4 = D->X(), y4 = D->Y();

    const double dx12 = x2 - x1;
    const double dy12 = y2 - y1;
    const double dx34 = x4 - x3;
    const double dy34 = y4 - y3;

    const double denom = dx12 * dy34 - dy12 * dx34;
    const double eps = 1.0e-12;

    // Parallel (including collinear) case
    if (std::abs(denom) < eps) {
        // Collinear check via orientation
        const Node::Pointer pA = A, pB = B, pC = C, pD = D;
        const double o1 = Orientation(pA, pB, pC);
        const double o2 = Orientation(pA, pB, pD);
        if (std::abs(o1) < eps || std::abs(o2) < eps) {
            // Overlapping or touching: return a shared endpoint if any
            if (OnSegment(pA, pC, pB)) { rIntersection[0]=x3; rIntersection[1]=y3; rIntersection[2]=0.0; return true; }
            if (OnSegment(pA, pD, pB)) { rIntersection[0]=x4; rIntersection[1]=y4; rIntersection[2]=0.0; return true; }
            if (OnSegment(pC, pA, pD)) { rIntersection[0]=x1; rIntersection[1]=y1; rIntersection[2]=0.0; return true; }
            if (OnSegment(pC, pB, pD)) { rIntersection[0]=x2; rIntersection[1]=y2; rIntersection[2]=0.0; return true; }
        }
        return false;
    }

    // Proper intersection: solve for t and u in A + t*(B-A) = C + u*(D-C)
    const double t = ((x3 - x1) * dy34 - (y3 - y1) * dx34) / denom;
    const double u = ((x3 - x1) * dy12 - (y3 - y1) * dx12) / denom;

    if (t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps) {
        rIntersection[0] = x1 + t * dx12;
        rIntersection[1] = y1 + t * dy12;
        // Preserve planar assumption; set Z from linear interpolation if available
        const double z1 = A->Z();
        const double z2 = B->Z();
        rIntersection[2] = z1 + t * (z2 - z1);
        return true;
    }
    return false;
}

/**
 * @brief Ray (AB) vs segment (CD) intersection in the plane, returns intersection point.
 *        The first primitive is treated as a half-infinite ray starting at A and going through B.
 */
bool SegmentsIntersectRay(
    const Node::Pointer& A, const Node::Pointer& B,
    const Node::Pointer& C, const Node::Pointer& D,
    array_1d<double,3>& rIntersection)
{
    // Vector representation in 2D
    const double x1 = A->X(), y1 = A->Y();
    const double x2 = B->X(), y2 = B->Y();
    const double x3 = C->X(), y3 = C->Y();
    const double x4 = D->X(), y4 = D->Y();

    const double dx12 = x2 - x1;
    const double dy12 = y2 - y1;
    const double dx34 = x4 - x3;
    const double dy34 = y4 - y3;

    const double denom = dx12 * dy34 - dy12 * dx34;
    const double eps = 1.0e-12;

    // Parallel (including collinear) case
    if (std::abs(denom) < eps) {
        // Collinear check via orientation
        const Node::Pointer pA = A, pB = B, pC = C, pD = D;
        const double o1 = Orientation(pA, pB, pC);
        const double o2 = Orientation(pA, pB, pD);
        if (std::abs(o1) < eps || std::abs(o2) < eps) {
            // Overlapping or touching: prefer the closest endpoint along the ray direction
            // Test C
            if (OnSegment(pA, pC, pB)) { rIntersection[0]=x3; rIntersection[1]=y3; rIntersection[2]=0.0; return true; }
            // Test D
            if (OnSegment(pA, pD, pB)) { rIntersection[0]=x4; rIntersection[1]=y4; rIntersection[2]=0.0; return true; }
        }
        return false;
    }

    // Proper intersection: solve for t and u in A + t*(B-A) = C + u*(D-C)
    const double t = ((x3 - x1) * dy34 - (y3 - y1) * dx34) / denom;
    const double u = ((x3 - x1) * dy12 - (y3 - y1) * dx12) / denom;

    // Finite-segment condition on AB to avoid passing beyond B: t in [0,1] (within tolerance)
    // Segment condition on CD: u in [0,1] within tolerance
    if (t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps) {
        rIntersection[0] = x1 + t * dx12;
        rIntersection[1] = y1 + t * dy12;
        // Preserve planar assumption; set Z from linear interpolation if available
        const double z1 = A->Z();
        const double z2 = B->Z();
        rIntersection[2] = z1 + t * (z2 - z1);
        return true;
    }
    return false;
}

/**
 * @brief Finds the closest node belonging to the target layer by following a prescribed direction.
 * @tparam TIsInnerLoop True when searching along the inner loop orientation.
 * @param rStartPoint Starting point used for the directional search.
 * @param rLayer Name of the skin layer to be matched.
 * @param rSkinSubModelPart Reference to the skin sub model part containing the candidate nodes.
 * @param rKnotSpanSizes Reference to the knot span sizes guiding the search.
 * @param rSkinConditionsPerSpan CSR matrix storing skin conditions per knot span.
 * @param rDirection Direction vector used for ray casting.
 * @return Identifier of the closest node fulfilling the criteria.
 */
template <bool TIsInnerLoop>
IndexType FindClosestNodeInLayerWithDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanIdsCSR& rSkinConditionsPerSpan,
    const Vector& rDirection);

/**
 * @brief Finds the closest node in the target layer along a given direction, skipping
 *        intersections that would form an angle of ~180 degrees with respect to a
 *        reference arm around a given center (used to avoid reflex angles across surrogate).
 * @tparam TIsInnerLoop True when searching along the inner loop orientation.
 * @param rStartPoint Starting point used for the directional search.
 * @param rLayer Name of the skin layer to be matched.
 * @param rSkinSubModelPart Reference to the skin sub model part containing the candidate nodes.
 * @param rKnotSpanSizes Reference to the knot span sizes guiding the search.
 * @param rSkinConditionsPerSpan CSR matrix storing skin conditions per knot span.
 * @param rDirection Direction vector used for ray casting (tangent direction).
 * @param rAngleCenter Center point about which the angle is measured.
 * @param rFixedArmPoint Point defining the fixed arm of the angle (typically p_first).
 * @param angle_eps Numerical tolerance for detecting 180-degree configurations.
 * @return Identifier of the closest node fulfilling the criteria, or max IndexType on failure.
 */
template <bool TIsInnerLoop>
IndexType FindClosestNodeInLayerWithDirectionAngleGuard(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanIdsCSR& rSkinConditionsPerSpan,
    const Vector& rDirection,
    const array_1d<double,3>& rAngleCenter,
    const array_1d<double,3>& rFixedArmPoint,
    const double angle_eps);

// starts from a coordinate and returns both endpoints of the intersected skin condition
template <bool TIsInnerLoop>
std::pair<IndexType, IndexType> FindClosestPairInLayerWithNormalDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanIdsCSR& rSkinConditionsPerSpan,
    const Vector& rDirection);

        
}; // Class SnakeGapSbmProcess

}  // namespace Kratos.
