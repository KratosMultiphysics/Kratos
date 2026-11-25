//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Andrea Gorgi
//                   Nicoló Antonelli
//
// System includes

// Project includes
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "snake_sbm_process.h"
#include "custom_utilities/create_breps_sbm_utilities.h"
#include "spaces/ublas_space.h"
#include "integration/integration_point_utilities.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"

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

    static constexpr double GeometricTolerance = 1.0e-12;

    using NodeType = Node;
    using ContainerNodeType = PointerVector<Node>;
    using ContainerEmbeddedNodeType = PointerVector<Point>;
    using BrepCurveOnSurfaceType = BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType>;
    using BrepSurfaceOnVolumeType = BrepSurfaceOnVolume<ContainerNodeType, true, ContainerNodeType>;
    using ElementsContainerType = ModelPart::ElementsContainerType;

    using GeometryType = Geometry<NodeType>;
    using GeometriesArrayType = GeometryType::GeometriesArrayType;
    using IntegrationPointsArrayType = GeometryType::IntegrationPointsArrayType;
    using CoordinatesArrayType = GeometryType::CoordinatesArrayType;
    using PropertiesPointerType = Properties::Pointer;
    using BrepCurveType = BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType>;
    
    using NurbsCurveGeometryType = NurbsCurveGeometry<3, PointerVector<Node>>;
    using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<NodeType>>;
    using NurbsVolumeType = NurbsVolumeGeometry<PointerVector<NodeType>>;
    using NodePointerVector = GlobalPointersVector<NodeType>;
    using ConditionType = Condition;
    using ConditionPointerType = ConditionType::Pointer;
    using ConditionPointerContainerType = std::vector<ConditionPointerType>;
    using NodePointerContainerType = std::vector<NodeType::Pointer>;
    using NodeBinsType = BinsDynamic<3, NodeType, NodePointerContainerType>;

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
     * @brief Configure for condition-based bins (AABB per condition geometry).
     * @details Bounding boxes, centers and distances use current nodal coordinates.
     * The gap modeler builds and searches these bins during geometry setup, before
     * the structural solve, using the newly created auxiliary skin conditions.
     * The bins must be rebuilt if those nodes move before a subsequent search.
     */
    struct ConditionConfigure
    {
        static constexpr auto Epsilon = std::numeric_limits<double>::epsilon();
        static constexpr auto Dimension = 3;

        using PointType = Point;
        using ObjectType = ConditionType;
        using PointerType = ConditionPointerType;

        using ObjectContainerType = ConditionPointerContainerType;
        using ContainerType = ObjectContainerType;
        using ResultContainerType = ObjectContainerType;
        using IteratorType = ContainerType::iterator;
        using ResultIteratorType = ResultContainerType::iterator;
        using DistanceIteratorType = std::vector<double>::iterator;

        using CoordinateType = double;
        using CoordinateArray = Tvector<CoordinateType, Dimension>;

        static void CalculateBoundingBox(
            const PointerType& rObject,
            PointType& rLowPoint,
            PointType& rHighPoint);

        static void CalculateBoundingBox(
            const PointerType& rObject,
            PointType& rLowPoint,
            PointType& rHighPoint,
            const double Radius);

        static void CalculateCenter(
            const PointerType& rObject,
            PointType& rCentralPoint);

        static bool Intersection(
            const PointerType& rObj_1,
            const PointerType& rObj_2);

        static bool Intersection(
            const PointerType& rObj_1,
            const PointerType& rObj_2,
            const double Radius);

        static bool IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint);

        static bool IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint,
            const double Radius);

        static void Distance(
            const PointerType& rObj_1,
            const PointerType& rObj_2,
            double& distance);

    };

    /**
     * @brief CSR container storing skin node and condition bins per knot span.
     */
    struct KnotSpanSkinBinsCSR
    {
        std::size_t NumberOfSpansX = 0;
        std::size_t NumberOfSpansY = 0;
        std::size_t NumberOfSpansZ = 1;
        double MinU = 0.0, MaxU = 0.0;
        double MinV = 0.0, MaxV = 0.0;
        double MinW = 0.0, MaxW = 0.0;
        double SpanSizeX = 0.0, SpanSizeY = 0.0, SpanSizeZ = 0.0;

        // In 3D the matrix column is the flattened (y,z) span index.
        SparseMatrixType Occupancy;

        struct CellBins
        {
            NodePointerContainerType Nodes;
            ConditionPointerContainerType Conditions;
            std::unique_ptr<NodeBinsType> pNodeBins;
            BinsObjectDynamic<ConditionConfigure> ConditionBins;
            bool HasNodeBins = false;
            bool HasConditionBins = false;
        };

        std::vector<CellBins> CellBinsByNnz;
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
    static std::size_t FindNnzIndex(const SparseMatrixType& A, std::size_t i, std::size_t j);

    static std::size_t FlattenSpanColumnIndex(
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const std::size_t SpanIndexY,
        const std::size_t SpanIndexZ = 0);

    static std::size_t FindSpanNnzIndex(
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const std::size_t SpanIndexX,
        const std::size_t SpanIndexY,
        const std::size_t SpanIndexZ = 0);


    /**
     * @brief Packs per-cell node identifiers into the CSR payload.
     * @param tmp Temporary storage containing unsorted node identifiers per cell.
     * @param out Reference to the CSR container receiving the compacted data.
     * @details The entries are sorted and deduplicated prior to being appended to the pool.
     */
    static void CommitPayload(std::vector<std::vector<IndexType>>& tmp,
                                    KnotSpanIdsCSR& out);

    /**
     * @brief Provides the node identifiers associated with a CSR entry index.
     * @param S Reference to the CSR structure storing the identifiers.
     * @param k Flat CSR index representing the queried cell.
     * @return Lightweight view on the node identifiers stored in the entry.
     */
    static IdsView CellIdsByK(const KnotSpanIdsCSR& S, const std::size_t k);

    /**
     * @brief Retrieves the node identifiers stored for a specific cell location.
     * @param S Reference to the CSR structure storing the identifiers.
     * @param i Row index of the inspected cell.
     * @param j Column index of the inspected cell.
     * @return Lightweight view on the node identifiers, or an empty view when the cell is empty.
     */
    static IdsView CellIds(const KnotSpanIdsCSR& S, const std::size_t i, const std::size_t j);

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
//-------------------------------
    void Execute() override
    {   
        CreateSbmExtendedGeometries();
    };
    
    /**
     * @brief Prepares the process before the simulation starts.
     */
//-------------------------------
    void ExecuteInitialize() override
    {
        SnakeSbmProcess::CreateTheSnakeCoordinates(true);
    };

    /**
     * @brief Initializes data required at the beginning of the solution step.
     */
//-------------------------------
    void ExecuteInitializeSolutionStep() override
    {};

    /**
     * @brief Executes tasks required before entering the solution loop.
     */
//-------------------------------
    void ExecuteBeforeSolutionLoop() override
    {};

    /**
     * @brief Builds a 2D CSR matrix linking skin nodes and conditions to knot span bins.
     * @param rSkinSubModelPart Reference to the skin model part.
     * @param rSurrogateSubModelPart Reference to the surrogate model part used for binning.
     * @return CSR structure collecting node and condition bins per knot span.
     */
    KnotSpanSkinBinsCSR CreateSkinBinsPerKnotSpanMatrix2D(
        const ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart) const;



/**
 * @brief Aggregates integration settings shared across gap SBM computations.
 */
struct IntegrationParameters
{
    std::size_t  NumberOfShapeFunctionsDerivatives;
    IntegrationInfo CurveIntegrationInfo;
    const Vector&           rKnotSpanSizes;
    const KnotSpanSkinBinsCSR* pSkinBinsPerSpan = nullptr;
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

/**
 * @brief Useful struct to store the result of projections
 * 
 */
struct ProjectionResult
{
    IndexType ProjectionId = std::numeric_limits<IndexType>::max();
    CoordinatesArrayType ProjectionPoint = ZeroVector(3);
    bool SuccessfulProjection = false;
};
    
private:

    ModelPart* mpGapElementsSubModelPart = nullptr;
    ModelPart* mpGapConditionsSubModelPart = nullptr;
    ModelPart* mpGapInterfaceSubModelPart = nullptr;
    std::string mGapElementName;
    std::string mGapConditionName;
    std::string mGapInterfaceConditionName;
    std::string mBackgroundSurfaceGeometryName = "IgaModelPart_snake_background_surface";
    std::size_t mInternalDivisions;
    double mGapRelativeToleranceForSubdivisions = 0.1;
    double mMinSpanSize = 0.0;
    std::size_t mNumberOfInterpolationLevels = 3; // number of additional levels of interpolation for the gap elements (in addition to the corner nodes)
    std::size_t mGapApproximationOrder = 1; // Order of quadrature of the Gap elemennts and conditions
    std::size_t mGapInterpolationOrder; // Order of the Gap elements
    std::string mGapSbmType; // "default" (uses ordered points of the curve) or
                             // "interpolation" (interpolates the points on the curve -- no assumption on their order)
                             // "sbm" (uses SBM on top of the "interpolation" approach to further optimize the gap geometry)

    /**
     * @brief Creates the gap SBM geometries for the configured model part.
     */
    void CreateSbmExtendedGeometries();

    /**
     * @brief Creates the gap SBM geometries for the configured model part. 2D Case
     */
    void CreateSbmExtendedGeometries2D();

    /**
     * @brief Creates the gap SBM geometries for the configured model part. 3D Case
     */
    void CreateSbmExtendedGeometries3D();

    /**
     * @brief Builds the gap SBM geometries for a given skin/surrogate pair.
     * @tparam TIsInnerLoop True when processing the inner loop orientation.
     * @param rSkinSubModelPart Reference to the target skin sub model part.
     * @param rSurrogateSubModelPart Reference to the surrogate sub model part.
     */
    template <bool TIsInnerLoop>
    void CreateSbmExtendedGeometries2D(
        ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart);

    /// Create the BREP midpoint geometries and gap segments for one surrogate loop.
    template <bool TIsInnerLoop>
    void CreateGapAlongSurrogateBreps2D(
        ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart,
        IndexType StartingBrepId,
        std::size_t SizeSurrogateLoop,
        std::size_t NumberOfShapeFunctionsDerivatives,
        const Vector& rKnotSpanSizes,
        const KnotSpanSkinBinsCSR& rSkinBinsPerKnotSpan,
        NurbsSurfaceType::Pointer pNurbsSurface);

    /// Build a projection chain and its gap quadrature for one surrogate BREP.
    template <bool TIsInnerLoop>
    void CreateGapForSurrogateSegment2D(
        const Node::Pointer& rpSurrogateNode1,
        const Node::Pointer& rpSurrogateNode2,
        ModelPart& rSkinSubModelPart,
        const array_1d<double, 3>& rSurrogateTangent,
        const array_1d<double, 3>& rSurrogateNormal,
        IntegrationParameters& rIntegrationParameters,
        const KnotSpanSkinBinsCSR& rSkinBinsPerKnotSpan,
        NurbsSurfaceType::Pointer pNurbsSurface,
        GeometryType::Pointer pSurrogateMiddleGeometry);

    /// Create the lateral interface conditions for one surrogate loop.
    void CreateLateralInterfaceConditions2D(
        ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart,
        IndexType FirstConditionId,
        IndexType LastConditionId,
        std::size_t NumberOfShapeFunctionsDerivatives,
        const Vector& rKnotSpanSizes);

    /**
     * @brief Builds the gap SBM geometries for a given skin/surrogate pair.
     * @tparam TIsInnerLoop True when processing the inner loop orientation.
     * @param rSkinSubModelPart Reference to the target skin sub model part.
     * @param rSurrogateSubModelPart Reference to the surrogate sub model part.
     */
    template <bool TIsInnerLoop>
    void CreateSbmExtendedGeometries3D(
        ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart);

    array_1d<double, 3> CalculateAxisAlignedBrepFaceNormal3D(
        const CoordinatesArrayType& rVertex00,
        const CoordinatesArrayType& rVertex01,
        const CoordinatesArrayType& rVertex10,
        const CoordinatesArrayType& rVertex11,
        const NodeType& rFirstSurrogateNode,
        const double MaxKnotSpanSize) const;

    /**
     * @brief Creates the projection surface between a surrogate edge and its projected skin edge.
     * @param pSurrogateNode1 First surrogate edge node.
     * @param pSurrogateNode2 Second surrogate edge node.
     * @param rSkinSubModelPart Skin model part storing the projected nodes.
     * @return Degree-one NURBS surface representing the current straight-edge Coons patch.
     */
    NurbsSurfaceType::Pointer CreateProjectionSurfaceFromSurrogateSegment(
        const Node::Pointer& pSurrogateNode1,
        const Node::Pointer& pSurrogateNode2,
        const ModelPart& rSkinSubModelPart) const;

    /**
     * @brief Creates the skin-side Coons surface bounded by the top edges of the lateral projection surfaces.
     * @param pProjectionSurface00_01 Lateral projection surface over edge pNode00-pNode01.
     * @param pProjectionSurface00_10 Lateral projection surface over edge pNode00-pNode10.
     * @param pProjectionSurface01_11 Lateral projection surface over edge pNode01-pNode11.
     * @param pProjectionSurface10_11 Lateral projection surface over edge pNode10-pNode11.
     * @param ReverseSurfaceNormal If true, swaps the surface parametric directions to reverse its normal.
     * @return Degree-one NURBS surface representing the current skin-side Coons approximation.
     */
    NurbsSurfaceType::Pointer CreateGapCoonsSurface(
        const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
        const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
        const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
        const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
        const bool ReverseSurfaceNormal = false) const;

    /**
     * @brief Creates the gap Coons volume bounded by the surrogate face, skin face, and lateral projection faces.
     * @param pNode00 Surrogate face node at (u,v,w) = (0,0,0).
     * @param pNode01 Surrogate face node at (u,v,w) = (0,1,0).
     * @param pNode10 Surrogate face node at (u,v,w) = (1,0,0).
     * @param pNode11 Surrogate face node at (u,v,w) = (1,1,0).
     * @param pProjectionSurface00_01 Lateral projection surface over edge pNode00-pNode01.
     * @param pProjectionSurface00_10 Lateral projection surface over edge pNode00-pNode10.
     * @param pProjectionSurface01_11 Lateral projection surface over edge pNode01-pNode11.
     * @param pProjectionSurface10_11 Lateral projection surface over edge pNode10-pNode11.
     * @param pTopSkinFace Top skin-side Coons surface representing the gap closure surface.
     * @return Degree-one NURBS volume representing the current Coons-volume approximation.
     */
    NurbsVolumeType::Pointer CreateGapCoonsVolume(
        const Node::Pointer& pNode00,
        const Node::Pointer& pNode01,
        const Node::Pointer& pNode10,
        const Node::Pointer& pNode11,
        const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
        const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
        const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
        const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
        const NurbsSurfaceType::Pointer& pTopSkinFace) const;

    /**
     * @brief Creates the gap skin conditions on the Coons skin surface and the corresponding gap volume elements.
     * @param pNode00 Surrogate face node at (u,v,w) = (0,0,0).
     * @param pNode01 Surrogate face node at (u,v,w) = (0,1,0).
     * @param pNode10 Surrogate face node at (u,v,w) = (1,0,0).
     * @param pNode11 Surrogate face node at (u,v,w) = (1,1,0).
     * @param pProjectionSurface00_01 Lateral projection surface over edge pNode00-pNode01.
     * @param pProjectionSurface00_10 Lateral projection surface over edge pNode00-pNode10.
     * @param pProjectionSurface01_11 Lateral projection surface over edge pNode01-pNode11.
     * @param pProjectionSurface10_11 Lateral projection surface over edge pNode10-pNode11.
     * @param rSurrogateBrepMiddleGeometry Surrogate reference geometry attached to the created entities.
     * @param rIntegrationParameters Integration settings and knot span metadata.
     * @param rSkinSubModelPart Skin model part containing the projected nodes.
     */
    void CreateGapElementsAndConditions(
        const Node::Pointer& pNode00,
        const Node::Pointer& pNode01,
        const Node::Pointer& pNode10,
        const Node::Pointer& pNode11,
        const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
        const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
        const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
        const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
        const GeometryType::Pointer& rSurrogateBrepMiddleGeometry,
        const IntegrationParameters& rIntegrationParameters,
        ModelPart& rSkinSubModelPart) const;

    /// Checks whether a node or condition belongs to the requested layer.
    static bool NodeHasRequestedLayer(const Node& rNode, const std::string& rLayer);
    static bool ConditionMatchesRequestedLayer(const Condition& rCondition, const std::string& rLayer);

    /// Searches for a skin node along one direction.
    std::pair<IndexType, double> SearchNodeWithDirection(
        const array_1d<double, 3>& rSearchDirection,
        const CoordinatesArrayType& rStartPoint,
        const std::string& rLayer,
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const double ReferenceSpanSize,
        const std::size_t BaseIndexX,
        const std::size_t BaseIndexY,
        const std::vector<IndexType>* pConsideredNodeIds) const;

    /// Searches for a skin segment along one direction.
    void SearchPairWithDirection(
        const array_1d<double, 3>& rSearchDirection,
        const CoordinatesArrayType& rStartPoint,
        const std::string& rLayer,
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const double ReferenceSpanSize,
        const std::size_t BaseIndexX,
        const std::size_t BaseIndexY,
        const bool ReverseNormalFilter,
        const bool IgnoreNormalFilter,
        std::pair<IndexType, IndexType>& rBestPair,
        double& rBestIntersectionDistance,
        std::vector<IndexType>& rClosestIntersectionCandidateIds) const;

    /// Selects the closest admissible skin node for a surrogate node.
    IndexType SelectSkinProjectionCandidate2D(
        const Node::Pointer& pSurrogateNode,
        const ModelPart& rSkinSubModelPart,
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const double SpanSizeX,
        const double SpanSizeY) const;

    /// Recursively projects subdivision points onto the skin.
    template <bool TIsInnerLoop>
    void ComputeRecursiveInterpolation(
        const std::size_t Depth,
        const std::size_t LeftIndex,
        const std::size_t RightIndex,
        const IndexType LeftNodeId,
        const IndexType RightNodeId,
        const ModelPart& rSkinSubModelPart,
        const std::string& rLayerName,
        const std::string& rConditionName,
        const IntegrationParameters& rIntegrationParameters,
        std::vector<IndexType>& rInterpolationProjectionIds) const;

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
        ModelPart& rSkinSubModelPart);

    /**
     * @brief Adds extra gap elements for curve intersections between surrogate1-skin1 and skin1-skin2.
     * @tparam TIsInnerLoop True when processing inner loop quadrature information.
     * @param rIntegrationParameters Reference to the integration parameters container.
     * @param pNurbsSurface Pointer to the NURBS surface used for projections.
     * @param rSurrogate1Skin1 Reference to the brep curve connecting surrogate_1 and skin_1.
     * @param rSkin1Skin2 Reference to the brep curve connecting skin_1 and skin_2.
     * @param rSurrogate1 Coordinates of the first surrogate point.
     * @param rSkin1 Coordinates of the first skin point.
     * @param rNeighbourGeometries Reference geometries used for the surrogate association.
     * @param rElementIdCounter Reference to the identifier counter used for the elements.
     */
    template <bool TIsInnerLoop>
    void CreateAdditionalGapElementsFromCurveIntersections(
        IntegrationParameters& rIntegrationParameters,
        NurbsSurfaceType::Pointer& pNurbsSurface,
        const BrepCurveType& rSurrogate1Skin1,
        const BrepCurveType& rSkin1Skin2,
        const array_1d<double,3>& rSurrogate1,
        const array_1d<double,3>& rSkin1,
        const std::vector<Geometry<Node>::Pointer>& rNeighbourGeometries,
        std::size_t& rElementIdCounter) const;

    /**
     * @brief Populates the surrogate-to-skin projection mapping.
     * @tparam TIsInnerLoop True when working with the inner loop orientation.
     * @param rSurrogateSubModelPart Reference to the surrogate sub model part.
     * @param rSkinSubModelPart Reference to the skin sub model part.
     * @param rSkinBinsPerSpan Reference to the skin CSR structure.
     */
    template <bool TIsInnerLoop>
    void SetSurrogateToSkinProjections2D(
        const ModelPart& rSurrogateSubModelPart,
        const ModelPart& rSkinSubModelPart,
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan);

    /**
     * @brief Populates the surrogate-to-skin projection mapping in 3D.
     * @tparam TIsInnerLoop True when working with the inner loop orientation.
     * @param rSurrogateSubModelPart Reference to the surrogate sub model part.
     * @param rSkinSubModelPart Reference to the skin sub model part.
     * @param rSkinBinsPerSpan Reference to the skin CSR structure.
     */
    template <bool TIsInnerLoop>
    void SetSurrogateToSkinProjections3D(
        const ModelPart& rSurrogateSubModelPart,
        ModelPart& rSkinSubModelPart,
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan);

    /**
     * @brief Checks that surrogate boundary vertices lie on the same skin layer.
     *        Otherwise search in the between the interface vertices for a surrogate vertex that does. 
     * @param rSkinSubModelPart Reference to the skin sub model part containing the layer information.
     * @param pSurrogateNode1 Pointer to the first surrogate vertex to be validated.
     * @param pSurrogateNode2 Pointer to the second surrogate vertex to be validated.
     */
    std::vector<IndexType> AssessProjectionsFeasibility(
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
        const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries,
        const double CharacteristicLength) const;

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
        const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries,
        const double CharacteristicLength) const;

    /**
     * @brief Compute Pojections of surrogate midpoint on the skin and determine the closest projection point and its id
     * 
     * @tparam TIsInnerLoop 
     * @param rSkinSubModelPart 
     * @param rSurrogateMidPoint 
     * @param LeftProjectionId 
     * @param RightProjectionId 
     * @param rCommonLayerName 
     * @param rNormalCond 
     * @param rIntegrationParameters 
     * @param rSkinConditionsPerKnotSpan 
     * @param rConsideredNodeIds Skin nodes already used by the current subdivision.
     * @return ProjectionResult 
     */
    template <bool TIsInnerLoop>
    ProjectionResult ComputeSingleProjection(
        const ModelPart& rSkinSubModelPart,
        const CoordinatesArrayType& rSurrogateMidPoint,
        const IndexType LeftProjectionId,
        const IndexType RightProjectionId,
        const std::string& rCommonLayerName,
        const std::string& rCommonConditionName,
        const array_1d<double,3>& rNormalCond,
        const IntegrationParameters& rIntegrationParameters,
        const KnotSpanSkinBinsCSR& rSkinBinsPerKnotSpan,
        const std::vector<IndexType>& rConsideredNodeIds) const;
    
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
        const Vector& rActiveRangeKnotVector);
    
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

    /**
     * @brief Returns Gauss points mapped on the gap Coons volume.
     * @param Order Quadrature order used in each parametric direction.
     * @param rGapVolume NURBS volume representing the current gap Coons volume.
     * @return Array of integration points with Jacobian-weighted volume weights.
     */
    IntegrationPointsArrayType CreateCoonsVolumeGaussPoints(
        std::size_t Order,
        const NurbsVolumeType& rGapVolume) const;

    // --- static utility helpers ----------------------------------

    /**
     * @brief Computes the global coordinates of a point on a Brep curve at the requested parameter.
     * @param rCurve Reference to the Brep curve geometry.
     * @param T Parametric coordinate where the curve is evaluated.
     * @return Global coordinates corresponding to the evaluated parameter.
     */
    static array_1d<double,3> EvaluateGlobalCurvePoint(
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
     * @param is_xi_derivative True to differentiate with respect to xi, false for eta.
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
        const bool                   is_xi_derivative,
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
     * @param is_xi_derivative True to differentiate with respect to xi, false for eta.
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
        const bool                       is_xi_derivative,
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
     * @brief Computes the characteristic length associated with the gap quadrilateral defined by the four reference points.
     * @param rSurrogatePoint1 Coordinates of the first surrogate node.
     * @param rSurrogatePoint2 Coordinates of the second surrogate node.
     * @param rSkinPoint1 Coordinates of the first skin node.
     * @param rSkinPoint2 Coordinates of the second skin node.
     * @return Radius obtained from the farthest vertex with respect to the centroid.
     */
    double CalculateGapElementCharacteristicLength(
        const array_1d<double,3>& rSurrogatePoint1,
        const array_1d<double,3>& rSurrogatePoint2,
        const array_1d<double,3>& rSkinPoint1,
        const array_1d<double,3>& rSkinPoint2) const;

    /**
     * @brief Computes a characteristic length for the 3D gap volume from its eight corner points.
     * @param rNode00 Surrogate face node at (u,v) = (0,0).
     * @param rNode01 Surrogate face node at (u,v) = (0,1).
     * @param rNode10 Surrogate face node at (u,v) = (1,0).
     * @param rNode11 Surrogate face node at (u,v) = (1,1).
     * @param rProjectedPoint00 Projected skin point of rNode00.
     * @param rProjectedPoint01 Projected skin point of rNode01.
     * @param rProjectedPoint10 Projected skin point of rNode10.
     * @param rProjectedPoint11 Projected skin point of rNode11.
     * @return Radius from the centroid to the farthest volume corner.
     */
    double CalculateGapVolumeElementCharacteristicLength(
        const NodeType& rNode00,
        const NodeType& rNode01,
        const NodeType& rNode10,
        const NodeType& rNode11,
        const array_1d<double,3>& rProjectedPoint00,
        const array_1d<double,3>& rProjectedPoint01,
        const array_1d<double,3>& rProjectedPoint10,
        const array_1d<double,3>& rProjectedPoint11) const;

    /**
     * @brief Reconstructs the ordered node path between two skin nodes by traversing skin conditions.
     * @param rSkinSubModelPart Skin model part containing the loop conditions.
     * @param StartNodeId Identifier of the starting node.
     * @param EndNodeId Identifier of the ending node.
     * @param ReverseConditionOrientation When true, traverses conditions against their stored orientation.
     * @param rLayerName Optional layer filter used to restrict the traversed conditions.
     * @param rConditionName Optional condition-name filter used to restrict the traversed conditions.
     * @return Ordered node ids from StartNodeId to EndNodeId, both included.
     */
    std::vector<IndexType> CollectOrderedSkinNodePath(
        const ModelPart& rSkinSubModelPart,
        IndexType StartNodeId,
        IndexType EndNodeId,
        bool ReverseConditionOrientation,
        const std::string& rLayerName = "",
        const std::string& rConditionName = "") const;


/**
 * @brief Collects skin nodes between two indices and projects them to the UV space.
 * @tparam TIsInnerLoop True when traversing the inner loop orientation.
 * @param rSkinSubModelPart Reference to the skin sub model part considered.
 * @param rSurface Reference to the surface used for the projection.
 * @param id_node_1 Identifier of the starting skin node.
 * @param id_node_2 Identifier of the ending skin node.
 * @return Ordered list of projected UV coordinates.
 */
//-------------------------------
template <bool TIsInnerLoop>
std::vector<array_1d<double,3>> CollectSkinUVBetween(
    const ModelPart& rSkinSubModelPart,
    const NurbsSurfaceType& rSurface,
    IndexType id_node_1,
    IndexType id_node_2,
    const std::string& rLayerName = "",
    const std::string& rConditionName = "") const;



/**
 * @brief Computes the chord-length parametrization on [0,1] for UV samples.
 * @param UV Collection of UV points describing the sampled curve.
 * @return Normalized cumulative arc-length parameters.
 */
//-------------------------------
std::vector<double> ChordLengthParams01_UV(const std::vector<array_1d<double,3>>& UV) const;



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
//-------------------------------
void BernsteinBasis(const int p, const double t, std::vector<double>& B) const;

/**
 * @brief Builds a Bézier knot vector with multiplicity p+1 at the boundaries.
 * @param p Polynomial degree of the Bézier curve.
 * @param rKnots Vector receiving the generated knot values.
 */
//-------------------------------
void BuildBezierKnots(const int p, Vector& rKnots) const;

/**
 * @brief Solves a symmetric positive definite linear system using a simple Cholesky factorization.
 * @param A Reference to the system matrix, overwritten with its factorization.
 * @param b Right-hand side vector, overwritten with the solution.
 * @return True when the system is successfully solved.
 */
//-------------------------------
bool SolveSPD(std::vector<std::vector<double>>& A, std::vector<double>& b) const;

/**
 * @brief Creates a UV Bézier curve of degree p from a set of control points.
 * @param CtrlsUV Control points expressed in UV coordinates.
 * @param p Polynomial degree of the resulting curve.
 * @return Pointer to the constructed Bézier curve geometry.
 */
//-------------------------------
typename NurbsCurveGeometryType::Pointer MakeBezierUV_FromControls(
    const std::vector<array_1d<double,3>>& CtrlsUV, const int p) const;

/**
 * @brief Performs a least-squares Bézier fit of degree p for UV samples with fixed endpoints.
 * @param UVsamples Ordered UV samples, inclusive of both endpoints.
 * @param p Polynomial degree of the fitted curve.
 * @param ridge Ridge regularization magnitude applied to the normal equations.
 * @return Pointer to the fitted Bézier curve geometry.
 */
//-------------------------------
typename NurbsCurveGeometryType::Pointer FitBezierUV_LS_Generic(
    const std::vector<array_1d<double,3>>& UVsamples, // ordered, inclusive
    const int p,
    const double ridge = GeometricTolerance) const;

/**
 * @brief Reverses the orientation of a UV Bézier curve.
 * @param p_forward Pointer to the curve whose orientation is reversed.
 * @return Pointer to a curve geometry with inverted control point order.
 */
//-------------------------------
typename NurbsCurveGeometryType::Pointer ReverseBezierUV_Generic(
    const typename NurbsCurveGeometryType::Pointer& p_forward) const;

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
//-------------------------------
template <bool TIsInnerLoop>
typename NurbsCurveGeometryType::Pointer FitUV_BetweenSkinNodes_Generic(
    const ModelPart& rSkinSubModelPart,
    const NurbsSurfaceType& rSurface,
    IndexType id_node_1,
    IndexType id_node_2,
    const int p,
    const double ridge = GeometricTolerance,
    const std::string& rLayerName = "",
    const std::string& rConditionName = "") const;

/**
 * @brief Finds the closest node belonging to the target layer by following a prescribed direction.
 * @tparam TIsInnerLoop True when searching along the inner loop orientation.
 * @param rStartPoint Starting point used for the directional search.
 * @param rLayer Name of the skin layer to be matched.
 * @param rSkinSubModelPart Reference to the skin sub model part containing the candidate nodes.
 * @param rKnotSpanSizes Reference to the knot span sizes guiding the search.
 * @param rSkinBinsPerSpan CSR matrix storing skin nodes and conditions per knot span.
 * @param rDirection Direction vector used for ray casting.
 * @param pConsideredNodeIds Nodes already used in the reconstruction, to prefer other candidates.
 * @return Identifier of the closest node fulfilling the criteria.
 */
template <bool TIsInnerLoop>
IndexType FindClosestNodeInLayerWithDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
    const Vector& rDirection,
    const std::vector<IndexType>* pConsideredNodeIds = nullptr) const;


/**
 * @brief 
 * 
 * @tparam TIsInnerLoop 
 * @param rStartPoint 
 * @param rLayer 
 * @param rSkinSubModelPart 
 * @param rKnotSpanSizes 
 * @param rSkinBinsPerSpan
 * @param rDirection 
 * @param pClosestIntersectionCandidateIds Endpoints of conditions tied at the selected intersection.
 * @param pConsideredNodeIds Nodes already used, to prefer an intersection with new candidates.
 * @param ReverseNormalFilter Reverse the skin-normal test for vertex projection corrections.
 * @return std::pair<SnakeGapSbmProcess::IndexType, SnakeGapSbmProcess::IndexType> 
 */
template <bool TIsInnerLoop>
std::pair<SnakeGapSbmProcess::IndexType, SnakeGapSbmProcess::IndexType> FindClosestPairInLayerWithNormalDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
    const Vector& rDirection,
    std::vector<IndexType>* pClosestIntersectionCandidateIds = nullptr,
    const std::vector<IndexType>* pConsideredNodeIds = nullptr,
    const bool ReverseNormalFilter = false,
    const bool IgnoreNormalFilter = false) const;
    
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



}; // Class SnakeGapSbmProcess

}  // namespace Kratos.
