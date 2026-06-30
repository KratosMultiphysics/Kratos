//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes
#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

// Project includes
#include "custom_processes/snake_gap_sbm_process.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) SnakeGapSbm3DUtilities
{
public:
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
    using VolumeQuadratureConsumerType = std::function<void(
        GeometriesArrayType&,
        const std::vector<Geometry<Node>::Pointer>&,
        double)>;
    
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

    using NodePointerType = Node::Pointer;
    using NeighbourGeometriesVectorType = std::vector<Geometry<Node>::Pointer>;
    using SkinEdgeKey = std::pair<IndexType, IndexType>;

    struct SkinEdgeKeyHash
    {
        std::size_t operator()(const SkinEdgeKey& rKey) const
        {
            const std::size_t first_hash = std::hash<IndexType>{}(rKey.first);
            const std::size_t second_hash = std::hash<IndexType>{}(rKey.second);
            return first_hash ^ (second_hash + 0x9e3779b9 + (first_hash << 6) + (first_hash >> 2));
        }
    };

    using SkinEdgeControlPointMap = std::unordered_map<
        SkinEdgeKey,
        std::vector<NodePointerType>,
        SkinEdgeKeyHash>;

    struct SkinTopFaceKey
    {
        std::array<IndexType, 3> NodeIds = {{0, 0, 0}};

        SkinTopFaceKey() = default;

        SkinTopFaceKey(
            const IndexType Id0,
            const IndexType Id1,
            const IndexType Id2)
        {
            NodeIds = {{Id0, Id1, Id2}};
            std::sort(NodeIds.begin(), NodeIds.end());
        }

        bool operator<(const SkinTopFaceKey& rOther) const
        {
            return NodeIds < rOther.NodeIds;
        }

        bool operator==(const SkinTopFaceKey& rOther) const
        {
            return NodeIds == rOther.NodeIds;
        }
    };

    struct SkinTopFaceKeyHash
    {
        std::size_t operator()(const SkinTopFaceKey& rKey) const
        {
            std::size_t seed = 0;
            for (const auto id : rKey.NodeIds) {
                const std::size_t value_hash = std::hash<IndexType>{}(id);
                seed ^= value_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    struct SpanKey3D
    {
        int I = 0;
        int J = 0;
        int K = 0;

        bool operator<(const SpanKey3D& rOther) const
        {
            return std::tie(I, J, K) < std::tie(rOther.I, rOther.J, rOther.K);
        }

        bool operator==(const SpanKey3D& rOther) const
        {
            return I == rOther.I && J == rOther.J && K == rOther.K;
        }
    };

    using ActiveSpanPredicateType = std::function<bool(const SpanKey3D&)>;

    enum class GapSpanType : int
    {
        Undefined = 0,
        Type1 = 1,
        Type2 = 2,
        Type3 = 3
    };

    struct GridPointKey3D
    {
        int I = 0;
        int J = 0;
        int K = 0;

        bool operator<(const GridPointKey3D& rOther) const
        {
            return std::tie(I, J, K) < std::tie(rOther.I, rOther.J, rOther.K);
        }

        bool operator==(const GridPointKey3D& rOther) const
        {
            return I == rOther.I && J == rOther.J && K == rOther.K;
        }
    };

    struct KnotSpanGridInfo
    {
        std::size_t NumberOfSpansX = 0;
        std::size_t NumberOfSpansY = 0;
        std::size_t NumberOfSpansZ = 0;

        double MinU = 0.0;
        double MaxU = 0.0;
        double MinV = 0.0;
        double MaxV = 0.0;
        double MinW = 0.0;
        double MaxW = 0.0;

        double SpanSizeX = 0.0;
        double SpanSizeY = 0.0;
        double SpanSizeZ = 0.0;
    };

    /**
     * @brief Configure for condition-based bins (AABB per condition geometry).
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

        static inline void CalculateBoundingBox(
            const PointerType& rObject,
            PointType& rLowPoint,
            PointType& rHighPoint)
        {
            const auto& r_geometry = rObject->GetGeometry();
            if (r_geometry.size() == 0) {
                rLowPoint = PointType(0.0, 0.0, 0.0);
                rHighPoint = rLowPoint;
                return;
            }
            rLowPoint = r_geometry[0];
            rHighPoint = r_geometry[0];
            for (IndexType i = 1; i < r_geometry.size(); ++i) {
                const auto& r_point = r_geometry[i];
                for (std::size_t d = 0; d < Dimension; ++d) {
                    rLowPoint[d] = std::min(rLowPoint[d], r_point[d]);
                    rHighPoint[d] = std::max(rHighPoint[d], r_point[d]);
                }
            }
        }

        static inline void CalculateBoundingBox(
            const PointerType& rObject,
            PointType& rLowPoint,
            PointType& rHighPoint,
            const double Radius)
        {
            CalculateBoundingBox(rObject, rLowPoint, rHighPoint);
            for (std::size_t d = 0; d < Dimension; ++d) {
                rLowPoint[d] -= Radius;
                rHighPoint[d] += Radius;
            }
        }

        static inline void CalculateCenter(
            const PointerType& rObject,
            PointType& rCentralPoint)
        {
            rCentralPoint = rObject->GetGeometry().Center();
        }

        static inline bool Intersection(
            const PointerType& rObj_1,
            const PointerType& rObj_2)
        {
            PointType low_1, high_1, low_2, high_2;
            CalculateBoundingBox(rObj_1, low_1, high_1);
            CalculateBoundingBox(rObj_2, low_2, high_2);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high_1[d] < low_2[d] - Epsilon || low_1[d] > high_2[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

        static inline bool Intersection(
            const PointerType& rObj_1,
            const PointerType& rObj_2,
            const double Radius)
        {
            PointType low_1, high_1, low_2, high_2;
            CalculateBoundingBox(rObj_1, low_1, high_1, Radius);
            CalculateBoundingBox(rObj_2, low_2, high_2, Radius);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high_1[d] < low_2[d] - Epsilon || low_1[d] > high_2[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

        static inline bool IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint)
        {
            PointType low, high;
            CalculateBoundingBox(rObject, low, high);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high[d] < rLowPoint[d] - Epsilon || low[d] > rHighPoint[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

        static inline bool IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint,
            const double Radius)
        {
            PointType low, high;
            CalculateBoundingBox(rObject, low, high, Radius);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high[d] < rLowPoint[d] - Epsilon || low[d] > rHighPoint[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

        static inline void Distance(
            const PointerType& rObj_1,
            const PointerType& rObj_2,
            double& distance)
        {
            PointType c1, c2;
            CalculateCenter(rObj_1, c1);
            CalculateCenter(rObj_2, c2);
            double pwdDistance = 0.0;
            for (std::size_t d = 0; d < Dimension; ++d) {
                const double delta = c1[d] - c2[d];
                pwdDistance += delta * delta;
            }
            distance = std::sqrt(pwdDistance);
        }

        static inline double GetObjectRadius(
            const PointerType&,
            const double)
        {
            return 0.0;
        }
    };

    struct KnotSpanSkinCellData
    {
        std::vector<NodePointerType> Nodes;
        std::vector<ConditionPointerType> Conditions;
        std::unique_ptr<NodeBinsType> pNodeBins;
        BinsObjectDynamic<ConditionConfigure> ConditionBins;
        bool HasNodeBins = false;
        bool HasConditionBins = false;
    };

    struct KnotSpanSkinBinsCSR
    {
        KnotSpanGridInfo GridInfo;
        CompressedMatrix Occupancy;
        std::vector<KnotSpanSkinCellData> CellDataByNnz;
    };

    struct SurrogateFaceData
    {
        const Condition* pCondition = nullptr;

        std::array<NodePointerType, 4> Nodes;
        std::array<GridPointKey3D, 4> GridNodes;

        SpanKey3D ActiveSpan;
        SpanKey3D ExternalSpan;

        int NormalAxis = -1;
        int NormalSign = 0;

        array_1d<double, 3> Normal;
    };

    struct ExternalSpanData
    {
        SpanKey3D Key;
        GapSpanType Type = GapSpanType::Undefined;

        IndexType ProjectionNodeId = 0;
        NodePointerType pProjectionNode;

        std::vector<IndexType> AdjacentSurrogateConditionIds;
        std::vector<SpanKey3D> AdjacentActiveSpans;

        bool HasProjectionNode() const
        {
            return ProjectionNodeId != 0 && pProjectionNode != nullptr;
        }
    };

    //--------------------------------------------
    struct CanonicalFaceKey3D
    {
        std::array<const NodeType*, 3> NodePointers = {{nullptr, nullptr, nullptr}};
        std::array<IndexType, 3> NodeIds = {{0, 0, 0}};

        bool operator==(const CanonicalFaceKey3D& rOther) const
        {
            if (NodePointers[0] != nullptr ||
                rOther.NodePointers[0] != nullptr) {
                return NodePointers == rOther.NodePointers;
            }

            return NodeIds == rOther.NodeIds;
        }
    };

    struct CanonicalFaceKey3DHasher
    {
        std::size_t operator()(const CanonicalFaceKey3D& rKey) const
        {
            std::size_t seed = 0;

            if (rKey.NodePointers[0] != nullptr) {
                for (const auto p_node : rKey.NodePointers) {
                    seed ^= std::hash<const NodeType*>{}(p_node)
                        + 0x9e3779b9
                        + (seed << 6)
                        + (seed >> 2);
                }
            } else {
                for (const auto id : rKey.NodeIds) {
                    seed ^= std::hash<IndexType>{}(id)
                        + 0x9e3779b9
                        + (seed << 6)
                        + (seed >> 2);
                }
            }

            return seed;
        }
    };

    struct LateralFaceOccurrence
    {
        int GapType = 0;
        IndexType SurrogateConditionId = 0;
        SpanKey3D ExternalSpan;

        GeometryType::Pointer pGeometry;
        Geometry<Node>::Pointer pNeighbourGeometry;
        NodePointerType pOppositeNode;
        std::array<NodePointerType, 3> FaceNodes;

        bool HasType1NeighbourPath = false;
        bool HasNeighbourActiveSpan = false;
        SpanKey3D NeighbourActiveSpan;
    };

    struct CurvedEdgeData
    {
        SkinEdgeKey Key;
        std::vector<NodePointerType> LinearControlNodes;
        std::vector<NodePointerType> CurvedControlNodes;
        std::vector<NodePointerType> CurrentControlNodes;
        NodePointerType pQuadraticInterpolationNode;
        double Alpha = 0.0;
        bool IsCurved = false;
    };

    struct CurvedTopFaceData
    {
        SkinTopFaceKey Key;
        std::array<NodePointerType, 3> CornerNodes;
        std::vector<NodePointerType> CurrentControlNodes;
        GeometryType::Pointer pSurfaceGeometry;
        double Alpha = 0.0;
        bool IsCurved = false;
    };

    using CurvedEdgeRegistry = std::unordered_map<
        SkinEdgeKey,
        CurvedEdgeData,
        SkinEdgeKeyHash>;
    using CurvedTopFaceRegistry = std::unordered_map<
        SkinTopFaceKey,
        CurvedTopFaceData,
        SkinTopFaceKeyHash>;

    using LateralFaceRegistry = std::unordered_map<
        CanonicalFaceKey3D,
        std::vector<LateralFaceOccurrence>,
        CanonicalFaceKey3DHasher>;

    struct BrepPatchData
    {
        GeometryType::Pointer pBrepGeometry;
        BrepSurfaceOnVolumeType::Pointer pBrepSurface;

        array_1d<double, 3> Vertex00 = ZeroVector(3);
        array_1d<double, 3> Vertex01 = ZeroVector(3);
        array_1d<double, 3> Vertex10 = ZeroVector(3);
        array_1d<double, 3> Vertex11 = ZeroVector(3);

        array_1d<double, 3> MiddlePoint = ZeroVector(3);
        array_1d<double, 3> MiddlePointLocalCoordinates = ZeroVector(3);
    };

    struct Type1VolumeQuadratureData
    {
        GeometriesArrayType VolumeQuadraturePointGeometries;
        std::vector<Geometry<Node>::Pointer> NeighbourGeometries;

        double CharacteristicLength = 0.0;

        std::array<NodePointerType, 4> BaseNodes;
        NodePointerType pProjectionNode;

        IndexType ProjectionNodeId = 0;
        IndexType SurrogateConditionId = 0;
        SpanKey3D ExternalSpan;
    };

    struct Type1CreationSummary
    {
        std::size_t NumberOfPyramids = 0;
        std::size_t NumberOfLateralFaces = 0;
        std::size_t NumberOfOpenFaces = 0;
        std::size_t NumberOfClosedFaces = 0;
        std::size_t NumberOfNonManifoldFaces = 0;
    };



    struct LateralSurfaceQuadratureData
    {
        GeometriesArrayType SurfaceQuadraturePointGeometries;
        std::vector<Geometry<Node>::Pointer> NeighbourGeometries;

        double CharacteristicLength = 0.0;

        int GapType = 0;
        IndexType SurrogateConditionId = 0;
        SpanKey3D ExternalSpan;
    };
    //--------------------------------------------
    // TYPE 2 ELEMENTS
    struct SurrogateEdgeKey3D
    {
        std::array<IndexType, 2> NodeIds = {0, 0};

        SurrogateEdgeKey3D() = default;

        SurrogateEdgeKey3D(
            const IndexType Id0,
            const IndexType Id1)
        {
            NodeIds[0] = std::min(Id0, Id1);
            NodeIds[1] = std::max(Id0, Id1);
        }

        bool operator==(const SurrogateEdgeKey3D& rOther) const
        {
            return NodeIds == rOther.NodeIds;
        }
    };

    struct SurrogateEdgeKey3DHasher
    {
        std::size_t operator()(const SurrogateEdgeKey3D& rKey) const
        {
            const std::size_t h0 = std::hash<IndexType>{}(rKey.NodeIds[0]);
            const std::size_t h1 = std::hash<IndexType>{}(rKey.NodeIds[1]);

            return h0 ^ (h1 + 0x9e3779b9 + (h0 << 6) + (h0 >> 2));
        }
    };
    struct Type1LateralFaceData
    {
        SurrogateEdgeKey3D EdgeKey;

        NodePointerType pEdgeNode0;
        NodePointerType pEdgeNode1;
        NodePointerType pProjectionNode;

        GridPointKey3D EdgeGridNode0;
        GridPointKey3D EdgeGridNode1;

        NurbsSurfaceType::Pointer pSurfaceGeometry;
        Geometry<Node>::Pointer pNeighbourGeometry;

        IndexType SurrogateConditionId = 0;
        SpanKey3D ActiveSpan;
        SpanKey3D ExternalSpan;

        bool IsClosedByType2 = false;
    };

    using Type1LateralFaceContainer = std::vector<Type1LateralFaceData>;

    using Type1LateralFacesByEdgeMap = std::unordered_map<
        SurrogateEdgeKey3D,
        std::vector<std::size_t>,
        SurrogateEdgeKey3DHasher>;

    using ExternalSpanDataMap = std::map<SpanKey3D, ExternalSpanData>;

    struct Type1CreationResult
    {
        Type1CreationSummary Summary;
        std::vector<Type1VolumeQuadratureData> VolumeQuadratureDataList;

        Type1LateralFaceContainer Type1LateralFaces;
        Type1LateralFacesByEdgeMap Type1LateralFacesByEdge;
    };

    //-------------------------------------------------
    // type 2 struct
    struct Type2VolumeQuadratureData
    {
        GeometriesArrayType VolumeQuadraturePointGeometries;
        std::vector<Geometry<Node>::Pointer> NeighbourGeometries;

        double CharacteristicLength = 0.0;

        SurrogateEdgeKey3D EdgeKey;

        NodePointerType pEdgeNode0;
        NodePointerType pEdgeNode1;
        NodePointerType pProjectionNode0;
        NodePointerType pProjectionNode1;

        IndexType ProjectionNodeId0 = 0;
        IndexType ProjectionNodeId1 = 0;

        SpanKey3D ExternalSpan0;
        SpanKey3D ExternalSpan1;

        IndexType FirstType1LateralFaceIndex = 0;
        IndexType SecondType1LateralFaceIndex = 0;

        bool RequiresQuadrature = false;
    };

    struct Type2OpenFaceData
    {
        NodePointerType pSurrogateNode;
        NodePointerType pProjectionNode0;
        NodePointerType pProjectionNode1;
        NodePointerType pOppositeNode;

        GeometryType::Pointer pSurfaceGeometry;

        std::vector<Geometry<Node>::Pointer> NeighbourGeometries;

        SurrogateEdgeKey3D ParentEdgeKey;
        SpanKey3D ExternalSpan0;
        SpanKey3D ExternalSpan1;
    };

    struct Type2CreationSummary
    {
        std::size_t NumberOfCandidateEdges = 0;
        std::size_t NumberOfCreatedVolumes = 0;
        std::size_t NumberOfSkippedEdges = 0;
        std::size_t NumberOfOpenFaces = 0;
    };

    struct Type2CreationResult
    {
        Type2CreationSummary Summary;
        std::vector<Type2VolumeQuadratureData> VolumeQuadratureDataList;
        std::vector<Type2OpenFaceData> OpenFaceDataList;
    };
    //-------------------------------------------------

    //-------------------------------------------------
    // type 3 structs
    struct Type3VolumeQuadratureData
    {
        GeometriesArrayType VolumeQuadraturePointGeometries;
        std::vector<Geometry<Node>::Pointer> NeighbourGeometries;

        double CharacteristicLength = 0.0;

        NodePointerType pSurrogateNode;
        NodePointerType pProjectionNode0;
        NodePointerType pProjectionNode1;
        NodePointerType pProjectionNode2;

        IndexType SurrogateNodeId = 0;
        IndexType ProjectionNodeId0 = 0;
        IndexType ProjectionNodeId1 = 0;
        IndexType ProjectionNodeId2 = 0;
        IndexType CornerProjectionNodeId = 0;

        bool RequiresQuadrature = false;
    };

    struct Type3CreationSummary
    {
        std::size_t NumberOfCandidateOpenFaces = 0;
        std::size_t NumberOfCreatedVolumes = 0;
        std::size_t NumberOfSkippedFaces = 0;
        std::size_t NumberOfCreatedCornerProjectionNodes = 0;
        std::size_t NumberOfReusedCornerProjectionNodes = 0;
    };

    struct Type3CreationResult
    {
        Type3CreationSummary Summary;
        std::vector<Type3VolumeQuadratureData> VolumeQuadratureDataList;
    };
    //-------------------------------------------------

    explicit SnakeGapSbm3DUtilities(const int EchoLevel = 0);

    KnotSpanGridInfo CreateKnotSpanGridInfo(
        const ModelPart& rBackgroundModelPart) const;

    KnotSpanSkinBinsCSR BuildSkinBinsPerKnotSpan3D(
        const ModelPart& rSkinSubModelPart,
        const KnotSpanGridInfo& rGridInfo) const;

    std::vector<SurrogateFaceData> BuildSurrogateFaceDataVector(
        const ModelPart& rSurrogateSubModelPart,
        const KnotSpanGridInfo& rGridInfo) const;

    std::set<SpanKey3D> ExtractActiveSpans(
        const std::vector<SurrogateFaceData>& rSurrogateFaces) const;

    ExternalSpanDataMap ClassifyExternalSpans(
        const std::set<SpanKey3D>& rActiveSpans,
        const std::vector<SurrogateFaceData>& rSurrogateFaces,
        const KnotSpanGridInfo& rGridInfo,
        const ActiveSpanPredicateType& rIsSpanActive) const;

    void AssignProjectionNodesToExternalSpans(
        ExternalSpanDataMap& rExternalSpans,
        ModelPart& rSkinSubModelPart,
        const KnotSpanSkinBinsCSR& rSkinBins,
        const KnotSpanGridInfo& rGridInfo) const;

    void AddInteriorSkinNodesForLargeTriangles(
        ModelPart& rSkinSubModelPart,
        const KnotSpanGridInfo& rGridInfo) const;

    ExternalSpanDataMap InitializeExternalSpanData(
        ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart,
        const ActiveSpanPredicateType& rIsSpanActive) const;

    std::string SpanToString(const SpanKey3D& rSpan) const;

    /**
     * @brief Builds a 3D CSR matrix linking skin nodes and conditions to knot span bins.
     * @param rSkinSubModelPart Reference to the skin model part.
     * @param rSurrogateSubModelPart Reference to the surrogate model part used for binning.
     * @return CSR structure collecting node and condition bins per knot span.
     */
    KnotSpanSkinBinsCSR BuildSkinBinsPerKnotSpan3D(
        const ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart) const;

    void SetGapApproximationOrder(
        const std::size_t GapApproximationOrder);

    void SetStoreGapDebugGeometries(
        const bool StoreGapDebugGeometries);

    void SetUsePyramidQuadratureForType1(
        const bool UsePyramidQuadratureForType1);

    void SetUseTetraQuadratureForLinearType2Type3(
        const bool UseTetraQuadratureForLinearType2Type3);

    void SetUseAnisotropicQuadratureForCurvedType2(
        const bool UseAnisotropicQuadratureForCurvedType2);

    void SetUseTriangleRadialQuadratureForCurvedType3(
        const bool UseTriangleRadialQuadratureForCurvedType3);

    
    // type 1 element creation----------------------------
    void ClearLateralFaceRegistry();

    const LateralFaceRegistry& GetLateralFaceRegistry() const;
    
    Type1CreationResult CreateType1GapGeometries(
        ModelPart& rRootModelPart,
        const ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart,
        const ExternalSpanDataMap& rExternalSpans,
        const std::size_t GapVolumeIntegrationOrder);

    
    std::vector<LateralSurfaceQuadratureData> CreateOpenLateralSurfaceQuadratureData(
        const std::size_t IntegrationOrder,
        const std::size_t NumberOfShapeFunctionsDerivatives) const;
    
    std::vector<LateralSurfaceQuadratureData> CreateInterfaceLateralSurfaceQuadratureData(
        const std::size_t IntegrationOrder,
        const std::size_t NumberOfShapeFunctionsDerivatives) const;
    //-----------------------------------------------------

    // type 2 element creation----------------------------
    Type2CreationResult CreateType2GapGeometries(
        ModelPart& rRootModelPart,
        ModelPart& rSkinSubModelPart,
        Type1CreationResult& rType1CreationResult,
        const ExternalSpanDataMap& rExternalSpans,
        const KnotSpanGridInfo& rGridInfo,
        const std::size_t IntegrationOrder,
        const std::size_t NumberOfShapeFunctionsDerivatives);

    Type3CreationResult CreateType3GapGeometries(
        ModelPart& rRootModelPart,
        ModelPart& rSkinSubModelPart,
        const ExternalSpanDataMap& rExternalSpans,
        const KnotSpanGridInfo& rGridInfo,
        const Type2CreationResult& rType2CreationResult,
        const std::size_t IntegrationOrder,
        const std::size_t NumberOfShapeFunctionsDerivatives);

    void FinalizeType2AndType3GapGeometries(
        ModelPart& rRootModelPart,
        ModelPart& rSkinSubModelPart,
        Type2CreationResult& rType2CreationResult,
        Type3CreationResult& rType3CreationResult,
        const KnotSpanGridInfo& rGridInfo,
        const std::size_t IntegrationOrder,
        const std::size_t NumberOfShapeFunctionsDerivatives,
        const VolumeQuadratureConsumerType& rVolumeQuadratureConsumer =
            VolumeQuadratureConsumerType());

    NurbsSurfaceType::Pointer CreateCollapsedTriangleSurface(
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2) const;

    GeometryType::Pointer GetOrCreateLateralFaceSurface(
        ModelPart& rDebugModelPart,
        IndexType& rNextGeometryId,
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const Geometry<Node>::Pointer& pNeighbourGeometry);

    void RegisterLateralFaceOccurrence(
        const int GapType,
        const IndexType SurrogateConditionId,
        const SpanKey3D& rExternalSpan,
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const NodePointerType& pOppositeNode,
        const GeometryType::Pointer& pSurfaceGeometry,
        const Geometry<Node>::Pointer& pNeighbourGeometry,
        const bool HasType1NeighbourPath,
        const bool HasNeighbourActiveSpan,
        const SpanKey3D& rNeighbourActiveSpan);

private:
    std::size_t ComputeSpanCount(
        const double DomainLength,
        const double SpanSize,
        const char* pDirection,
        const char* pCaller) const;

    std::size_t ComputeSpanIndex(
        const double Coordinate,
        const double MinValue,
        const double MaxValue,
        const double SpanSize,
        const std::size_t NumberOfSpans,
        const double Tolerance,
        const char* pCaller) const;

    int ComputeGridLineIndex(
        const double Coordinate,
        const double MinValue,
        const double SpanSize,
        const std::size_t NumberOfSpans) const;

    std::vector<std::size_t> ComputeNodeSpanIndices(
        const double Coordinate,
        const double MinValue,
        const double MaxValue,
        const double SpanSize,
        const std::size_t NumberOfSpans,
        const double Tolerance,
        const double BoundaryTolerance,
        const char* pCaller) const;

    std::pair<std::size_t, std::size_t> ComputeSpanIndexRange(
        const double MinCoordinate,
        const double MaxCoordinate,
        const double MinValue,
        const double MaxValue,
        const double SpanSize,
        const std::size_t NumberOfSpans,
        const double Tolerance,
        const double BoundaryTolerance,
        const char* pCaller) const;

    std::size_t FlattenSpanColumnIndex(
        const KnotSpanGridInfo& rGridInfo,
        const std::size_t J,
        const std::size_t K) const;

    std::size_t FindSpanNnzIndex(
        const KnotSpanSkinBinsCSR& rSkinBins,
        const SpanKey3D& rSpan) const;

    bool IsSpanInsideDomain(
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo) const;

    array_1d<double, 3> SpanCenter(
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo) const;

    void ConditionBoundingBox(
        const Geometry<Node>& rGeometry,
        array_1d<double, 3>& rMinPoint,
        array_1d<double, 3>& rMaxPoint) const;

    GridPointKey3D ComputeGridPointKey(
        const array_1d<double, 3>& rPoint,
        const KnotSpanGridInfo& rGridInfo) const;

    SurrogateFaceData CreateSurrogateFaceData(
        const Condition& rCondition,
        const KnotSpanGridInfo& rGridInfo) const;

    int ComputeFaceNormalAxis(
        const array_1d<double, 3>& rNormal) const;

    int ComputeFaceNormalSign(
        const array_1d<double, 3>& rNormal,
        const int NormalAxis) const;

    SpanKey3D ComputeExternalSpanFromSurrogateFace(
        const SurrogateFaceData& rFaceData,
        const KnotSpanGridInfo& rGridInfo) const;

    SpanKey3D ComputeActiveSpanFromSurrogateFace(
        const SurrogateFaceData& rFaceData) const;

    void AddUniqueActiveSpan(
        std::vector<SpanKey3D>& rSpans,
        const SpanKey3D& rSpan) const;

    void AddUniqueConditionId(
        std::vector<IndexType>& rConditionIds,
        const IndexType ConditionId) const;

    void RegisterExternalSpanCandidate(
        ExternalSpanDataMap& rExternalSpans,
        const SpanKey3D& rExternalSpan,
        const GapSpanType Type,
        const SpanKey3D& rAdjacentActiveSpan,
        const IndexType AdjacentSurrogateConditionId) const;

    NodePointerType FindOrCreateProjectionNodeInSpan(
        ModelPart& rSkinSubModelPart,
        const SpanKey3D& rSpan,
        const array_1d<double, 3>& rReferencePoint,
        const KnotSpanSkinBinsCSR& rSkinBins,
        const KnotSpanGridInfo& rGridInfo) const;

    void PrintExternalSpanSummary(
        const ExternalSpanDataMap& rExternalSpans) const;

    
    bool FindTriangleBoxIntersectionPoint(
        const Geometry<Node>& rGeometry,
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo,
        array_1d<double, 3>& rIntersectionPoint) const
    {
        if (!IsSpanInsideDomain(rSpan, rGridInfo) ||
            rGeometry.PointsNumber() < 3) {
            return false;
        }

        array_1d<double, 3> span_min = ZeroVector(3);
        array_1d<double, 3> span_max = ZeroVector(3);
    
        span_min[0] = rGridInfo.MinU + static_cast<double>(rSpan.I) * rGridInfo.SpanSizeX;
        span_max[0] = span_min[0] + rGridInfo.SpanSizeX;
    
        span_min[1] = rGridInfo.MinV + static_cast<double>(rSpan.J) * rGridInfo.SpanSizeY;
        span_max[1] = span_min[1] + rGridInfo.SpanSizeY;
    
        span_min[2] = rGridInfo.MinW + static_cast<double>(rSpan.K) * rGridInfo.SpanSizeZ;
        span_max[2] = span_min[2] + rGridInfo.SpanSizeZ;
    
        constexpr double tolerance = 1.0e-12;

        auto clip_triangle_to_box = [&](
            const array_1d<double, 3>& rPoint0,
            const array_1d<double, 3>& rPoint1,
            const array_1d<double, 3>& rPoint2,
            array_1d<double, 3>& rPointInsideBox) -> bool
        {
            std::vector<array_1d<double, 3>> polygon = {
                rPoint0,
                rPoint1,
                rPoint2};

            auto clip_with_plane = [&](
                const IndexType Axis,
                const double Bound,
                const bool KeepGreater)
            {
                if (polygon.empty()) {
                    return;
                }

                std::vector<array_1d<double, 3>> clipped_polygon;
                clipped_polygon.reserve(polygon.size() + 1);

                auto is_inside = [&](const array_1d<double, 3>& rPoint) {
                    if (KeepGreater) {
                        return rPoint[Axis] >= Bound - tolerance;
                    }

                    return rPoint[Axis] <= Bound + tolerance;
                };

                auto intersection_point = [&](
                    const array_1d<double, 3>& rStart,
                    const array_1d<double, 3>& rEnd) -> array_1d<double, 3>
                {
                    const double denominator = rEnd[Axis] - rStart[Axis];
                    if (std::abs(denominator) <= tolerance) {
                        return rStart;
                    }

                    const double t = std::max(
                        0.0,
                        std::min(1.0, (Bound - rStart[Axis]) / denominator));

                    array_1d<double, 3> point = rStart;
                    noalias(point) += t * (rEnd - rStart);
                    return point;
                };

                for (std::size_t i = 0; i < polygon.size(); ++i) {
                    const auto& r_current_point = polygon[i];
                    const auto& r_next_point =
                        polygon[(i + 1) % polygon.size()];

                    const bool current_is_inside =
                        is_inside(r_current_point);
                    const bool next_is_inside =
                        is_inside(r_next_point);

                    if (current_is_inside && next_is_inside) {
                        clipped_polygon.push_back(r_next_point);
                    } else if (current_is_inside && !next_is_inside) {
                        clipped_polygon.push_back(intersection_point(
                            r_current_point,
                            r_next_point));
                    } else if (!current_is_inside && next_is_inside) {
                        clipped_polygon.push_back(intersection_point(
                            r_current_point,
                            r_next_point));
                        clipped_polygon.push_back(r_next_point);
                    }
                }

                polygon = std::move(clipped_polygon);
            };

            for (IndexType axis = 0; axis < 3; ++axis) {
                clip_with_plane(axis, span_min[axis], true);
                clip_with_plane(axis, span_max[axis], false);
            }

            if (polygon.empty()) {
                return false;
            }

            rPointInsideBox = ZeroVector(3);
            for (const auto& r_polygon_point : polygon) {
                noalias(rPointInsideBox) += r_polygon_point;
            }
            rPointInsideBox /= static_cast<double>(polygon.size());

            return IsPointInsideBox(
                rPointInsideBox,
                span_min,
                span_max,
                tolerance);
        };

        const array_1d<double, 3> point_0 = rGeometry[0].Coordinates();
        const array_1d<double, 3> point_1 = rGeometry[1].Coordinates();
        const array_1d<double, 3> point_2 = rGeometry[2].Coordinates();

        if (clip_triangle_to_box(
                point_0,
                point_1,
                point_2,
                rIntersectionPoint)) {
            return true;
        }

        if (rGeometry.PointsNumber() == 4) {
            const array_1d<double, 3> point_3 =
                rGeometry[3].Coordinates();

            if (clip_triangle_to_box(
                    point_0,
                    point_2,
                    point_3,
                    rIntersectionPoint)) {
                return true;
            }
        }

        return false;
    }

    void ClampPointInsideSpan(
        array_1d<double, 3>& rPoint,
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo) const
    {
        array_1d<double, 3> span_min = ZeroVector(3);
        array_1d<double, 3> span_max = ZeroVector(3);
    
        span_min[0] = rGridInfo.MinU + static_cast<double>(rSpan.I) * rGridInfo.SpanSizeX;
        span_max[0] = span_min[0] + rGridInfo.SpanSizeX;
    
        span_min[1] = rGridInfo.MinV + static_cast<double>(rSpan.J) * rGridInfo.SpanSizeY;
        span_max[1] = span_min[1] + rGridInfo.SpanSizeY;
    
        span_min[2] = rGridInfo.MinW + static_cast<double>(rSpan.K) * rGridInfo.SpanSizeZ;
        span_max[2] = span_min[2] + rGridInfo.SpanSizeZ;
    
        constexpr double relative_margin = 1.0e-12;
    
        const double margin_x = relative_margin * rGridInfo.SpanSizeX;
        const double margin_y = relative_margin * rGridInfo.SpanSizeY;
        const double margin_z = relative_margin * rGridInfo.SpanSizeZ;
    
        rPoint[0] = std::max(span_min[0] + margin_x, std::min(span_max[0] - margin_x, rPoint[0]));
        rPoint[1] = std::max(span_min[1] + margin_y, std::min(span_max[1] - margin_y, rPoint[1]));
        rPoint[2] = std::max(span_min[2] + margin_z, std::min(span_max[2] - margin_z, rPoint[2]));
    }

    std::size_t GetNextAuxiliarySkinNodeId(
        const ModelPart& rSkinSubModelPart) const
    {
        IndexType next_node_id = 1;
    
        const ModelPart& r_root_model_part = rSkinSubModelPart.GetRootModelPart();
    
        for (const auto& r_node : r_root_model_part.Nodes()) {
            next_node_id = std::max(next_node_id, r_node.Id() + 1);
        }
    
        return next_node_id;
    }

    // -----------------------------------------------------------------
    void SpanBox(
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo,
        array_1d<double, 3>& rBoxMin,
        array_1d<double, 3>& rBoxMax) const
    {
        KRATOS_ERROR_IF_NOT(IsSpanInsideDomain(rSpan, rGridInfo))
            << "[SnakeGapSbm3DUtilities::SpanBox] Span "
            << SpanToString(rSpan) << " is outside the knot-span grid.\n";
    
        rBoxMin[0] = rGridInfo.MinU + static_cast<double>(rSpan.I) * rGridInfo.SpanSizeX;
        rBoxMax[0] = rBoxMin[0] + rGridInfo.SpanSizeX;
    
        rBoxMin[1] = rGridInfo.MinV + static_cast<double>(rSpan.J) * rGridInfo.SpanSizeY;
        rBoxMax[1] = rBoxMin[1] + rGridInfo.SpanSizeY;
    
        rBoxMin[2] = rGridInfo.MinW + static_cast<double>(rSpan.K) * rGridInfo.SpanSizeZ;
        rBoxMax[2] = rBoxMin[2] + rGridInfo.SpanSizeZ;
    }

    bool IsPointInsideBox(
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rBoxMin,
        const array_1d<double, 3>& rBoxMax,
        const double Tolerance) const
    {
        return rPoint[0] >= rBoxMin[0] - Tolerance && rPoint[0] <= rBoxMax[0] + Tolerance &&
               rPoint[1] >= rBoxMin[1] - Tolerance && rPoint[1] <= rBoxMax[1] + Tolerance &&
               rPoint[2] >= rBoxMin[2] - Tolerance && rPoint[2] <= rBoxMax[2] + Tolerance;
    }

    bool SpanHasSkinNode(
        const ModelPart& rSkinSubModelPart,
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo) const
    {
        if (!IsSpanInsideDomain(rSpan, rGridInfo)) {
            return false;
        }
    
        array_1d<double, 3> box_min = ZeroVector(3);
        array_1d<double, 3> box_max = ZeroVector(3);
        SpanBox(rSpan, rGridInfo, box_min, box_max);
    
        const double tolerance =
            1.0e-10 * std::max({
                rGridInfo.SpanSizeX,
                rGridInfo.SpanSizeY,
                rGridInfo.SpanSizeZ});
    
        for (const auto& r_node : rSkinSubModelPart.Nodes()) {
            if (IsPointInsideBox(
                    r_node.Coordinates(),
                    box_min,
                    box_max,
                    tolerance)) {
                return true;
            }
        }
    
        return false;
    }

    bool FindConditionSpanIntersectionPoint(
        const Geometry<Node>& rGeometry,
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo,
        array_1d<double, 3>& rPointInsideSpan) const
    {
        if (!IsSpanInsideDomain(rSpan, rGridInfo)) {
            return false;
        }
    
        if (rGeometry.PointsNumber() < 3) {
            return false;
        }
    
        array_1d<double, 3> box_min = ZeroVector(3);
        array_1d<double, 3> box_max = ZeroVector(3);
        SpanBox(rSpan, rGridInfo, box_min, box_max);
    
        const double tolerance =
            1.0e-10 * std::max({
                rGridInfo.SpanSizeX,
                rGridInfo.SpanSizeY,
                rGridInfo.SpanSizeZ});
    
        auto segment_box_intersection_point = [&](
            const array_1d<double, 3>& rA,
            const array_1d<double, 3>& rB,
            array_1d<double, 3>& rPointInsideBox) -> bool
        {
            double t_min = 0.0;
            double t_max = 1.0;
    
            const array_1d<double, 3> direction = rB - rA;
    
            for (IndexType axis = 0; axis < 3; ++axis) {
                if (std::abs(direction[axis]) <= tolerance) {
                    if (rA[axis] < box_min[axis] - tolerance ||
                        rA[axis] > box_max[axis] + tolerance) {
                        return false;
                    }
                } else {
                    double t1 = (box_min[axis] - rA[axis]) / direction[axis];
                    double t2 = (box_max[axis] - rA[axis]) / direction[axis];
    
                    if (t1 > t2) {
                        std::swap(t1, t2);
                    }
    
                    t_min = std::max(t_min, t1);
                    t_max = std::min(t_max, t2);
    
                    if (t_min > t_max + tolerance) {
                        return false;
                    }
                }
            }
    
            const double t = 0.5 * (t_min + t_max);
            rPointInsideBox = rA + t * direction;
    
            return IsPointInsideBox(rPointInsideBox, box_min, box_max, tolerance);
        };
    
        auto is_point_inside_triangle = [&](
            const array_1d<double, 3>& rPoint,
            const array_1d<double, 3>& rA,
            const array_1d<double, 3>& rB,
            const array_1d<double, 3>& rC) -> bool
        {
            const array_1d<double, 3> v0 = rC - rA;
            const array_1d<double, 3> v1 = rB - rA;
            const array_1d<double, 3> v2 = rPoint - rA;
    
            const double dot00 = inner_prod(v0, v0);
            const double dot01 = inner_prod(v0, v1);
            const double dot02 = inner_prod(v0, v2);
            const double dot11 = inner_prod(v1, v1);
            const double dot12 = inner_prod(v1, v2);
    
            const double denominator = dot00 * dot11 - dot01 * dot01;
    
            if (std::abs(denominator) <= std::numeric_limits<double>::epsilon()) {
                return false;
            }
    
            const double inv_denominator = 1.0 / denominator;
            const double u = (dot11 * dot02 - dot01 * dot12) * inv_denominator;
            const double v = (dot00 * dot12 - dot01 * dot02) * inv_denominator;
    
            return u >= -1.0e-10 &&
                   v >= -1.0e-10 &&
                   u + v <= 1.0 + 1.0e-10;
        };
    
        auto triangle_box_intersection_point = [&](
            const array_1d<double, 3>& rA,
            const array_1d<double, 3>& rB,
            const array_1d<double, 3>& rC,
            array_1d<double, 3>& rPoint) -> bool
        {
            // 1. If a triangle vertex is inside the span box, use it.
            if (IsPointInsideBox(rA, box_min, box_max, tolerance)) {
                rPoint = rA;
                return true;
            }
    
            if (IsPointInsideBox(rB, box_min, box_max, tolerance)) {
                rPoint = rB;
                return true;
            }
    
            if (IsPointInsideBox(rC, box_min, box_max, tolerance)) {
                rPoint = rC;
                return true;
            }
    
            // 2. Check triangle edges against the span box.
            if (segment_box_intersection_point(rA, rB, rPoint)) {
                return true;
            }
    
            if (segment_box_intersection_point(rB, rC, rPoint)) {
                return true;
            }
    
            if (segment_box_intersection_point(rC, rA, rPoint)) {
                return true;
            }
    
            // 3. Check whether a box corner projects inside the triangle.
            array_1d<double, 3> edge_ab = rB - rA;
            array_1d<double, 3> edge_ac = rC - rA;

            const array_1d<double, 3> normal =
                MathUtils<double>::CrossProduct(edge_ab, edge_ac);

            const double normal_norm = norm_2(normal);
    
            if (normal_norm <= std::numeric_limits<double>::epsilon()) {
                return false;
            }
    
            const array_1d<double, 3> unit_normal = normal / normal_norm;
    
            std::array<array_1d<double, 3>, 8> box_corners = {{
                array_1d<double, 3>{box_min[0], box_min[1], box_min[2]},
                array_1d<double, 3>{box_max[0], box_min[1], box_min[2]},
                array_1d<double, 3>{box_min[0], box_max[1], box_min[2]},
                array_1d<double, 3>{box_max[0], box_max[1], box_min[2]},
                array_1d<double, 3>{box_min[0], box_min[1], box_max[2]},
                array_1d<double, 3>{box_max[0], box_min[1], box_max[2]},
                array_1d<double, 3>{box_min[0], box_max[1], box_max[2]},
                array_1d<double, 3>{box_max[0], box_max[1], box_max[2]}
            }};
    
            for (const auto& r_corner : box_corners) {
                const double signed_distance = inner_prod(unit_normal, r_corner - rA);
                const array_1d<double, 3> projected_point =
                    r_corner - signed_distance * unit_normal;
    
                if (IsPointInsideBox(projected_point, box_min, box_max, tolerance) &&
                    is_point_inside_triangle(projected_point, rA, rB, rC)) {
                    rPoint = projected_point;
                    return true;
                }
            }
    
            return false;
        };
    
        const array_1d<double, 3> p0 = rGeometry[0].Coordinates();
        const array_1d<double, 3> p1 = rGeometry[1].Coordinates();
        const array_1d<double, 3> p2 = rGeometry[2].Coordinates();
    
        if (triangle_box_intersection_point(p0, p1, p2, rPointInsideSpan)) {
            return true;
        }
    
        if (rGeometry.PointsNumber() == 4) {
            const array_1d<double, 3> p3 = rGeometry[3].Coordinates();
    
            if (triangle_box_intersection_point(p0, p2, p3, rPointInsideSpan)) {
                return true;
            }
        }
    
        return false;
    }
    //--------------------------------------------------------------


    void EnsureSkinNodesInExternalSpans(
        ModelPart& rSkinSubModelPart,
        const ExternalSpanDataMap& rExternalSpans,
        const KnotSpanGridInfo& rGridInfo) const;

    // type 1 element creation methods
    Vector CreateOpenUnitKnotVectorDegree1() const;

    Vector CreateOpenUnitKnotVector(
        const std::size_t PolynomialDegree) const;

    CanonicalFaceKey3D MakeCanonicalFaceKey3D(
        const IndexType NodeId0,
        const IndexType NodeId1,
        const IndexType NodeId2) const;

    CanonicalFaceKey3D MakeCanonicalFaceKey3D(
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2) const;

    IndexType GetNextGeometryId(
        const ModelPart& rRootModelPart) const;

    ModelPart& GetOrCreateSubModelPart(
        ModelPart& rRootModelPart,
        const std::string& rSubModelPartName) const;

    std::size_t ComputeNumberOfShapeFunctionsDerivatives(
        const ModelPart& rIgaModelPart) const;
    
    bool IsInnerSurrogateLoop(
        const ModelPart& rSurrogateSubModelPart) const;
    
    IndexType ComputeStartingBrepId(
        const ModelPart& rIgaModelPart,
        const ModelPart& rSurrogateSubModelPart) const;
    
    std::size_t ComputeBrepPatchLoopSize(
        const ModelPart& rSurrogateSubModelPart) const;
    
    BrepPatchData ComputeBrepPatchData(
        const ModelPart& rIgaModelPart,
        const IndexType BrepId) const;
    
    std::vector<BrepPatchData> BuildBrepPatchDataVector(
        const ModelPart& rIgaModelPart,
        const ModelPart& rSurrogateSubModelPart) const;
    
    const BrepPatchData& FindBrepPatchMatchingCondition(
        const Condition& rCondition,
        const std::vector<BrepPatchData>& rBrepPatchDataList) const;
    
    Geometry<Node>::Pointer CreateSurrogateFaceNeighbourGeometry(
        const BrepPatchData& rBrepPatchData,
        const std::size_t NumberOfShapeFunctionsDerivatives) const;

    NurbsSurfaceType::Pointer CreateType1CollapsedLateralCoonsSurface(
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pApexNode) const;

    NurbsVolumeType::Pointer CreateType1CollapsedPyramidVolume(
        const SurrogateFaceData& rFaceData,
        const NodePointerType& pApexNode) const;

    void SetType1GeometryData(
        Geometry<Node>& rGeometry,
        const IndexType ProjectionNodeId,
        const Geometry<Node>::Pointer& pNeighbourGeometry,
        const std::string& rIdentifier) const;

    void CreateAndTagSurfaceQuadraturePointGeometries(
        const NurbsSurfaceType::Pointer& pSurfaceGeometry,
        const Geometry<Node>::Pointer& pNeighbourGeometry) const;

    IntegrationPointsArrayType CreateCoonsVolumeGaussPoints(
        const std::size_t IntegrationOrder,
        const NurbsVolumeType& rGapVolume) const;

    IntegrationPointsArrayType CreateCoonsVolumeGaussPoints(
        const std::size_t IntegrationOrderU,
        const std::size_t IntegrationOrderV,
        const std::size_t IntegrationOrderW,
        const NurbsVolumeType& rGapVolume) const;
    
    GeometriesArrayType CreateAndTagVolumeQuadraturePointGeometries(
        const NurbsVolumeType::Pointer& pVolumeGeometry,
        const Geometry<Node>::Pointer& pNeighbourGeometry,
        const std::size_t IntegrationOrder,
        const std::size_t NumberOfShapeFunctionsDerivatives) const;

    GeometriesArrayType CreateAndTagAnisotropicVolumeQuadraturePointGeometries(
        const NurbsVolumeType::Pointer& pVolumeGeometry,
        const Geometry<Node>::Pointer& pNeighbourGeometry,
        const std::size_t IntegrationOrderU,
        const std::size_t IntegrationOrderV,
        const std::size_t IntegrationOrderW,
        const std::size_t NumberOfShapeFunctionsDerivatives) const;

    GeometriesArrayType CreateDirectMappedVolumeQuadraturePointGeometries(
        const Geometry<Node>::Pointer& pVolumeGeometry,
        const IntegrationPointsArrayType& rReferenceIntegrationPoints) const;

    GeometriesArrayType CreateNativeVolumeQuadraturePointGeometries(
        const Geometry<Node>::Pointer& pVolumeGeometry,
        const std::size_t IntegrationOrder) const;

    GeometriesArrayType CreateAndTagPyramidVolumeQuadraturePointGeometries(
        const std::array<NodePointerType, 4>& rBaseNodes,
        const NodePointerType& pApexNode,
        const std::size_t IntegrationOrder) const;

    GeometriesArrayType CreateAndTagLinearTetraVolumeQuadraturePointGeometries(
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const NodePointerType& pNode3,
        const std::size_t IntegrationOrder) const;

    GeometriesArrayType CreateAndTagCurvedType2QuadraturePointGeometries(
        const NurbsVolumeType::Pointer& pVolumeGeometry,
        const std::size_t IntegrationOrder) const;

    GeometriesArrayType CreateAndTagCurvedType3TriangleRadialQuadraturePointGeometries(
        const NodePointerType& pSurrogateNode,
        const GeometryType::Pointer& pTopSurfaceGeometry,
        const std::size_t IntegrationOrder) const;

    double CalculateType1CharacteristicLength(
        const SurrogateFaceData& rFaceData,
        const NodePointerType& pApexNode) const;

    void RegisterType1LateralFace(
        const SurrogateFaceData& rFaceData,
        const IndexType LocalEdgeIndex,
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pApexNode,
        const NurbsSurfaceType::Pointer& pLateralSurface,
        const Geometry<Node>::Pointer& pNeighbourGeometry);

    Type1CreationSummary ComputeType1CreationSummary(
        const std::size_t NumberOfPyramids) const;

    void PrintType1CreationSummary(
        const Type1CreationSummary& rSummary) const;

    void AddNeighbourGeometry(
        Geometry<Node>& rGeometry,
        const Geometry<Node>::Pointer& pNeighbourGeometry) const;

    void AddUniqueNeighbourGeometry(
        Geometry<Node>& rGeometry,
        const Geometry<Node>::Pointer& pCandidateNeighbourGeometry) const;


    IntegrationPointsArrayType CreateCoonsSurfaceGaussPoints(
        const std::size_t IntegrationOrder,
        const NurbsSurfaceType& rGapSurface) const;
    
    GeometriesArrayType CreateSurfaceQuadraturePointGeometries(
        const GeometryType::Pointer& pSurfaceGeometry,
        const std::size_t IntegrationOrder,
        const std::size_t NumberOfShapeFunctionsDerivatives) const;
    
    double CalculateSurfaceCharacteristicLength(
        const Geometry<Node>& rGeometry) const;
    
    bool ContainsNeighbourGeometry(
        const std::vector<Geometry<Node>::Pointer>& rNeighbourGeometries,
        const Geometry<Node>::Pointer& pCandidateNeighbourGeometry) const;
    
    void AddUniqueNeighbourGeometry(
        std::vector<Geometry<Node>::Pointer>& rNeighbourGeometries,
        const Geometry<Node>::Pointer& pCandidateNeighbourGeometry) const;
    
    std::vector<Geometry<Node>::Pointer> CollectUniqueNeighbourGeometries(
        const std::vector<LateralFaceOccurrence>& rOccurrences) const;
    
    //-------------------------------------------------------
    // type 2 elements
    NurbsVolumeType::Pointer CreateType2CollapsedEdgeVolume(
        const NodePointerType& pEdgeNode0,
        const NodePointerType& pEdgeNode1,
        const std::vector<NodePointerType>& rSkinEdgeControlNodes) const;

    NurbsVolumeType::Pointer CreateType2LinearCollapsedEdgeVolume(
        const NodePointerType& pEdgeNode0,
        const NodePointerType& pEdgeNode1,
        const NodePointerType& pSkinNode0,
        const NodePointerType& pSkinNode1) const;

    struct Type3CollapsedCornerData
    {
        NurbsVolumeType::Pointer pVolume;
        NurbsSurfaceType::Pointer pTopSurface;
        std::array<NodePointerType, 3> ProjectionNodes;
    };

    Type3CollapsedCornerData CreateType3CollapsedCornerVolume(
        const ModelPart& rSkinSubModelPart,
        const double MaximumProjectionDistance,
        const NodePointerType& pSurrogateNode,
        const NodePointerType& pProjectionNode0,
        const NodePointerType& pProjectionNode1,
        const NodePointerType& pProjectionNode2) const;

    Type3CollapsedCornerData CreateType3LinearCollapsedCornerVolume(
        const NodePointerType& pSurrogateNode,
        const NodePointerType& pProjectionNode0,
        const NodePointerType& pProjectionNode1,
        const NodePointerType& pProjectionNode2,
        const bool CreateTopSurface = true) const;

    Type3CollapsedCornerData CreateType3CollapsedCornerVolumeFromTopControlNodes(
        const NodePointerType& pSurrogateNode,
        const std::array<NodePointerType, 3>& rProjectionNodes,
        const std::vector<NodePointerType>& rTopControlNodes,
        const bool CreateTopSurface = true) const;

    void OrientTriangleFaceAwayFromOppositeNode(
        NodePointerType& rpNode0,
        NodePointerType& rpNode1,
        NodePointerType& rpNode2,
        const NodePointerType& pOppositeNode) const;

    void RegisterType3LateralFace(
        ModelPart& rDebugModelPart,
        IndexType& rNextGeometryId,
        const IndexType SurrogateConditionId,
        const SpanKey3D& rExternalSpan,
        NodePointerType pNode0,
        NodePointerType pNode1,
        NodePointerType pNode2,
        const NodePointerType& pOppositeNode,
        const Geometry<Node>::Pointer& pNeighbourGeometry,
        const bool HasType1NeighbourPath,
        const bool HasNeighbourActiveSpan,
        const SpanKey3D& rNeighbourActiveSpan);

    void CheckType3QuadraturePointGeometries(
        const GeometriesArrayType& rQuadraturePointGeometries,
        const NodePointerType& pSurrogateNode,
        const NodePointerType& pProjectionNode0,
        const NodePointerType& pProjectionNode1,
        const NodePointerType& pProjectionNode2) const;

    bool HaveCommonActiveSpan(
        const std::vector<SpanKey3D>& rFirst,
        const std::vector<SpanKey3D>& rSecond,
        const std::vector<SpanKey3D>& rThird) const;

    std::array<SpanKey3D, 4> GetSurroundingSpansOfGridEdge(
        const GridPointKey3D& rNode0,
        const GridPointKey3D& rNode1) const;

    bool IsValidSpan(
        const SpanKey3D& rSpan,
        const KnotSpanGridInfo& rGridInfo) const;

    bool AreFaceAdjacentSpans(
        const SpanKey3D& rA,
        const SpanKey3D& rB) const;

    double MaximumSkinProjectionDistance(
        const KnotSpanGridInfo& rGridInfo) const;

    bool ProjectPointToSkinBoundaryAlongDirection(
        const ModelPart& rSkinSubModelPart,
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rDirection,
        const double MaximumDistance,
        array_1d<double, 3>& rProjectionPoint) const;

    bool ProjectPointToClosestSkinBoundary(
        const ModelPart& rSkinSubModelPart,
        const array_1d<double, 3>& rPoint,
        array_1d<double, 3>& rProjectionPoint,
        double& rProjectionDistance) const;

    std::vector<NodePointerType> GetOrCreateQuadraticSkinEdgeControlNodes(
        ModelPart& rSkinSubModelPart,
        const NodePointerType& pSkinNode0,
        const NodePointerType& pSkinNode1,
        const double MaximumProjectionDistance,
        const std::string& rContext);

    NodePointerType CreateOrReuseAuxiliaryControlNode(
        ModelPart& rSkinSubModelPart,
        const array_1d<double, 3>& rPoint) const;

    std::vector<NodePointerType> CreateStraightSkinEdgeControlNodes(
        ModelPart& rSkinSubModelPart,
        const NodePointerType& pSkinNode0,
        const NodePointerType& pSkinNode1) const;

    std::vector<NodePointerType> GetOrCreateFinalSkinEdgeControlNodes(
        ModelPart& rSkinSubModelPart,
        const NodePointerType& pSkinNode0,
        const NodePointerType& pSkinNode1,
        const double MaximumProjectionDistance,
        const std::string& rContext);

    void UpdateLateralFaceRegistryGeometries(
        ModelPart& rRootModelPart,
        ModelPart& rSkinSubModelPart);

    NurbsSurfaceType::Pointer CreateCollapsedTriangleSurfaceWithCurvedEdge(
        const NodePointerType& pApexNode,
        const std::vector<NodePointerType>& rSkinEdgeControlNodes) const;

    NurbsSurfaceType::Pointer CreateCollapsedTriangleSurfaceWithCurvedTop(
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2) const;

    const std::vector<NodePointerType>* FindCachedSkinEdgeControlNodes(
        const NodePointerType& pNode0,
        const NodePointerType& pNode1) const;

    struct SkinProjectionTriangleData
    {
        array_1d<double, 3> A = ZeroVector(3);
        array_1d<double, 3> B = ZeroVector(3);
        array_1d<double, 3> C = ZeroVector(3);
        array_1d<double, 3> Center = ZeroVector(3);
        double Radius = 0.0;
    };

    int mEchoLevel = 0;
    std::size_t mGapApproximationOrder = 1;
    bool mStoreGapDebugGeometries = false;
    bool mUsePyramidQuadratureForType1 = true;
    bool mUseTetraQuadratureForLinearType2Type3 = true;
    bool mUseAnisotropicQuadratureForCurvedType2 = true;
    bool mUseTriangleRadialQuadratureForCurvedType3 = true;
    SkinEdgeControlPointMap mSkinEdgeControlNodes;
    std::set<SkinEdgeKey> mLinearSkinEdges;
    CurvedEdgeRegistry mCurvedEdgeRegistry;
    CurvedTopFaceRegistry mCurvedTopFaceRegistry;
    std::vector<SkinProjectionTriangleData> mSkinProjectionTriangleData;
    LateralFaceRegistry mLateralFaceRegistry;
};

} // namespace Kratos
