//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//                   Nicolo' Antonelli
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
#include <cmath>
#include <limits>
#include <memory>

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

    class Gap3DData;

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

    static std::size_t FlattenSpanColumnIndex(
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const std::size_t SpanIndexY,
        const std::size_t SpanIndexZ = 0)
    {
        return SpanIndexY + rSkinBinsPerSpan.NumberOfSpansY * SpanIndexZ;
    }

    static std::size_t FindSpanNnzIndex(
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const std::size_t SpanIndexX,
        const std::size_t SpanIndexY,
        const std::size_t SpanIndexZ = 0)
    {
        return FindNnzIndex(
            rSkinBinsPerSpan.Occupancy,
            SpanIndexX,
            FlattenSpanColumnIndex(rSkinBinsPerSpan, SpanIndexY, SpanIndexZ));
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
        SnakeSbmProcess::CreateTheSnakeCoordinates(true);
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
    bool SuccessfullProjection = false;
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
    std::size_t mNumberOfInterpolationLevels = 3; // number of additional levels of interpolation for the gap elements (in addition to the corner nodes)
    std::size_t mGapApproximationOrder = 1; // Order of quadrature of the Gap elemennts and conditions
    std::size_t mGapInterpolationOrder; // Order of the Gap elements
    std::string mGapSbmType; // "default" (uses ordered points of the curve) or
                             // "interpolation" (interpolates the points on the curve -- no assumption on their order)
                             // "sbm" (uses SBM on top of the "interpolation" approach to further optimize the gap geometry)
    double mMinSpanSize = 0.0;
    bool mUseForMultipatch = false;
    bool mUseForLocalRefinement = false;

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
    std::vector<IndexType> AssestProjectionsFeasibility(
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

    /**
     * @brief Returns Gauss points mapped on the gap Coons volume.
     * @param Order Quadrature order used in each parametric direction.
     * @param rGapVolume NURBS volume representing the current gap Coons volume.
     * @return Array of integration points with Jacobian-weighted volume weights.
     */
    IntegrationPointsArrayType CreateCoonsVolumeGaussPoints(
        std::size_t Order,
        const NurbsVolumeType& rGapVolume) const;

    /**
     * @brief Prints diagnostic information about the orientation of a 3D gap Coons volume.
     * @return True if the determinant scan identifies a locally folded volume.
     * @param rGapVolume Gap volume to inspect.
     * @param rNode00 Surrogate face node at (u,v,w) = (0,0,0).
     * @param rNode01 Surrogate face node at (u,v,w) = (0,1,0).
     * @param rNode10 Surrogate face node at (u,v,w) = (1,0,0).
     * @param rNode11 Surrogate face node at (u,v,w) = (1,1,0).
     * @param rProjectedPoint00 Projected skin point at (u,v,w) = (0,0,1).
     * @param rProjectedPoint01 Projected skin point at (u,v,w) = (0,1,1).
     * @param rProjectedPoint10 Projected skin point at (u,v,w) = (1,0,1).
     * @param rProjectedPoint11 Projected skin point at (u,v,w) = (1,1,1).
     * @param pProjectionSurface00_01 Lateral projection surface over edge rNode00-rNode01.
     * @param pProjectionSurface00_10 Lateral projection surface over edge rNode00-rNode10.
     * @param pProjectionSurface01_11 Lateral projection surface over edge rNode01-rNode11.
     * @param pProjectionSurface10_11 Lateral projection surface over edge rNode10-rNode11.
     * @param rIntegrationPoints Integration points used for determinant scanning.
     */
    bool CheckGapVolumeOrientation(
        const NurbsVolumeType& rGapVolume,
        const NodeType& rNode00,
        const NodeType& rNode01,
        const NodeType& rNode10,
        const NodeType& rNode11,
        const array_1d<double,3>& rProjectedPoint00,
        const array_1d<double,3>& rProjectedPoint01,
        const array_1d<double,3>& rProjectedPoint10,
        const array_1d<double,3>& rProjectedPoint11,
        const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
        const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
        const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
        const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
        const IntegrationPointsArrayType& rIntegrationPoints) const;

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
template <bool TIsInnerLoop>
std::vector<array_1d<double,3>> CollectSkinUVBetween(
    const ModelPart& rSkinSubModelPart,
    const NurbsSurfaceType& rSurface,
    IndexType id_node_1,
    IndexType id_node_2,
    const std::string& rLayerName = "",
    const std::string& rConditionName = "") const
{
    std::vector<array_1d<double,3>> uv_pts;
    uv_pts.reserve(rSkinSubModelPart.NumberOfConditions());

    const auto ordered_node_ids = CollectOrderedSkinNodePath( //FIXME: go back to previous version
        rSkinSubModelPart,
        id_node_1,
        id_node_2,
        TIsInnerLoop,
        rLayerName,
        rConditionName);

    for (const auto id : ordered_node_ids) {
        const auto& p_node = rSkinSubModelPart.pGetNode(id);
        array_1d<double,3> uv = p_node->Coordinates();
        uv_pts.push_back(uv);
    }
    return uv_pts;
}



/**
 * @brief Computes the chord-length parametrization on [0,1] for UV samples.
 * @param UV Collection of UV points describing the sampled curve.
 * @return Normalized cumulative arc-length parameters.
 */
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
void BernsteinBasis(const int p, const double t, std::vector<double>& B) const;

/**
 * @brief Builds a Bézier knot vector with multiplicity p+1 at the boundaries.
 * @param p Polynomial degree of the Bézier curve.
 * @param rKnots Vector receiving the generated knot values.
 */
void BuildBezierKnots(const int p, Vector& rKnots) const;

/**
 * @brief Solves a symmetric positive definite linear system using a simple Cholesky factorization.
 * @param A Reference to the system matrix, overwritten with its factorization.
 * @param b Right-hand side vector, overwritten with the solution.
 * @return True when the system is successfully solved.
 */
bool SolveSPD(std::vector<std::vector<double>>& A, std::vector<double>& b) const;

/**
 * @brief Creates a UV Bézier curve of degree p from a set of control points.
 * @param CtrlsUV Control points expressed in UV coordinates.
 * @param p Polynomial degree of the resulting curve.
 * @return Pointer to the constructed Bézier curve geometry.
 */
typename NurbsCurveGeometryType::Pointer MakeBezierUV_FromControls(
    const std::vector<array_1d<double,3>>& CtrlsUV, const int p) const;

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
    const double ridge = 1e-12) const;

/**
 * @brief Reverses the orientation of a UV Bézier curve.
 * @param p_forward Pointer to the curve whose orientation is reversed.
 * @return Pointer to a curve geometry with inverted control point order.
 */
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
template <bool TIsInnerLoop>
typename NurbsCurveGeometryType::Pointer FitUV_BetweenSkinNodes_Generic(
    const ModelPart& rSkinSubModelPart,
    const NurbsSurfaceType& rSurface,
    IndexType id_node_1,
    IndexType id_node_2,
    const int p,
    const double ridge = 1e-12,
    const std::string& rLayerName = "",
    const std::string& rConditionName = "") const
{
    auto UV = CollectSkinUVBetween<TIsInnerLoop>(
        rSkinSubModelPart,
        rSurface,
        id_node_1,
        id_node_2,
        rLayerName,
        rConditionName);
    KRATOS_ERROR_IF(UV.size() < (std::size_t)(p+1))
        << "FitUV_BetweenSkinConditions_Generic: not enough samples for degree p=" << p << ".\n";

    // Ensure endpoints match exactly the first and last
    UV.front() = UV.front();
    UV.back()  = UV.back();

    return FitBezierUV_LS_Generic(UV, p, ridge);
}

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

void CreateInnerSkinLocalRefinementCouplingConditions(
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    NurbsSurfaceType::Pointer& pNurbsSurface);


}; // Class SnakeGapSbmProcess

class KRATOS_API(IGA_APPLICATION) SnakeGapSbmProcess::Gap3DData
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

    using SkinEdgeNormalMap = std::unordered_map<
        SkinEdgeKey,
        array_1d<double, 3>,
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

    enum class VolumeQuadratureRule
    {
        Pyramid,
        LinearTetrahedron,
        CurvedType2,
        CurvedType3TriangleRadial,
        NurbsVolume
    };

    struct VolumeGeometryData
    {
        NurbsVolumeType::Pointer pVolumeGeometry;
        GeometryType::Pointer pTopSurfaceGeometry;
        std::vector<NodePointerType> RuleNodes;
        std::vector<Geometry<Node>::Pointer> NeighbourGeometries;

        double CharacteristicLength = 0.0;
        std::size_t IntegrationOrder = 1;
        std::size_t NumberOfShapeFunctionsDerivatives = 1;
        int GapType = 0;
        VolumeQuadratureRule QuadratureRule = VolumeQuadratureRule::NurbsVolume;
    };

    struct PhysicalQuadraturePoint
    {
        array_1d<double, 3> Coordinates = ZeroVector(3);
        double Weight = 0.0;
    };

    struct MomentFittingSettings
    {
        std::size_t PolynomialDegree = 1;
        int TotalDegree = -1;
        bool UseTensorProductBasis = true;
        int TensorProductDegree = 3;
        std::size_t CandidateIntegrationOrder = 3;
        std::size_t MaximumCandidateIntegrationOrder = 5;
        double RankTolerance = 1.0e-12;
        double ResidualTolerance = 1.0e-8;
        double MaximumAbsoluteWeightRatio = 100.0;
        std::size_t MaximumSelectedPointMultiplier = 2;
    };

    struct MomentFittedVolumeQuadratureResult
    {
        GeometriesArrayType QuadraturePointGeometries;
        std::size_t NumberOfVolumeGeometries = 0;
        std::size_t NumberOfCandidatePoints = 0;
        std::size_t NumberOfSelectedPoints = 0;
        std::size_t NumberOfMoments = 0;
        std::size_t EstimatedRank = 0;
        std::size_t NumberOfReferencePoints = 0;
        double RelativeMomentResidual = 0.0;
        double FittedVolume = 0.0;
        double ReferenceVolume = 0.0;
        double AbsoluteWeightRatio = 0.0;
        double LastAttemptAbsoluteWeightRatio = 0.0;
        double MinimumWeight = 0.0;
        double MaximumWeight = 0.0;
        array_1d<double, 4> ReferenceLinearMoments = ZeroVector(4);
        array_1d<double, 4> FittedLinearMoments = ZeroVector(4);
        bool UsedFallback = false;
    };

    struct Type1VolumeQuadratureData : public VolumeGeometryData
    {
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
        Geometry<Node>::Pointer pSurfaceGeometry;
        std::vector<Geometry<Node>::Pointer> NeighbourGeometries;

        double CharacteristicLength = 0.0;

        int GapType = 0;
        bool IsType3TopFace = false;
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
    struct Type2VolumeQuadratureData : public VolumeGeometryData
    {
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
    struct Type3VolumeQuadratureData : public VolumeGeometryData
    {
        NodePointerType pSurrogateNode;
        NodePointerType pProjectionNode0;
        NodePointerType pProjectionNode1;
        NodePointerType pProjectionNode2;

        IndexType SurrogateNodeId = 0;
        IndexType SurrogateConditionId = 0;
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

    explicit Gap3DData(const int EchoLevel = 0);

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

    void SetInitialNurbsSkinModelPart(
        const ModelPart& rInitialSkinModelPart);

    void SetMaximumNurbsProjectionDistance(
        const double MaximumProjectionDistance);

    bool TryProjectPointToInitialNurbsSkin(
        const array_1d<double, 3>& rPoint,
        array_1d<double, 3>& rProjectionPoint,
        double& rProjectionDistance,
        array_1d<double, 3>* pProjectionNormal = nullptr,
        const array_1d<double, 3>* pInitialGuess = nullptr) const;

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

    GeometriesArrayType CreateAdaptiveInterfaceSurfaceQuadrature(
        const GeometryType::Pointer& pSurfaceGeometry,
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
        const std::size_t NumberOfShapeFunctionsDerivatives);

    MomentFittedVolumeQuadratureResult CreateMomentFittedVolumeQuadrature(
        const std::vector<const VolumeGeometryData*>& rVolumeGeometries,
        const MomentFittingSettings& rSettings) const;

    GeometriesArrayType CreateOriginalVolumeQuadrature(
        const std::vector<const VolumeGeometryData*>& rVolumeGeometries) const;

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

    /**
     * @brief Build the surrogate face-centre quadrature geometries used as
     * NEIGHBOUR_GEOMETRIES by the type1/type2/type3 workflow.
     *
     * The returned geometry is deliberately not the linear surrogate
     * condition geometry: it carries the NURBS control-point connectivity and
     * shape functions evaluated at the BREP face centre.
     */
    std::unordered_map<IndexType, Geometry<Node>::Pointer>
    CreateSurrogateFaceCentreNeighbourGeometries(
        const ModelPart& rSurrogateSubModelPart) const;

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
            << "[SnakeGapSbmProcess::Gap3DData::SpanBox] Span "
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

    GeometriesArrayType CreateOriginalVolumeQuadrature(
        const VolumeGeometryData& rVolumeGeometry) const;

    std::vector<PhysicalQuadraturePoint>
    CreateOriginalPhysicalVolumeQuadrature(
        const VolumeGeometryData& rVolumeGeometry) const;

    std::vector<PhysicalQuadraturePoint>
    CreateMappedPhysicalVolumeQuadrature(
        const Geometry<Node>::Pointer& pVolumeGeometry,
        const IntegrationPointsArrayType& rReferenceIntegrationPoints) const;

    std::vector<PhysicalQuadraturePoint>
    CreateNativePhysicalVolumeQuadrature(
        const Geometry<Node>::Pointer& pVolumeGeometry,
        const std::size_t IntegrationOrder) const;

    std::vector<PhysicalQuadraturePoint>
    CreateCurvedType3TriangleRadialPhysicalQuadrature(
        const NodePointerType& pSurrogateNode,
        const GeometryType::Pointer& pTopSurfaceGeometry,
        const std::size_t IntegrationOrder) const;

    GeometriesArrayType MaterializePhysicalVolumeQuadrature(
        const std::vector<PhysicalQuadraturePoint>& rPhysicalPoints) const;

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

    bool ProjectPointToClosestNurbsSurface(
        const array_1d<double, 3>& rPoint,
        array_1d<double, 3>& rProjectionPoint,
        double& rProjectionDistance,
        const bool IncludeKnotSpanBoundaryLines = false,
        const double MaximumProjectionDistance =
            std::numeric_limits<double>::max(),
        array_1d<double, 3>* pProjectionNormal = nullptr,
        const array_1d<double, 3>* pInitialGuess = nullptr) const;

    bool ProjectPointToInitialNurbsSkinAlongDirection(
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rDirection,
        const double MaximumDistance,
        array_1d<double, 3>& rProjectionPoint) const;

    void RebuildSkinEdgeAverageNormalMap(
        const ModelPart& rSkinSubModelPart);

    NodePointerType FindClosestSkinMeshNode(
        const ModelPart& rSkinSubModelPart,
        const array_1d<double, 3>& rPoint,
        double& rDistance) const;

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

    bool HasGeometricEdgeCurvature(
        const std::vector<NodePointerType>& rLinearControlNodes,
        const std::vector<NodePointerType>& rCurrentControlNodes) const;

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

    struct NurbsSurfaceProjectionData
    {
        array_1d<double, 3> ControlPointMin = ZeroVector(3);
        array_1d<double, 3> ControlPointMax = ZeroVector(3);
    };

    int mEchoLevel = 0;
    std::size_t mGapApproximationOrder = 1;
    bool mStoreGapDebugGeometries = false;
    bool mUsePyramidQuadratureForType1 = true;
    bool mUseTetraQuadratureForLinearType2Type3 = true;
    bool mUseAnisotropicQuadratureForCurvedType2 = true;
    bool mUseTriangleRadialQuadratureForCurvedType3 = true;
    CurvedEdgeRegistry mCurvedEdgeRegistry;
    CurvedTopFaceRegistry mCurvedTopFaceRegistry;
    std::vector<SkinProjectionTriangleData> mSkinProjectionTriangleData;
    std::vector<NurbsSurfaceType::Pointer> mInitialNurbsSkinSurfaces;
    std::vector<NurbsSurfaceProjectionData> mInitialNurbsSkinProjectionData;
    double mMaximumNurbsProjectionDistance =
        std::numeric_limits<double>::max();
    SkinEdgeNormalMap mSkinEdgeAverageNormals;
    LateralFaceRegistry mLateralFaceRegistry;
};


}  // namespace Kratos.
