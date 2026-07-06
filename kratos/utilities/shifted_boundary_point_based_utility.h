//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/model.h"
#include "geometries/plane_3d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/key_hash.h"
#include <cstddef>

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace ShiftedBoundaryUtilityInternals {

    void VoigtTransformForProduct(
        const array_1d<double,3>& rVector,
        BoundedMatrix<double,2,3>& rVoigtMatrix);

    void VoigtTransformForProduct(
        const array_1d<double,3>& rVector,
        BoundedMatrix<double,3,6>& rVoigtMatrix);

    template <std::size_t TDim>
    void CalculateStrainMatrix(
        const Matrix& rDN_DX,
        const std::size_t& NumNodes,
        Matrix& rB);

    // ─────────────────────────────────────────────────────────────────────────────
    // Bounding box
    // ─────────────────────────────────────────────────────────────────────────────
    struct BBox {
        double min[3], max[3];

        BBox() {
            for (std::size_t i = 0; i < 3; ++i) {
                min[i] = 1e30; max[i] = -1e30;
            }
        }

        /**
         * @brief Calculate min and max coordinates from element geometry enlarging by a tolerance
         *
         * @param rGeom
         * @param Tolerance
         */
        BBox(const Element::GeometryType& rGeom, const double Tolerance=1e-12) {
            min[0] = max[0] = rGeom[0].X();
            min[1] = max[1] = rGeom[0].Y();
            min[2] = max[2] = rGeom[0].Z();
            for (std::size_t i = 1; i < rGeom.PointsNumber(); ++i) {
                min[0] = std::min(min[0], rGeom[i].X());
                min[1] = std::min(min[1], rGeom[i].Y());
                min[2] = std::min(min[2], rGeom[i].Z());
                max[0] = std::max(max[0], rGeom[i].X());
                max[1] = std::max(max[1], rGeom[i].Y());
                max[2] = std::max(max[2], rGeom[i].Z());
            }
            for (std::size_t d = 0; d < 3; ++d) {
                min[d] -= Tolerance; max[d] += Tolerance;
            }
        }

        void expand(double x, double y, double z) {
            min[0]=std::min(min[0],x); max[0]=std::max(max[0],x);
            min[1]=std::min(min[1],y); max[1]=std::max(max[1],y);
            min[2]=std::min(min[2],z); max[2]=std::max(max[2],z);
        }

        void expand(const BBox& o) {
            for (int i = 0; i < 3; ++i) {
                min[i]=std::min(min[i],o.min[i]);
                max[i]=std::max(max[i],o.max[i]);
            }
        }

        bool overlaps(const BBox& o) const {
            return max[0]>=o.min[0] && min[0]<=o.max[0] &&
                max[1]>=o.min[1] && min[1]<=o.max[1] &&
                max[2]>=o.min[2] && min[2]<=o.max[2];
        }

        double half_surface_area() const {
            double dx=max[0]-min[0], dy=max[1]-min[1], dz=max[2]-min[2];
            return dx*dy + dy*dz + dz*dx;
        }

        int longest_axis() const {
            double d[3]={max[0]-min[0], max[1]-min[1], max[2]-min[2]};
            return (d[0]>=d[1] && d[0]>=d[2]) ? 0 : (d[1]>=d[2] ? 1 : 2);
        }

        double centroid(int axis) const { return 0.5*(min[axis]+max[axis]); }
    };

    // ─────────────────────────────────────────────────────────────────────────────
    // Flat BVH node – 32 bytes, fits nicely in cache lines
    // ─────────────────────────────────────────────────────────────────────────────
    struct BVHNode {
        BBox  bbox;            // 24 bytes
        int   left_or_first;   //  4 bytes  (inner: left child index; leaf: first prim)
        int   right_or_count;  //  4 bytes  (inner: right child index; leaf: -count)

        bool is_leaf() const { return right_or_count <= 0; }
        int  prim_begin() const { return left_or_first; }
        int  prim_count() const { return -right_or_count; }
    };

    // ─────────────────────────────────────────────────────────────────────────────
    // Surface-mesh BVH  (built once, queried read-only)
    // ─────────────────────────────────────────────────────────────────────────────
    class ElementBVH {
    public:
        std::vector<BVHNode>   nodes;
        std::vector<int>       prim_ids;   // reordered primitive indices
        std::vector<BBox>      prim_boxes; // bbox per primitive (original order)

        static constexpr int MAX_LEAF_SIZE = 4;   // tune: 4-8 is usually optimal
        static constexpr int SAH_BINS      = 8;   // SAH binning resolution

        void build(const std::vector<BBox>& boxes);

        /**
         * @brief Query: append all primitive indices whose bbox overlaps query_box
         * Thread-safe (read-only after build).
         * @param query_box
         * @param results
         */
        void query(const BBox& query_box, std::vector<std::size_t>& results);

        /**
         * @brief update prim_boxes from new geometry, then refit the tree bottom-up
         * O(N), no reallocation; called every time geometry changes (e.g. moving mesh)
         * @param updated_boxes
         */
        void refit(const std::vector<BBox>& updated_boxes);

    private:
        /**
         * @brief Recursive build with SAH binning
         *
         * @param first
         * @param count
         * @return node_idx
         */
        int build_recursive(int first, int count);

        /**
         * @brief Recursive refit
         * Since tree is small after reuse iterative reverse traversal is cleaner and branchless
         * @param idx
         */
        void refit_node(int idx);
    };

}  // namespace ShiftedBoundaryUtilityInternals

///@}
///@name Kratos Classes
///@{

/**
 * @brief Utilities for the SBM boundary conditions via extension operators
 * This class provides the utilities for the calculation of the extension operator
 * in the Shifted Boundary Method Without Taylor Expansions allowing for a discontinuous interface (thin-walled structure),
 * which is defined by points on the surface of the interface (e.g. IGA integration points).
 */
class KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility
{
public:

    ///@name Type Definitions
    ///@{

    enum class ExtensionOperator
    {
        MLS,
        RBF
    };

    // variable types
    using NodeType = ModelPart::NodeType;
    using ElementType = Element;
    using GeometryType = ModelPart::GeometryType;

    // function types
    using MeshlessShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;
    using ElementSizeFunctionType = std::function<double(const GeometryType&)>;

    // set and vector types
    using NodesSetType = std::unordered_set<NodeType::Pointer, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;
    using ElementsSetType = std::unordered_set<ElementType::Pointer, SharedPointerHasher<ElementType::Pointer>, SharedPointerComparator<ElementType::Pointer>>;
    using CloudDataVectorType = std::vector<std::pair<NodeType::Pointer, double>>;
    using SkinPointsDataVectorType = std::vector<std::tuple< array_1d<double,3>, array_1d<double,3>, std::size_t >>; // vector of position, area normal and ID of skin points
    using EdgesVectorType = std::vector<std::pair<std::size_t, std::size_t>>;

    // map types
    using SkinPointsToElementsMapType = std::unordered_map<ElementType::Pointer, SkinPointsDataVectorType, SharedPointerHasher<ElementType::Pointer>, SharedPointerComparator<ElementType::Pointer>>;
    using SidesVectorToElementsMapType = std::unordered_map<ElementType::Pointer, Vector, SharedPointerHasher<ElementType::Pointer>, SharedPointerComparator<ElementType::Pointer>>;
    using AverageSkinToElementsMapType = std::unordered_map<ElementType::Pointer, std::pair<array_1d<double,3>, array_1d<double,3>>, SharedPointerHasher<ElementType::Pointer>, SharedPointerComparator<ElementType::Pointer>>;
    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ShiftedBoundaryPointBasedUtility
    KRATOS_CLASS_POINTER_DEFINITION(ShiftedBoundaryPointBasedUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Standard constructor
    /// @param rModel Model container
    /// @param ThisParameters Parameters object encapsulating the settings
    ShiftedBoundaryPointBasedUtility(
        Model& rModel,
        Parameters ThisParameters);

    /// Copy constructor.
    ShiftedBoundaryPointBasedUtility(ShiftedBoundaryPointBasedUtility const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ShiftedBoundaryPointBasedUtility& operator=(ShiftedBoundaryPointBasedUtility const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    //TODO
    void ResetFlags();

    /** TODO
     * @brief Set BOUNDARY flags for elements that are intersected by the tessellated skin geometry.
     * element BOUNDARY : elements in which a part of the surface is located (split/ intersected elements)
     */
    void FindElementsAtTessellatedBoundary(double Tolerance = 1e-10);

    //TODO
    void UpdateBoundaryElements(double Tolerance = 1e-10);

    //TODO
    void MapSkinPointsToElements();

    //TODO
    // Elements that were found to be touched or intersected by the tessellated skin geometry or containing skin integration points are flagged as BOUNDARY elements.
    void FlagBoundaryElements();

    // TODO Set the corresponding flags at the interface elements
    // element INTERFACE : elements owning the surrogate boundary nodes adjacent to an deactivated BOUNDARY element
    // Makes all elements surrounding a BOUNDARY element INTERFACE elements.
    // To be done after all necessary elements have been declared SBM_BOUNDARY (--> after FindElementsAtTessellatedBoundary and MapSkinPointsToElements)
    void FlagInterfaceElements();

    // TODO node ACTIVE : nodes that belong to the elements to be assembled (all nodes as the interface is discontinuous)
    // element ACTIVE : elements which are not BOUNDARY (the ones to be assembled)
    // Deactivates BOUNDARY elements and nodes that are not part of any ACTIVE element anymore.
    // If DeactivateUnstableClusters is true, than FindAndDeactivateUnstableClusters() is called.
    // To be done after all necessary elements have been declared SBM_BOUNDARY (--> after FindElementsAtTessellatedBoundary and MapSkinPointsToElements)
    void DeactivateElementsAndNodes(const bool DeactivateUnstableClusters);

    //TODO
    void CalculateAndAddSkinIntegrationPointConditions();

    //TODO
    // Calculate positive and negative side pressure and velocity and traction at the integration points of the skin geometry - for mpSkinPointsSubModelPart.
    // Results are stored in POSITIVE_FACE_PRESSURE, NEGATIVE_FACE_PRESSURE, POSITIVE_FACE_FLUID_VELOCITY, NEGATIVE_FACE_FLUID_VELOCITY, TRACTION_FROM_FLUID_PRESSURE, TRACTION_FROM_FLUID_STRESS, DRAG_FORCE.
    //TODO
    void CalculateVariablesAtSkinPoints();

    //TODO
    // Calls CalculateVariablesAtSkinPoints() and then calculates the values at the nodes of the skin model part using the skin points shape function values - for mpSkinDiscSubModelPart.
    // Results are stored in POSITIVE_FACE_PRESSURE, NEGATIVE_FACE_PRESSURE, POSITIVE_FACE_FLUID_VELOCITY, NEGATIVE_FACE_FLUID_VELOCITY.
    void CalculateVariablesAtSkinPointsAndNodes();

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    const Parameters GetDefaultParameters() const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "ShiftedBoundaryPointBasedUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ShiftedBoundaryPointBasedUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;
    ModelPart* mpBoundarySubModelPart = nullptr;

    //NOTE that right now integration points of mpSkinDiskretSubModelPart are taken as skin points with area and normal for SkinPointsSubModelPart.
    ModelPart* mpSkinDiscSubModelPart = nullptr;
    ModelPart* mpSkinPointsSubModelPart = nullptr;
    std::string mSkinModelPartName = "";

    ShiftedBoundaryUtilityInternals::ElementBVH mSkinBVH;
    std::vector<ElementType::Pointer> mBVHIdxToSkinElementVector;

    ElementsSetType mBoundaryElementsSet;
    SkinPointsToElementsMapType mSkinPointsMap;

    SidesVectorToElementsMapType mSidesVectorMap;
    NodesCloudMapType mExtensionOperatorMap;

    ExtensionOperator mExtensionOperator;
    std::size_t mMLSExtensionOperatorOrder;

    const Condition* mpConditionPrototype;

    bool mPositiveSideIsEnclosed = false;
    bool mNegativeSideIsEnclosed = false;

    /// @brief Protected empty constructor for derived classes
    //ShiftedBoundaryPointBasedUtility() {}

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief TODO. This method requires the skin to be stored in mpSkinDiscSubModelPart as a Kratos model part with elements, integration points and area normals.
     *
     * @tparam TDim Working space dimension
     * @param rSkinPointsMap
     */
    template <std::size_t TDim>
    void MapSkinPointsToElementsTemplated(SkinPointsToElementsMapType& rSkinPointsMap);

    //TODO
    bool LocatePoint(
        ShiftedBoundaryUtilityInternals::ElementBVH& rCandidatesBvh,
        std::vector<ElementType::Pointer>& rIdxToPointer,
        array_1d<double,3>& rPointCoords,
        ElementType::Pointer& pElement);

    //TODO
    void FindAndDeactivateUnstableClusters();

    //TODO - move to separate process?
    std::vector<std::vector<ElementType::Pointer>> FindClusters();

    //TODO - move to separate process?
    EdgesVectorType GetActiveAdjacencyGraph();

    /**TODO*/
    void SetSidesVectorsAndSkinNormalsForSplitElements(
        const SkinPointsToElementsMapType& rSkinPointsMap,
        SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap);

    /** TODO
     * @brief Set the Extension Operators For Split Element Nodes object selecting a support cloud for each node of each split element and using MLS shape functions
     *
     * @param rExtensionOperatorMap
     */
    void SetExtensionOperatorsForSplitElementNodes(
        const SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap,
        NodesCloudMapType& rExtensionOperatorMap);

    //TODO
    /**
     * @brief Set the support cloud for the same side of the boundary of which the first nodes nodes are given.
     * For given nodes on one side of the boundary, this function creates a cloud of nodes around them on the same side of the boundary.
     * These first nodes could be the positive nodes of a BOUNDARY element in order to get a support cloud and extension basis for one of the element's negative nodes and vice versa.
     * The support cloud created by this function is to be used for calculating the MLS-based extension operator.
     * @param rSameSideNodes Pointer to the first nodes of the support cloud
     * TODO
     * @param rCloudNodes Vector containing pointers to the nodes of the cloud
     * @param rCloudCoordinates Matrix containing the coordinates of the nodes of the cloud
     */
    void SetLateralSupportCloud(
        const NodeType::Pointer pOtherSideNode,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates,
        const Kratos::Flags& rSearchSideFlag);

    //TODO without dot product
    void AddLateralSupportLayer(
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet);

    //TODO with dot product
    void AddLateralSupportLayer(
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet);

    /**
     * @brief Create a pointer vector of pointers to all the nodes affecting the respective side of a split element's boundary.
     * Create vectors to pointers of cloud nodes for a split element using all the nodes affecting the element.
     * Pointer Vectors are sorted by ID to properly get the extension operator data
     *
     * @param rElement
     * @param rSidesVector
     * @param rExtensionOperatorMap
     * @param rCloudNodeVectorPositiveSide
     * @param rCloudNodeVectorNegativeSide
     */
    void CreateCloudNodeVectorsForSplitElement(
        const ElementType& rElement,
        const Vector& rSidesVector,
        NodesCloudMapType& rExtensionOperatorMap,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide);

    /* TODO */
    void GetDataForSplitElementSkinPoint(
        const ElementType& rElement,
        const array_1d<double,3>& rIntPtCoordinates,
        Vector& rIntPtShapeFunctionValues,
        Matrix& rIntPtShapeFunctionDerivatives);

    /* TODO */
    bool AddIntegrationPointCondition(
        const ElementType& rElement,
        const Vector& rSidesVector,
        const double ElementSize,
        const array_1d<double,3>& rIntPtCoordinates,
        const array_1d<double,3>& rIntPtAreaNormal,
        NodesCloudMapType& rExtensionOperatorMap,
        const PointerVector<NodeType>& rCloudNodeVector,
        const Vector& rIntPtShapeFunctionValues,
        const Matrix& rIntPtShapeFunctionDerivatives,
        std::size_t& r_ConditionId,
        const bool ConsiderPositiveSide);

    // TODO
    bool SetEnclosedNodesPressure(
        ElementType& rElement,
        const Vector& rSidesVector,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal);

    //TODO
    // Calculate positive and negative side pressure at the nodes of the skin model part
    // Result is stored in POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE.
    //TODO
    template <std::size_t TDim>
    void CalculateVariablesAtSkinPointsTemplated();

    //TODO
    template <std::size_t TDim>
    void CalculateVariablesAtSkinPointsAndNodesTemplated();

    //TODO
    template <std::size_t TDim>
    bool CalculateUnknownsForBothSidesOfSplitElement(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns);

    /**
     * @brief Get the MLS shape functions factory object
     * This function returns a prototype for the MLS shape functions calculation
     * @return MLSShapeFunctionsFunctionType MLS shape functions call prototype
     */
    MeshlessShapeFunctionsFunctionType GetMLSShapeFunctionsFunction() const;

    /**
     * @brief Get the RBF shape functions factory object
     * This function returns a prototype for the RBF shape functions calculation
     * @return RBFShapeFunctionsFunctionType RBF shape functions call prototype
     */
    MeshlessShapeFunctionsFunctionType GetRBFShapeFunctionsFunction() const;

    /**
     * @brief Get the element size function object
     * This function returns a prototype for the element size calculation call
     * @param rGeometry Geometry to calculate the element size
     * @return ElementSizeFunctionType Element size calculation call
     */
    ElementSizeFunctionType GetElementSizeFunction(const GeometryType& rGeometry);

    /**
     * @brief Calculates the kernel function radius
     * For the given cloud of points coordinates, this function calculates the maximum distance between
     * the center of the kernel (origin) and the points. This is supposed to be used in the MLS
     * approximation calculation.
     * @param rCloudCoordinates Matrix containing the coordinates of the nodes of the cloud
     * @param rOrigin Coordinates of the point on which the kernel function is centered
     * @return double Kernel function radius
     */
    double CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin);

    /**
     * @brief Get the number of required points
     * For the MLS approximation, this function returns the minimum number of required points according to
     * dimension and order. If the cloud of points has less points than the value returned by this function
     * the case is ill-conditioned, meaning that the cloud needs to be enlarged.
     * @return std::size_t Number of required points
     */
    std::size_t GetRequiredNumberOfPoints();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class ShiftedBoundaryPointBasedUtility

}  // namespace Kratos.
