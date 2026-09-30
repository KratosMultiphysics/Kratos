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
#include "includes/model_part.h"
#include <cstddef>
#if defined(_MSC_VER) && !defined(_BACKUP_ITERATOR_DEBUG_LEVEL)
#define _BACKUP_ITERATOR_DEBUG_LEVEL 0
#endif
#include <span>
#include <string>
#include <vector>

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
    // Mesh BVH  (built once, queried read-only) - BVH: Bounding Volume Hierarchy
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

    // ─────────────────────────────────────────────────────────────────────────────
    // Node CSR  (built once, queried read-only) - CSR: Compressed Sparse Row
    // ─────────────────────────────────────────────────────────────────────────────
    struct NodeCSR {
        // maps global node Id → local index (0..n_local_nodes-1)
        std::unordered_map<std::size_t,std::size_t> global_to_local;
        // maps local index → node pointer
        std::vector<ModelPart::NodeType::Pointer> local_nodes;

        // CSR arrays for neighboring nodes and for active neighboring nodes only
        std::vector<std::size_t> neighbors; // flat neighbor list (local indices)
        std::vector<std::size_t> offset;    // offset[i]..offset[i+1] = neighbors of local node i
        std::vector<std::size_t> active_neighbors; // flat active neighbor list (local indices)
        std::vector<std::size_t> active_offset;    // offset[i]..offset[i+1] = active neighbors of local node i

        void build(const std::vector<ModelPart::ElementType::Pointer>& rElements)
        {
            // Temporary adjacency lists, a set is used to avoid duplicates
            std::vector<std::unordered_set<std::size_t>> adj;
            std::vector<std::unordered_set<std::size_t>> active_adj;

            // Collect unique nodes and add their adjacency, omitting deactivated elements, thus separating both sides
            for (auto p_elem : rElements) {
                const auto& r_geom  = p_elem->GetGeometry();
                const std::size_t n_pts = r_geom.PointsNumber();
                std::array<std::size_t, 4> local_indices;  // max 4 for tet

                for (std::size_t i = 0; i < n_pts; ++i) {
                    auto p_node = r_geom(i);
                    auto [it, inserted] = global_to_local.emplace(p_node->Id(), local_nodes.size());
                    if (inserted) {
                        local_nodes.push_back(p_node);
                        adj.emplace_back();   // add empty neighbor list for this node
                        active_adj.emplace_back();   // add empty neighbor list for this node
                    }
                    local_indices[i] = it->second;
                }

                // Add edges between all pairs of nodes in this element
                if (p_elem->Is(ACTIVE)) {
                    for (std::size_t i = 0; i < n_pts; ++i) {
                        for (std::size_t j = i+1; j < n_pts; ++j) {
                            adj[local_indices[i]].insert(local_indices[j]);
                            adj[local_indices[j]].insert(local_indices[i]);
                            active_adj[local_indices[i]].insert(local_indices[j]);
                            active_adj[local_indices[j]].insert(local_indices[i]);
                        }
                    }
                } else {
                    for (std::size_t i = 0; i < n_pts; ++i) {
                        for (std::size_t j = i+1; j < n_pts; ++j) {
                            adj[local_indices[i]].insert(local_indices[j]);
                            adj[local_indices[j]].insert(local_indices[i]);
                        }
                    }
                }
            }

            // Pack into CSR
            const std::size_t n_nodes = local_nodes.size();
            offset.resize(n_nodes + 1, 0);
            active_offset.resize(n_nodes + 1, 0);
            for (std::size_t i = 0; i < n_nodes; ++i) {
                offset[i + 1] = offset[i] + adj[i].size();
                active_offset[i + 1] = active_offset[i] + active_adj[i].size();
            }
            neighbors.resize(offset[n_nodes]);
            active_neighbors.resize(active_offset[n_nodes]);
            for (std::size_t i = 0; i < n_nodes; ++i) {
                std::copy(adj[i].begin(), adj[i].end(), neighbors.begin() + offset[i]);
                std::copy(active_adj[i].begin(), active_adj[i].end(), active_neighbors.begin() + active_offset[i]);
            }
        }

        bool get_neighbors(const std::size_t node_id, std::vector<ModelPart::NodeType::Pointer>& neighboring_nodes) const {
            // Get local index of the node if it is in the CSR
            auto it = global_to_local.find(node_id);
            if (it == global_to_local.end()) { return false; }
            const std::size_t local_idx = it-> second;

            // Get all neighboring node pointers of the given node
            for (std::size_t neigh_idx : std::span<const std::size_t>(neighbors.data() + offset[local_idx], offset[local_idx+1] - offset[local_idx])) {
                neighboring_nodes.push_back(local_nodes[neigh_idx]);
            }
            return true;
        }

        bool get_active_neighbors(const std::size_t node_id, std::vector<ModelPart::NodeType::Pointer>& neighboring_nodes) const {
            // Get local index of the node if it is in the CSR
            auto it = global_to_local.find(node_id);
            if (it == global_to_local.end()) { return false; }
            const std::size_t local_idx = it-> second;

            // Get all neighboring node pointers of the given node
            for (std::size_t neigh_idx : std::span<const std::size_t>(active_neighbors.data() + active_offset[local_idx], active_offset[local_idx+1] - active_offset[local_idx])) {
                neighboring_nodes.push_back(local_nodes[neigh_idx]);
            }
            return true;
        }
    };

}  // namespace ShiftedBoundaryUtilityInternals

///@}
///@name Kratos Classes
///@{

/**
 * @brief Utilities for the SBM boundary conditions via sampling-point extension operators
 * This class provides the utilities for the calculation of the extension operator
 * in the Shifted Boundary Method (SBM) without Taylor expansions allowing for a double-sided boundary (thin-walled structure).
 * Moving-least-squares (MLS) or radial basis function (RBF) can be choosen as extension operators to extrapolate the solution values.
 * The boundary is defined by points of the boundary surface (e.g. IGA integration points).
 * A discretization of the surface is still needed as well (1D surface for 2D background mesh and 2D surface for 3D background mesh).
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
    using NodesWeightsMapType = std::unordered_map<NodeType::Pointer, Vector, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;
    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ShiftedBoundaryPointBasedUtility
    KRATOS_CLASS_POINTER_DEFINITION(ShiftedBoundaryPointBasedUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Shifted Boundary Point Based Utility object (standard constructor)
     * Gets the model part and the parameters from the input and initializes the SBM utility.
     * The boundary sub model part and skin model parts are accessed by name and stored in mpBoundarySubModelPart, mpSkinDiscSubModelPart and mpSkinPointsSubModelPart.
     * mpBoundarySubModelPart should be empty. Boundary conditions as stated in "boundary_wall_condition_name" and stored in mpConditionPrototype
     * will be added to this sub model part when calling CalculateAndAddSkinIntegrationPointConditions().
     * mpSkinDiscSubModelPart should already contain the discretized skin geometry.
     * mpSkinPointsSubModelPart should contain the skin integration points before MapSkinPointsToElements() is called.
     * The order and type of the extension operator is stored in mMLSExtensionOperatorOrder and mExtensionOperator (default: linear MLS).
     * @param rModel Model container
     * @param ThisParameters Parameters object encapsulating the settings
     */
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

    /**
     * @brief Resets all SBM flags (to false) and reactivates all elements and nodes.
     * Initializes all elements and nodes to ACTIVE state and clears the SBM_BOUNDARY and SBM_INTERFACE flags.
     */
    void ResetFlags();

    /**
     * @brief Identifies elements intersected by the tessellated skin geometry as boundary elements and stores them in mBoundaryElementsSet.
     * Constructs a Bounding Volume Hierarchy (BVH) of the skin discretization (mSkinBVH) from mpSkinDiscSubModelPart and performs
     * a spatial search to identify background mesh elements that overlap with or contain the skin geometry. An element is treated
     * as intersected if the geometry's HasIntersection method returns true, which it should if any of the element's nodes or faces
     * are intersected by or coincide with the discretized skin geometry.
     * @param Tolerance Tolerance for coordinate comparisons in bounding boxes (default: 1e-10)
     */
    void FindElementsAtTessellatedBoundary(double Tolerance = 1e-10);

    /**
     * @brief Identifies boundary elements after mesh movement using BVH refitting and stores them in mBoundaryElementsSet.
     * Instead of building a new BVH, this method refits the existing tree mSkinBVH with updated
     * skin geometry coordinates and their bounding boxes. Then, it performs the same spatial search as
     * FindElementsAtTessellatedBoundary() to identify intersected elements.
     * This method is more efficient when the boundary moves but the topology remains unchanged.
     * @param Tolerance Tolerance for coordinate comparisons in bounding boxes (default: 1e-10)
     */
    void UpdateBoundaryElements(double Tolerance = 1e-10);

    /**
     * @brief Maps skin integration points to the background mesh elements containing them.
     * This method is a wrapper for the 2D and 3D templated implementations of MapSkinPointsToElementsTemplated().
     * It locates each skin point (commonly the integration points of the skin discretization) within
     * the background mesh elements. Elements containing skin integration points are stored in
     * mSkinPointsMap and mapped to their skin points data.
     * Must be called after FindElementsAtTessellatedBoundary().
     */
    void MapSkinPointsToElements();

    /**
     * @brief Marks all elements in mBoundaryElementsSet as SBM_BOUNDARY.
     * Flags all elements found to be intersected by the tessellated skin geometry or containing
     * skin integration points with the SBM_BOUNDARY flag. These are the elements in which the true
     * boundary is located (split/intersected elements) and will be deactivated from the main assembly.
     */
    void FlagBoundaryElements();

    /**
     * @brief Identifies elements as SBM_INTERFACE which are adjacent to SBM_BOUNDARY elements and therefore the surrogate boundary.
     * Flags all background elements that share a face with a SBM_BOUNDARY element as SBM_INTERFACE.
     * These interface elements contain nodes of the surrogate boundary.
     * Must be called after FindElementsAtTessellatedBoundary() and MapSkinPointsToElements() and can
     * be called once for all SBM utilities as it does not use the boundary elements set of the current utility.
     */
    void FlagInterfaceElements();

    /**
     * @brief Deactivates SBM_BOUNDARY elements and isolated nodes to not assemble their contributions to the system.
     * Marks elements flagged as SBM_BOUNDARY as NOT ACTIVE and deactivates nodes that no longer
     * belong to any ACTIVE element.
     * Optionally deactivates unstable element clusters calling FindAndDeactivateUnstableClusters().
     * Must be called after FlagBoundaryElements() and can be called once for all SBM utilities.
     * @param DeactivateUnstableClusters If true, additionally finds and deactivates unstable element clusters.
     */
    void DeactivateElementsAndNodes(const bool DeactivateUnstableClusters);

    /**
     * @brief Creates and adds SBM wall conditions at skin integration points.
     * Iterates over all boundary elements containing skin points and:
     * 1. Determines the positive and negative sides of the boundary for each node (SetSidesForSkinPointElements())
     * 2. Sets up extension operators (MLS or RBF) for extrapolating solution values (SetExtensionForSkinPointElements())
     * 3. Fixes the pressure of the first node in an enclosed region (FixPressureOfEnclosedNode())
     * 4. Calculates shape functions and derivatives at each skin point using extension operators and
     *    adds a wall conditions to mpBoundarySubModelPart (CalculateExtendedValuesAndAddSkinPointConditions())
     */
    void CalculateAndAddSkinIntegrationPointConditions();

    /**
     * @brief Frees the pressure constraint from a node (mEnclosedNodeId) if mEnclosedPressureIsSet is true.
     * If a domain is enclosed (mPositiveSideIsEnclosed or mNegativeSideIsEnclosed is true) and
     * a pressure constraint was applied in CalculateAndAddSkinIntegrationPointConditions()
     * to a node of ID mEnclosedNodeId in that region, then this function removes the pressure constraint.
     * This is necessary in case of a moving boundary.
     */
    void FreePressureOfEnclosedNode();

    /**
     * @brief Calculates pressure, velocity and traction on both sides of the boundary at skin integration points.
     * This method is a wrapper for the 2D and 3D templated implementations of CalculateVariablesAtSkinPointsTemplated().
     * The results are assigned to nodes in mpSkinPointsSubModelPart.
     */
    void CalculateVariablesAtSkinPoints();

    /**
     * @brief Calculates values at skin integration points and pressure and velocity at nodes of the discretized skin.
     * This method is a wrapper for the 2D and 3D templated implementations of CalculateVariablesAtSkinPointsAndNodesTemplated().
     * Values are assigned to skin points in mpSkinPointsSubModelPart and the nodes in mpSkinDiscSubModelPart.
     */
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

    // Background mesh
    ModelPart* mpModelPart = nullptr;
    // Model part containing the boundary conditions at the skin integration points
    ModelPart* mpBoundarySubModelPart = nullptr;

    std::string mSkinModelPartName = "";
    // Model part containing the discretized skin geometry
    ModelPart* mpSkinDiscSubModelPart = nullptr;
    // Model part containing the skin integration points as nodes
    ModelPart* mpSkinPointsSubModelPart = nullptr;

    // Bounding Volume Hierarchy (BVH) for skin discretization
    ShiftedBoundaryUtilityInternals::ElementBVH mSkinBVH;
    std::vector<ElementType::Pointer> mBVHIdxToSkinElementVector;

    // Set of boundary elements intersected by the skin geometry (SBM_BOUNDARY elements)
    ElementsSetType mBoundaryElementsSet;
    // Vector of boundary elements and several layers of elements surrounding them
    std::vector<ElementType::Pointer> mSurroundingElementsVector;

    // Map background mesh elements containing skin integration points and their skin point data (position, area normal, ID)
    SkinPointsToElementsMapType mSkinPointsMap;

    // Map background mesh elements containing skin integration points to a vector of signed values indicating which side of the boundary each node belongs to
    SidesVectorToElementsMapType mSidesVectorMap;

    // Map nodes of elements containing skin integration points to their extension (support cloud nodes and their weights)
    NodesCloudMapType mExtensionOperatorMap;            // support nodes from the other side in order to extrapolate the solution values of the other side
    NodesCloudMapType mExtensionOperatorSameSideMap;    // support nodes from their own side for a deactivated node

    // Extension operator type (MLS or RBF) and order for MLS
    ExtensionOperator mExtensionOperator;
    std::size_t mMLSExtensionOperatorOrder;

    // Boundary wall condition prototype used to create new conditions for skin integration points
    const Condition* mpConditionPrototype;

    // Bools indicating whether the positive or negative side of the boundary is enclosed and whether a pressure constraint has been applied to a node in an enclosed region
    bool mPositiveSideIsEnclosed = false;
    bool mNegativeSideIsEnclosed = false;
    bool mEnclosedPressureIsSet = false;
    std::size_t mEnclosedNodeId = 0;

    /// @brief Protected empty constructor for derived classes
    //ShiftedBoundaryPointBasedUtility() {}

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Maps skin integration points to the background mesh elements containing them.
     * The boundary elements in mBoundaryElementsSet and two layers of elements around them (equals 1 complete layer of nodes) are used as
     * background element candidates for locating the skin integration points. Those elements are stored in mSurroundingElementsVector for
     * CalculateAndAddSkinIntegrationPointConditions.
     * A BVH is built for the background element candidates and used to efficiently locate each skin point in mpSkinPointsSubModelPart (LocatePoint()).
     * This method requires the skin integration points to be stored in mpSkinPointsSubModelPart as nodes with coordinates and area normals.
     * After sorting, elements containing skin integration points are stored in mSkinPointsMap and mapped to their skin points data (positions, area normals, IDs).
     * Must be called after FindElementsAtTessellatedBoundary().
     * @tparam TDim Working space dimension (2 or 3)
     * @param rSkinPointsMap Output map associating background elements to their contained skin integration points
     */
    template <std::size_t TDim>
    void MapSkinPointsToElementsTemplated(SkinPointsToElementsMapType& rSkinPointsMap);

    /**
     * @brief Locates a point in the mesh using BVH spatial indexing.
     * Performs point location within the background mesh using a prebuilt Bounding Volume Hierarchy (BVH).
     * Finds the element containing the given point coordinates by first querying the BVH for candidate
     * elements, then performing IsInside() checks.
     * @param rCandidatesBvh BVH structure for candidate element selection
     * @param rIdxToPointer Vector mapping BVH indices to element pointers
     * @param rPointCoords Coordinates of the point to locate
     * @param pElement Output pointer to the element containing the point (if found)
     * @return bool True if point was successfully located in an element, false otherwise
     */
    bool LocatePoint(
        ShiftedBoundaryUtilityInternals::ElementBVH& rCandidatesBvh,
        std::vector<ElementType::Pointer>& rIdxToPointer,
        array_1d<double,3>& rPointCoords,
        ElementType::Pointer& pElement);

    /**
     * @brief Identifies and deactivates unstable clusters of background elements.
     * Finds connected components of background elements (FindClusters()) and deactivates all clusters
     * except for the largest one. This prevents unstable small element groups from
     * remaining in the simulation.
     * TODO: deactivate clusters that do not have a fixed DOF instead?
     */
    void FindAndDeactivateUnstableClusters();

    /**
     * @brief Identifies connected clusters of ACTIVE background mesh elements.
     * Groups ACTIVE elements into connected clusters based on adjacency (GetActiveAdjacencyGraph())
     * using Union-Find. Elements are considered connected if they share an edge.
     * Returns a vector of clusters, where each cluster is a vector of element pointers.
     * TODO: move to separate process?
     * @return std::vector<std::vector<ElementType::Pointer>> Vector of element clusters
     */
    std::vector<std::vector<ElementType::Pointer>> FindClusters();

    /**
     * @brief Constructs an adjacency graph of ACTIVE elements of mpModelPart.
     * Creates a vector of edges representing connectivity between ACTIVE elements.
     * Two elements are connected if they share an edge (GetValue(NEIGHBOUR_ELEMENTS)).
     * Used for cluster identification and determining element connectivity.
     * @return EdgesVectorType Vector of element Id() pairs representing edges in the adjacency graph
     */
    EdgesVectorType GetActiveAdjacencyGraph();

    /**
     * @brief Determines and stores which side of the boundary each node of an element with skin points belongs to in mSidesVectorMap.
     * For each boundary element containing skin points from mSkinPointsMap, determines the positive/negative side
     * classification for each of its nodes. Creates a vector (one entry per node) where positive
     * values indicate positive side and negative values indicate negative side. The vector is stored in mSidesVectorMap.
     * For the classification, all skin integration points around a node contribute a weighted side vote based on their area and distance to the node.
     * The side for each node is then determined by the sign of the sum of weighted votes and afterwards corrected using the sides of active neighboring nodes.
     * Then, if a boundary node has no skin points around it, its side is determined by the sides of its active neighbors.
     * The flag SBM_BOUNDARY is set for positive side nodes and the flag SBM_INTERFACE is set for negative side nodes, which is required for SetLateralSupportCloud().
     * The method also computes the average skin point position and normal for each boundary element.
     * @param rSurroundingNodesCSR CSR providing node connectivity of the boundary nodes and at least one layer around it
     * @param rAvgSkinMap Output map storing average position and normal for each element with skin points
     */
    void SetSidesForSkinPointElements(
        const ShiftedBoundaryUtilityInternals::NodeCSR& rSurroundingNodesCSR,
        AverageSkinToElementsMapType& rAvgSkinMap);

    /**
     * @brief Fixes the pressure DOF of a node in an enclosed volume.
     * For enclosed volumes, fixes the pressure of an element's node to zero if it is on the enclosed side.
     * This removes one degree of freedom to ensure the problem is well-posed (pressure is defined up to a constant).
     * @param rElement Boundary element being checked
     * @param rSidesVector Vector indicating positive/negative side classification for each node
     * @return bool True if pressure was successfully fixed, false otherwise
     */
    bool FixPressureOfEnclosedNode(
        ElementType& rElement,
        const Vector& rSidesVector);

    /**
     * @brief Creates extension operators for nodes of elements containing skin points using MLS or RBF approximation.
     * For each node of an element with skin points in mSidesVectorMap, builds a support cloud from neighboring nodes on the opposite
     * side of the boundary calling SetLateralSupportCloud. These clouds are used to compute MLS or RBF shape functions
     * for support node weights in order to extrapolate solution values to the node.
     * For inactive nodes, the support cloud is built from neighboring nodes on the same side of the boundary.
     * The support nodes and weights for a node are stored in mExtensionOperatorMap if the extension is for
     * the other side of the boundary and in mExtensionOperatorSameSideMap if the extension is for the same side of the boundary.
     * @param rSurroundingNodesCSR CSR providing node connectivity of the boundary nodes and at least one layer around it
     * @param rAvgSkinMap Map containing average boundary position and normal for each element containing skin points
     */
    void SetExtensionForSkinPointElements(
        const ShiftedBoundaryUtilityInternals::NodeCSR& rSurroundingNodesCSR,
        AverageSkinToElementsMapType& rAvgSkinMap);

    /**
     * @brief Constructs a support cloud of nodes for pBaseNode on the side of the boundary given by rSearchSideFlag.
     * For a given node (pBaseNode), this function collects a cloud of nodes neighboring it on the specified side of the boundary.
     * The base node could be a positive side node of a BOUNDARY element in order to get a support cloud and extension basis
     * for the negative side collecting support nodes on the negative side and vice versa.
     * It could also be a positive inactive node of a BOUNDARY element in order to get a support cloud for the positive side
     * collecting support nodes on the positive side.
     * The support cloud created by this function is to be used for calculating an extension operator.
     * @param pBaseNode The node for which to build the support cloud
     * @param rAvgSkinPosition Average position of the boundary in the element
     * @param rAvgSkinNormal Average normal vector of the boundary in the element
     * @param rSurroundingNodesCSR CSR providing node connectivity of the boundary nodes and at least one layer around it
     * @param rCloudNodes Output vector of node pointers in the support cloud
     * @param rCloudCoordinates Output matrix of coordinates for cloud nodes
     * @param rSearchSideFlag Flag indicating which side of boundary to search (positive/negative)
     */
    void SetLateralSupportCloud(
        const NodeType::Pointer pBaseNode,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const ShiftedBoundaryUtilityInternals::NodeCSR& rSurroundingNodesCSR,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates,
        const Kratos::Flags& rSearchSideFlag);

    /**
     * @brief Adds neighboring nodes to the support cloud (expanding outward without directional filtering).
     * Expands the support cloud by adding all neighbors of the current layer of nodes without crossing the boundary (as SBM_BOUNDARY elements).
     * Here, neighboring nodes are all nodes of neighboring elements (GetValue(NEIGHBOUR_ELEMENTS)).
     * @param PreviousLayerNodes Nodes added in the previous expansion step
     * @param CurrentLayerNodes Output vector to store newly added nodes for this layer
     * @param SupportNodesSet Set tracking all nodes already added to the cloud
     */
    void AddLateralSupportLayer(
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet);

    /**
     * @brief Adds neighboring nodes to the support cloud (expanding in normal direction of the boundary using a surrounding node CSR).
     * Expands the support cloud by adding active neighbors of current layer nodes from the surrounding node CSR
     * in normal direction of the averaged boundary by computing a dot product.
     * @param rAvgSkinPosition Average position of the boundary
     * @param rAvgSkinNormal Average normal vector of the boundary (used to filter nodes in normal direction of the boundary)
     * @param rSurroundingNodesCSR CSR providing node connectivity of the boundary nodes and at least one layer around it
     * @param PreviousLayerNodes Nodes added in the previous expansion step
     * @param CurrentLayerNodes Output vector to store newly added nodes for this layer
     * @param SupportNodesSet Set tracking all nodes already added to the cloud
     */
    void AddLateralSupportLayer(
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const ShiftedBoundaryUtilityInternals::NodeCSR& rSurroundingNodesCSR,
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet);

    /**
     * @brief Calculates extended shape function values and derivatives and adds SBM boundary conditions for all skin integration points.
     * Requires mSkinPointsMap and mSidesVectorMap. For each element containing skin points extension operators are collected
     * (CreateCloudVectorsForSkinPointElement()).
     * Then, for each skin integration point:
     * 1. The velocity is obtained from the skin point node in mpSkinPointsSubModelPart as EMBEDDED_VELOCITY.
     * 2. The background element's shape function values for the skin point (GetDataForSkinPointInElement()) are multiplied by the extension weights.
     * 3. For valid extensions and the sides of the boundary that are not enclosed, a boundary condition is added (AddSkinPointCondition()).
     */
    void CalculateExtendedValuesAndAddSkinPointConditions();

    /**
     * @brief Collects support nodes and weights for both sides of a boundary element separately.
     * Extracts MLS/RBF support nodes and weights for each node of the element on both positive
     * and negative sides. Determine whether extension of the element is valid for each side independently based on whether the necessary
     * extensions were found in mExtensionOperatorMap. Pointers to support nodes are stored for both sides separately.
     * Weights are stored for both sides separately in a map for each support node as a Vector of weights (one weight per node of the element).
     * The output vectors and maps are used to calculate extended shape function values and derivatives at skin integration points.
     * @param rElement The boundary element being processed
     * @param rSidesVector Vector indicating positive/negative side classification for each node
     * @param rCloudNodeVectorPositiveSide Output vector of support nodes for positive side
     * @param rCloudNodeVectorNegativeSide Output vector of support nodes for negative side
     * @param rCloudWeightsMapPositiveSide Output map of MLS/RBF weights for positive side nodes
     * @param rCloudWeightsMapNegativeSide Output map of MLS/RBF weights for negative side nodes
     * @param rPositiveExtensionValid Output flag indicating if positive side extension is valid
     * @param rNegativeExtensionValid Output flag indicating if negative side extension is valid
     */
    void CreateCloudVectorsForSkinPointElement(
        const ElementType& rElement,
        const Vector& rSidesVector,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide,
        NodesWeightsMapType& rCloudWeightsMapPositiveSide,
        NodesWeightsMapType& rCloudWeightsMapNegativeSide,
        bool& rPositiveExtensionValid,
        bool& rNegativeExtensionValid);

    /**
     * @brief Computes shape function values and derivatives at a skin point within an element.
     * Calculates the standard Lagrange shape function values and derivatives of a point
     * within the reference coordinate system of an element. Used to determine how the
     * element's nodal solution contributes to values at the skin point location.
     * @param rElement The element containing the skin point
     * @param rSkinPtCoordinates Physical coordinates of the skin point
     * @param rSkinPtShapeFunctionValues Output vector of shape function values
     * @param rSkinPtShapeFunctionDerivatives Output matrix of shape function derivatives
     */
    void GetDataForSkinPointInElement(
        const ElementType& rElement,
        const array_1d<double,3>& rSkinPtCoordinates,
        Vector& rSkinPtShapeFunctionValues,
        Matrix& rSkinPtShapeFunctionDerivatives);

    /**
     * @brief Creates and adds a boundary condition for a skin integration point.
     * Constructs a new condition at the specified skin point location using the extended shape function values
     * and derivatives for one side of the boundary. The condition is added to mpBoundarySubModelPart.
     * Normal, integration weight (area) and point velocity are assigned.
     * The condition ID is incremented for each new condition added.
     * @param rElement Background element containing the skin point
     * @param ElementSize Characteristic size of the backgroundelement
     * @param rSkinPtCoordinates Physical coordinates of the skin point
     * @param rSkinPtAreaNormal Area normal vector at the skin point
     * @param rSkinPtVelocity Velocity vector at the skin point
     * @param rCloudNodeVector Support nodes used for extension
     * @param rSkinPtShapeFunctionValuesExtended Extended shape function values
     * @param rSkinPtShapeFunctionDerivativesExtended Extended shape function derivatives
     * @param r_ConditionId Input/output counter for condition ID assignment
     * @return bool True if condition was successfully added, false if insufficient integration weight or zero normal
     */
    bool AddSkinPointCondition(
        const ElementType& rElement,
        const double ElementSize,
        const array_1d<double,3>& rSkinPtCoordinates,
        const array_1d<double,3>& rSkinPtAreaNormal,
        const array_1d<double,3>& rSkinPtVelocity,
        const PointerVector<NodeType>& rCloudNodeVector,
        const Vector& rSkinPtShapeFunctionValuesExtended,
        const Matrix& rSkinPtShapeFunctionDerivativesExtended,
        std::size_t& r_ConditionId);

    /**
     * @brief Calculates pressure, velocity and traction on both sides of the boundary at skin integration points.
     * It gets pressure and velocity values at each skin point (in mSkinPointsMap) for the positive and negative side
     * of the boundary (CalculateUnknownsForSkinPointElement()).
     * It also calculates the traction and drag forces contributions from both sides.
     * Results are stored in the following variables:
     * - POSITIVE_FACE_PRESSURE, NEGATIVE_FACE_PRESSURE
     * - POSITIVE_FACE_FLUID_VELOCITY, NEGATIVE_FACE_FLUID_VELOCITY
     * - TRACTION_FROM_FLUID_PRESSURE, TRACTION_FROM_FLUID_STRESS, DRAG_FORCE
     * The results are assigned to nodes in mpSkinPointsSubModelPart.
     * @tparam TDim Working space dimension (2 or 3)
     */
    template <std::size_t TDim>
    void CalculateVariablesAtSkinPointsTemplated();

    /**
     * @brief Calculates values at skin integration points and pressure and velocity at nodes of the discretized skin.
     * First calls CalculateVariablesAtSkinPoints() to evaluate variables at skin integration points,
     * then interpolates the positive and negative side pressure and velocity to the nodes of the skin discretization mesh
     * with shape function values at the skin integration points from querying mSkinBVH and using methods of the skin geometry.
     * Results at the skin nodes are stored in:
     * - POSITIVE_FACE_PRESSURE, NEGATIVE_FACE_PRESSURE
     * - POSITIVE_FACE_FLUID_VELOCITY, NEGATIVE_FACE_FLUID_VELOCITY
     * Values are assigned to skin points in mpSkinPointsSubModelPart and the nodes in mpSkinDiscSubModelPart.
     * @tparam TDim Working space dimension (2 or 3)
     */
    template <std::size_t TDim>
    void CalculateVariablesAtSkinPointsAndNodesTemplated();

    /**
     * @brief Calculates positive and negative side unknowns (velocity and pressure) of a boundary element (containing skin points).
     * Evaluate the unknowns for each node by getting the node's side from the side vector of the element (from mSidesVectorMap) and then
     * calculating the unknowns at the node from mExtensionOperatorSameSideMap if it's not ACTIVE
     * and by using mExtentionOperatorMap for the other side of the boundary.
     * @tparam TDim Working space dimension (2 or 3)
     * @param pElement Pointer to the boundary element
     * @param rPositiveSideUnknowns Output vector of pressure and velocity values on positive side
     * @param rNegativeSideUnknowns Output vector of pressure and velocity values on negative side
     * @return bool True if calculation was successful (sides vector and extensions were found)
     */
    template <std::size_t TDim>
    bool CalculateUnknownsForSkinPointElement(
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
