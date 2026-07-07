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

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/pointer_vector.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/exception.h"
#include "includes/kratos_flags.h"
#include "includes/lock_object.h"
#include "includes/node.h"
#include "includes/smart_pointers.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
#include "input_output/logger.h"
#include "processes/tetrahedra_mesh_edge_swapping_process.h"
#include "utilities/element_size_calculator.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/shifted_boundary_point_based_utility.h"
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "geometries/plane_3d.h"
#include <cstddef>
#include <numeric>
#include <omp.h>
#include <unordered_map>
#include <utility>
#include <vector>


namespace Kratos
{
    namespace ShiftedBoundaryUtilityInternals {

        // see FluidElementUtilities 2D
        void VoigtTransformForProduct(
            const array_1d<double,3>& rVector,
            BoundedMatrix<double,2,3>& rVoigtMatrix)
        {
            rVoigtMatrix.clear();

            rVoigtMatrix(0,0) = rVector(0);
            rVoigtMatrix(0,2) = rVector(1);
            rVoigtMatrix(1,1) = rVector(1);
            rVoigtMatrix(1,2) = rVector(0);
        }

        // see FluidElementUtilities 3D
        void VoigtTransformForProduct(
            const array_1d<double,3>& rVector,
            BoundedMatrix<double,3,6>& rVoigtMatrix)
        {
            rVoigtMatrix.clear();

            rVoigtMatrix(0,0) = rVector(0);
            rVoigtMatrix(0,3) = rVector(1);
            rVoigtMatrix(0,5) = rVector(2);
            rVoigtMatrix(1,1) = rVector(1);
            rVoigtMatrix(1,3) = rVector(0);
            rVoigtMatrix(1,4) = rVector(2);
            rVoigtMatrix(2,2) = rVector(2);
            rVoigtMatrix(2,4) = rVector(1);
            rVoigtMatrix(2,5) = rVector(0);
        }

        // see FluidElementUtilities 2D
        template <>
        void CalculateStrainMatrix<2>(
            const Matrix& rDN_DX,
            const std::size_t& NumNodes,
            Matrix& rB)
        {
            rB.clear();
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node ) {
                const std::size_t col = 3 * i_node;
                rB(0, col  ) = rDN_DX(i_node, 0);
                rB(1, col+1) = rDN_DX(i_node, 1);
                rB(2, col  ) = rDN_DX(i_node, 1);
                rB(2, col+1) = rDN_DX(i_node, 0);
            }
        }

        // see FluidElementUtilities 3D
        template <>
        void CalculateStrainMatrix<3>(
            const Matrix& rDN_DX,
            const std::size_t& NumNodes,
            Matrix& rB)
        {
            rB.clear();
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node ) {
                const std::size_t col = 4 * i_node;
                rB(0, col  ) = rDN_DX(i_node, 0);
                rB(1, col+1) = rDN_DX(i_node, 1);
                rB(2, col+2) = rDN_DX(i_node, 2);
                rB(3, col  ) = rDN_DX(i_node, 1);
                rB(3, col+1) = rDN_DX(i_node, 0);
                rB(4, col+1) = rDN_DX(i_node, 2);
                rB(4, col+2) = rDN_DX(i_node, 1);
                rB(5, col  ) = rDN_DX(i_node, 2);
                rB(5, col+2) = rDN_DX(i_node, 0);
            }
        }

        void ElementBVH::build(const std::vector<BBox>& boxes) {
            const int N = static_cast<int>(boxes.size());
            prim_boxes = boxes;
            prim_ids.resize(N);
            std::iota(prim_ids.begin(), prim_ids.end(), 0);

            nodes.clear();
            nodes.reserve(2 * N);           // upper bound on node count
            build_recursive(0, N);
        }

        void ElementBVH::query(const BBox& query_box, std::vector<std::size_t>& results) {
            if (nodes.empty()) return;
            // Iterative traversal with explicit stack (avoids recursion overhead)
            int stack[64];
            int top = 0;
            stack[top++] = 0;

            while (top > 0) {
                const BVHNode& n = nodes[stack[--top]];
                if (!n.bbox.overlaps(query_box)) continue;

                if (n.is_leaf()) {
                    const std::size_t end = n.prim_begin() + n.prim_count();
                    for (std::size_t i = n.prim_begin(); i < end; ++i)
                        results.push_back(prim_ids[i]);
                } else {
                    // Push both children; visit the one with smaller bbox first
                    // (minor heuristic – avoids checking obviously far child first)
                    stack[top++] = n.left_or_first;
                    stack[top++] = n.right_or_count;
                }
            }
        }

        void ElementBVH::refit(const std::vector<BBox>& updated_boxes) {
            // Step 1: update leaf-level primitive boxes
            prim_boxes = updated_boxes;

            // Step 2: refit all nodes bottom-up
            // The flat array is built in pre-order (parent before children),
            // so iterating in REVERSE order processes children before parents
            refit_node(0);
        }

        int ElementBVH::build_recursive(int first, int count) {
            int node_idx = static_cast<int>(nodes.size());
            nodes.emplace_back();
            BVHNode& node = nodes[node_idx];

            // Compute bounding box of all primitives in [first, first+count)
            node.bbox = BBox{};
            for (int i = first; i < first + count; ++i)
                node.bbox.expand(prim_boxes[prim_ids[i]]);

            // Leaf condition
            if (count <= MAX_LEAF_SIZE) {
                node.left_or_first  = first;
                node.right_or_count = -count;   // negative → leaf
                return node_idx;
            }

            // SAH split
            int   best_axis  = -1;
            int   best_split = -1;
            float best_cost  = std::numeric_limits<float>::max();

            for (int axis = 0; axis < 3; ++axis) {
                // Bin centroids
                float ax_min = node.bbox.min[axis];
                float ax_max = node.bbox.max[axis];
                if (ax_max - ax_min < 1e-12f) continue;

                struct Bin { BBox bbox; int count = 0; };
                std::array<Bin, SAH_BINS> bins;

                float inv_range = SAH_BINS / (ax_max - ax_min);
                for (int i = first; i < first + count; ++i) {
                    float c   = prim_boxes[prim_ids[i]].centroid(axis);
                    int   bin = std::min(static_cast<int>((c - ax_min) * inv_range),
                                        SAH_BINS - 1);
                    bins[bin].bbox.expand(prim_boxes[prim_ids[i]]);
                    bins[bin].count++;
                }

                // Sweep left→right, then right→left to evaluate split costs
                std::array<float, SAH_BINS-1> left_area,  right_area;
                std::array<int,   SAH_BINS-1> left_count, right_count;

                BBox  running_box; int running_cnt = 0;
                for (int b = 0; b < SAH_BINS-1; ++b) {
                    running_box.expand(bins[b].bbox);
                    running_cnt += bins[b].count;
                    left_area [b] = running_box.half_surface_area();
                    left_count[b] = running_cnt;
                }
                running_box = BBox{}; running_cnt = 0;
                for (int b = SAH_BINS-2; b >= 0; --b) {
                    running_box.expand(bins[b+1].bbox);
                    running_cnt += bins[b+1].count;
                    right_area [b] = running_box.half_surface_area();
                    right_count[b] = running_cnt;
                }

                float parent_area = node.bbox.half_surface_area();
                for (int b = 0; b < SAH_BINS-1; ++b) {
                    if (left_count[b] == 0 || right_count[b] == 0) continue;
                    float cost = (left_area[b]  * left_count[b] +
                                right_area[b] * right_count[b]) / parent_area;
                    if (cost < best_cost) {
                        best_cost  = cost;
                        best_axis  = axis;
                        best_split = b;
                    }
                }
            }

            // Fallback: split at median along longest axis
            if (best_axis == -1) {
                best_axis  = node.bbox.longest_axis();
                best_split = SAH_BINS / 2 - 1;
            }

            // Partition prim_ids in-place
            float ax_min   = node.bbox.min[best_axis];
            float ax_max   = node.bbox.max[best_axis];
            float inv_range= SAH_BINS / (ax_max - ax_min + 1e-30f);

            auto mid_it = std::partition(
                prim_ids.begin() + first,
                prim_ids.begin() + first + count,
                [&](int id) {
                    float c   = prim_boxes[id].centroid(best_axis);
                    int   bin = std::min(static_cast<int>((c - ax_min) * inv_range),
                                        SAH_BINS - 1);
                    return bin <= best_split;
                });

            int left_count_val = static_cast<int>(mid_it - (prim_ids.begin() + first));

            // Guard against degenerate splits (all in one side)
            if (left_count_val == 0 || left_count_val == count) {
                node.left_or_first  = first;
                node.right_or_count = -count;
                return node_idx;
            }

            // Build children (NB: node reference may be invalidated by push_back,
            // so we save what we need and re-fetch by index afterwards)
            int left_child  = build_recursive(first, left_count_val);
            int right_child = build_recursive(first + left_count_val,
                                            count  - left_count_val);

            // Re-fetch node (vector may have reallocated during recursion)
            nodes[node_idx].left_or_first  = left_child;
            nodes[node_idx].right_or_count = right_child;

            return node_idx;
        }

        void ElementBVH::refit_node(int idx) {
            BVHNode& node = nodes[idx];

            if (node.is_leaf()) {
                // Recompute bbox from primitives
                node.bbox = BBox{};
                const int end = node.prim_begin() + node.prim_count();
                for (int i = node.prim_begin(); i < end; ++i)
                    node.bbox.expand(prim_boxes[prim_ids[i]]);
            } else {
                refit_node(node.left_or_first);
                refit_node(node.right_or_count);
                // Parent bbox = union of children
                node.bbox = nodes[node.left_or_first].bbox;
                node.bbox.expand(nodes[node.right_or_count].bbox);
            }
        }
    }  // namespace ShiftedBoundaryUtilityInternals

    //----------------------------------------------------------------
    //     PUBLIC METHODS
    //----------------------------------------------------------------

    const Parameters ShiftedBoundaryPointBasedUtility::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "skin_model_part_name" : "",
            "boundary_sub_model_part_name" : "",
            "boundary_wall_condition_name" : "",
            "extension_operator_type" : "MLS",
            "mls_extension_operator_order" : 1,
            "enclosed_area" : "none",
            "use_tessellated_boundary" : true
        })" );

        return default_parameters;
    }

    ShiftedBoundaryPointBasedUtility::ShiftedBoundaryPointBasedUtility(
        Model& rModel,
        Parameters ThisParameters)
    {
        // Validate input settings with defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Get volume model part and check that the volume model part has elements
        const std::string model_part_name = ThisParameters["model_part_name"].GetString();
        mpModelPart = &rModel.GetModelPart(model_part_name);
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements()) << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;

        // Get and check boundary conditions sub model part
        const std::string boundary_sub_model_part_name = ThisParameters["boundary_sub_model_part_name"].GetString();
        mpBoundarySubModelPart = &rModel.GetModelPart(boundary_sub_model_part_name);
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", mpBoundarySubModelPart->NumberOfNodes() != 0) << "Provided SBM conditions model part has nodes." << std::endl;
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", mpBoundarySubModelPart->NumberOfElements() != 0) << "Provided SBM conditions model part has elements." << std::endl;
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", mpBoundarySubModelPart->NumberOfConditions() != 0) << "Provided SBM conditions model part has conditions." << std::endl;

        // Get and check sub model part for skin discretization
        mSkinModelPartName = ThisParameters["skin_model_part_name"].GetString();
        mpSkinDiscSubModelPart = &rModel.GetModelPart(mSkinModelPartName);
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", mpSkinDiscSubModelPart->NumberOfElements() == 0) << "Provided SBM skin discretization model part has no elements." << std::endl;

        // Get and check sub model part for skin points
        mpSkinPointsSubModelPart = &rModel.GetModelPart(mSkinModelPartName + "Points");
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", mpSkinPointsSubModelPart->NumberOfNodes() != 0) << "Provided SBM skin points model part has nodes." << std::endl;
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", mpSkinPointsSubModelPart->NumberOfElements() != 0) << "Provided SBM skin points model part has elements." << std::endl;
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", mpSkinPointsSubModelPart->NumberOfConditions() != 0) << "Provided SBM skin points model part has conditions." << std::endl;

        // Set the order of the MLS extension operator used in the MLS shape functions utility
        mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();
        // If true, the basis is created such that the surrogate boundary gradient is kept
        const std::string ext_op_type = ThisParameters["extension_operator_type"].GetString();
        if (ext_op_type == "MLS") {
            mExtensionOperator = ExtensionOperator::MLS;
        } else if (ext_op_type == "RBF") {
            mExtensionOperator = ExtensionOperator::RBF;
        } else {
            KRATOS_ERROR << "Wrong 'extension_operator_type' provided. Only 'MLS' and 'RBF' are supported by point-based SBM utility." << std::endl;
        }

        // Set the shifted-boundary condition prototype to be used in the condition creation
        const std::string wall_condition_name = ThisParameters["boundary_wall_condition_name"].GetString();
        KRATOS_ERROR_IF(wall_condition_name == "") << "SBM boundary wall condition has not been provided." << std::endl;
        mpConditionPrototype = &KratosComponents<Condition>::Get(wall_condition_name);

        // Set which side of the skin model part is an enclosed area
        // If a side is given as enclosed, then the pressure of the first node will be set to zero.
        const std::string enclosed_area = ThisParameters["enclosed_area"].GetString();
        if (enclosed_area != "positive" && enclosed_area != "negative" && enclosed_area != "none") {
            KRATOS_ERROR << "Unknown 'enclosed_area' keyword given. 'positive' or 'negative' side or 'none' are supported by point based shifted boundary interface utility." << std::endl;
        }
        if (enclosed_area == "positive") {
            mPositiveSideIsEnclosed = true;
        } else if (enclosed_area == "negative") {
            mNegativeSideIsEnclosed = true;
        }
    }

    void ShiftedBoundaryPointBasedUtility::ResetFlags()
    {
        // Activate all elements and initialize flags to false
        // NOTE Resetting the SBM flags will eliminate previously immersed model parts except for their wall conditions (which remain in mpBoundarySubModelPart)!
        block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
            rNode.Set(ACTIVE, true);             // Nodes that belong to the elements to be assembled
            rNode.Set(SBM_BOUNDARY, false);      // Nodes that belong to the support clouds of the positive side
            rNode.Set(SBM_INTERFACE, false);     // Nodes that belong to the support clouds of the negative side
            rNode.Set(MODIFIED, false);          // Nodes that have been relocated because of a small distance to the skin geometry
        });
        block_for_each(mpModelPart->Elements(), [](ElementType& rElement){
            rElement.Set(ACTIVE, true);          // Elements to be assembled
            rElement.Set(SBM_BOUNDARY, false);   // Elements in which the skin geometry is located (true boundary gamma) - tessellated skin and skin integration points.
            rElement.Set(SBM_INTERFACE, false);  // Elements owning the surrogate boundary nodes
        });

        KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "Boundary and interface flags were reset and all elements and nodes re-activated." << std::endl;
    }

    void ShiftedBoundaryPointBasedUtility::FindElementsAtTessellatedBoundary(double Tolerance)
    {
        // Check that the skin model part has elements
        KRATOS_ERROR_IF_NOT(mpSkinDiscSubModelPart->NumberOfElements())
            << "There are no elements in skin model part (boundary) '" << mpSkinDiscSubModelPart->FullName() << "'." << std::endl;

        // Build AABBs of the discretized skin geometry (single-threaded)
        const std::size_t n_skin_elements = mpSkinDiscSubModelPart->NumberOfElements();
        std::vector<ShiftedBoundaryUtilityInternals::BBox> skin_boxes(n_skin_elements);
        std::vector<ElementType::Pointer> idx_to_skin_element_pointer(n_skin_elements);
        idx_to_skin_element_pointer.reserve(n_skin_elements);
        std::size_t idx = 0;
        for (auto& rElement : mpSkinDiscSubModelPart->Elements()) {
            skin_boxes[idx] = ShiftedBoundaryUtilityInternals::BBox(rElement.GetGeometry(), Tolerance);
            idx_to_skin_element_pointer[idx] = &rElement;
            idx += 1;
        }

        // Build BVH (AABB tree) of the discretized skin geometry
        mSkinBVH.build(skin_boxes);  // O(N_skin * log N_skin)
        mBVHIdxToSkinElementVector = idx_to_skin_element_pointer;

        // Reset the set of boundary elements
        const std::size_t n_intersected_elements_estimated = mpSkinDiscSubModelPart->NumberOfElements() * mpModelPart->GetProcessInfo()[DOMAIN_SIZE] * 2;
        mBoundaryElementsSet.clear();
        mBoundaryElementsSet.reserve(n_intersected_elements_estimated);

        const std::size_t num_threads = ParallelUtilities::GetNumThreads();
        std::vector< std::vector<ElementType::Pointer> > local_intersected_elements(num_threads);
        for (auto& vec : local_intersected_elements) vec.reserve(n_intersected_elements_estimated/num_threads);

        // Parallel search over background mesh for elements, which are intersected by the discretized skin geometry
        // NOTE that we already treat elements as intersected if one of their nodes coincides with the discretized skin geometry
        block_for_each(mpModelPart->Elements(), [&](ElementType& rElement){
            auto& r_background_geom = rElement.GetGeometry();
            const ShiftedBoundaryUtilityInternals::BBox element_box(r_background_geom, Tolerance);
            std::vector<std::size_t> candidates;
            mSkinBVH.query(element_box, candidates);

            // Check the background element for intersections with any skin element candidate (skin element which were found in proximity by the BVH query)
            // NOTE that "HasIntersection" of the geometry should already find coinciding nodes and faces!
            for (std::size_t idx : candidates) {
                const ElementType::Pointer p_skin_elem = mBVHIdxToSkinElementVector[idx];
                auto& r_skin_geom = p_skin_elem->GetGeometry();
                if (r_background_geom.HasIntersection(r_skin_geom)) {
                    local_intersected_elements[omp_get_thread_num()].push_back(&rElement);
                    break;
                }
            }
        });

        // Merge intersected background elements that were found by each thread
        for (auto& vec : local_intersected_elements) {
            mBoundaryElementsSet.insert(vec.begin(), vec.end());
        }
        // Free space, which was previously allocated by reserve
        mBoundaryElementsSet.rehash(mBoundaryElementsSet.size());
    }

    void ShiftedBoundaryPointBasedUtility::UpdateBoundaryElements(double Tolerance)
    {
        // Build AABBs of the discretized skin geometry (single-threaded)
        const std::size_t n_skin_elements = mpSkinDiscSubModelPart->NumberOfElements();
        std::vector<ShiftedBoundaryUtilityInternals::BBox> skin_boxes(n_skin_elements);
        std::vector<ElementType::Pointer> idx_to_skin_element_pointer(n_skin_elements);
        idx_to_skin_element_pointer.reserve(n_skin_elements);
        std::size_t idx = 0;
        for (auto& rElement : mpSkinDiscSubModelPart->Elements()) {
            skin_boxes[idx] = ShiftedBoundaryUtilityInternals::BBox(rElement.GetGeometry(), Tolerance);
            idx_to_skin_element_pointer[idx] = &rElement;
            idx += 1;
        }

        // Refit BVH (AABB tree) with updated skin geometry (instead of rebuilding the BVH from scratch) whenever the skin geometry moves
        mSkinBVH.refit(skin_boxes);
        mBVHIdxToSkinElementVector = idx_to_skin_element_pointer;

        // Reset the set of boundary elements
        const std::size_t n_intersected_elements_estimated = mpSkinDiscSubModelPart->NumberOfElements() * mpModelPart->GetProcessInfo()[DOMAIN_SIZE] * 2;
        mBoundaryElementsSet.clear();
        mBoundaryElementsSet.reserve(n_intersected_elements_estimated);

        const std::size_t num_threads = ParallelUtilities::GetNumThreads();
        std::vector< std::vector<ElementType::Pointer> > local_intersected_elements(num_threads);
        for (auto& vec : local_intersected_elements) vec.reserve(n_intersected_elements_estimated/num_threads);

        // Parallel search over background mesh candidates for elements, which are intersected by the discretized skin geometry
        // NOTE that we already treat elements as intersected if one of their nodes coincides with the discretized skin geometry
        block_for_each(mpModelPart->Elements(), [&](ElementType& rElement){
            auto& r_background_geom = rElement.GetGeometry();
            const ShiftedBoundaryUtilityInternals::BBox element_box(r_background_geom, Tolerance);
            std::vector<std::size_t> candidates;
            mSkinBVH.query(element_box, candidates);

            // Check the background element for intersections with any skin element candidate (skin element which were found in proximity by the skin BVH query)
            // NOTE that "HasIntersection" of the geometry should already find coinciding nodes and faces!
            for (std::size_t idx : candidates) {
                const ElementType::Pointer p_skin_elem = mBVHIdxToSkinElementVector[idx];
                auto& r_skin_geom = p_skin_elem->GetGeometry();
                if (r_background_geom.HasIntersection(r_skin_geom)) {
                    local_intersected_elements[omp_get_thread_num()].push_back(&rElement);
                    break;
                }
            }
        });

        // Merge intersected background elements that were found by each thread
        for (auto& vec : local_intersected_elements) {
            mBoundaryElementsSet.insert(vec.begin(), vec.end());
        }
        // Free space, which was previously allocated by reserve
        mBoundaryElementsSet.rehash(mBoundaryElementsSet.size());
    }

    void ShiftedBoundaryPointBasedUtility::MapSkinPointsToElements()
    {
        mSkinPointsMap.clear();
        // Map skin points (boundary) to elements of volume model part
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        switch (n_dim) {
            case 2:
                MapSkinPointsToElementsTemplated<2>(mSkinPointsMap);
                break;
            case 3:
                MapSkinPointsToElementsTemplated<3>(mSkinPointsMap);
                break;
            default:
                KRATOS_ERROR << "Wrong domain size.";
        }
    }

    void ShiftedBoundaryPointBasedUtility::FlagBoundaryElements()
    {
        // Mark elements as SBM_BOUNDARY (Gamma) in which the tessellated skin geometry and/ or skin integration points were located
        // NOTE that SBM_BOUNDARY elements will get deactivated and that the BC inside boundary elements will be applied by means of the extension operators
        std::for_each(mBoundaryElementsSet.begin(), mBoundaryElementsSet.end(), [](ElementType::Pointer pElement) {
            pElement->Set(SBM_BOUNDARY, true);
        });
    }

    void ShiftedBoundaryPointBasedUtility::FlagInterfaceElements()
    {
        // Find the surrogate boundary elements and mark them as SBM_INTERFACE (gamma_tilde)
        // Note that we rely on the fact that the neighbors are sorted according to the faces
        LockObject mutex;
        block_for_each(mpModelPart->Elements(), [&mutex](ElementType& rElement){
            if (rElement.Is(SBM_BOUNDARY)) {
                //TODO for laplacian testing
                rElement.Set(BOUNDARY, true);
                const std::size_t n_faces = rElement.GetGeometry().FacesNumber();
                auto& r_neigh_elems = rElement.GetValue(NEIGHBOUR_ELEMENTS);
                for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
                    // If neighbour corresponding to the current face is not also SBM_BOUNDARY, it means that the current face is surrogate boundary (SBM_INTERFACE)
                    // Flag the current neighbour owning the surrogate face as SBM_INTERFACE
                    // The nodes will be flagged if required (MLS basis) when creating the cloud
                    auto p_neigh_elem = r_neigh_elems(i_face).get();
                    if (p_neigh_elem != nullptr) {
                        if (!p_neigh_elem->Is(SBM_BOUNDARY)) {
                            {
                                std::scoped_lock<LockObject> lock(mutex);
                                p_neigh_elem->Set(SBM_INTERFACE, true);
                                //TODO for laplacian testing
                                p_neigh_elem->Set(INTERFACE, true);
                            }
                        }
                    }
                }
            }
        });
        KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "Interface flags were set."  << std::endl;
    }

    void ShiftedBoundaryPointBasedUtility::DeactivateElementsAndNodes(
        const bool DeactivateUnstableClusters
    )
    {
        // Deactivate elements in which the (true) boundary is located
        block_for_each(mpModelPart->Elements(), [](ElementType& rElement){
            if (rElement.Is(SBM_BOUNDARY)) {
                rElement.Set(ACTIVE, false);
            }
        });
        KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "Boundary elements were deactivated." << std::endl;

        if (DeactivateUnstableClusters) {
            FindAndDeactivateUnstableClusters();
            KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "Unstable Clusters were deactivated." << std::endl;
        }

        // Deactivate nodes that do not belong to any active element (anymore)
        block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
            // Check whether any of the elements surrounding the node is active
            const auto& r_neigh_elems = rNode.GetValue(NEIGHBOUR_ELEMENTS);
            bool active_neighbor_found = false;
            for (std::size_t i_neigh = 0; i_neigh < r_neigh_elems.size(); ++i_neigh) {
                auto p_neigh = r_neigh_elems(i_neigh).get();
                if (p_neigh != nullptr) {
                    // Continue with next node if the current node is part of an active element
                    if (p_neigh->Is(ACTIVE)) {
                        active_neighbor_found = true;
                        break;
                    }
                }
            }
            // If there is no active element around the current node, we will deactivate it
            if (!active_neighbor_found) {
                rNode.Set(ACTIVE, false);
                // Set all dofs to zero
                for (auto& p_dof : rNode.GetDofs()) {
                    const auto& r_variable = p_dof->GetVariable();
                    rNode.FastGetSolutionStepValue(static_cast<const Variable<double>&>(r_variable)) = 0.0;
                }
            }
        });

        KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "Single nodes were deactivated." << std::endl;
    }

    void ShiftedBoundaryPointBasedUtility::CalculateAndAddSkinIntegrationPointConditions()
    {
        // Nodes are reset here to not be disturbed by other skin model parts - TODO rethink this
        for (auto& p_elem : mBoundaryElementsSet) {
            for (auto& r_node : p_elem->GetGeometry()) {
                r_node.Set(SBM_BOUNDARY, false);
                r_node.Set(SBM_INTERFACE, false);
            }
        }

        // Iterate over the elements which contain skin integration points to create a vector for each element defining both sides.
        // The resulting vector is as long as the number of nodes of the element and a positive value stands for the positive side of the boundary, a negative one for the negative side.
        // Also store the average position and average normal of the skin points located in the element
        // Nodes on the positive side will be declared SBM_BOUNDARY and nodes on the negative side SBM_INTERFACE which is used in the creation of the support clouds
        mSidesVectorMap.clear();
        AverageSkinToElementsMapType avg_skin_map;
        SetSidesForSkinPointElements(mSkinPointsMap, mSidesVectorMap, avg_skin_map);
        //KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "Sides vectors and skin normals were set." << std::endl;

        // Iterate over the elements which contain skin integration points to create an extension basis for each node of the element (MLS shape functions values for support cloud of node)
        // NOTE that no extension bases will be calculated and added for a node for which not a sufficient number of support nodes were found
        mExtensionOperatorMap.clear();
        SetExtensionForSkinPointElements(mSidesVectorMap, avg_skin_map, mExtensionOperatorMap);
        //KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "Extension operators were set." << std::endl;

        // Set the pressure of the first node of an enclosed volume to zero if one side is enclosed.
        auto skin_pt_element_iter = mSkinPointsMap.begin();
        bool enclosed_pressure_is_set = false;
        if (mPositiveSideIsEnclosed || mNegativeSideIsEnclosed) {
            while (!enclosed_pressure_is_set && skin_pt_element_iter != mSkinPointsMap.end()) {
                auto p_elem = skin_pt_element_iter->first;
                const array_1d<double,3> avg_skin_position = avg_skin_map[p_elem].first;
                const array_1d<double,3> avg_skin_normal = avg_skin_map[p_elem].second;
                enclosed_pressure_is_set = FixPressureOfEnclosedNode(*p_elem, mSidesVectorMap[p_elem], avg_skin_position, avg_skin_normal);
                skin_pt_element_iter++;
            }
        }

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        const auto p_element_size_func = GetElementSizeFunction(mpModelPart->ElementsBegin()->GetGeometry());

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Create the interface conditions
        //TODO: THIS CAN BE PARALLEL (WE JUST NEED TO MAKE CRITICAL THE CONDITION ID UPDATE)? And adding the condition????
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        std::size_t n_skin_points = 0;
        std::size_t n_skin_pt_conditions_added_pos = 0;
        std::size_t n_skin_pt_conditions_added_neg = 0;

        for (const auto& [p_element, skin_points_data_vector]: mSkinPointsMap) {
            const auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            n_skin_points += skin_points_data_vector.size();

            // For each side of the boundary separately (positive and negative side of gamma), create a pointer vector with all the nodes that affect that side of the current element
            // To be used in the creation of the condition. Positive side refers to adding the positive side's nodes of the element and the negative node's support cloud nodes.
            // NOTE that the obtained clouds are sorted by id to properly get the extension operator data //TODO: necessary? ids are compared anyway?!
            PointerVector<NodeType> cloud_nodes_vector_pos;
            PointerVector<NodeType> cloud_nodes_vector_neg;
            CreateCloudVectorsForSkinPointElement(*p_element, mSidesVectorMap[p_element], mExtensionOperatorMap, cloud_nodes_vector_pos, cloud_nodes_vector_neg);

            // Calculate parent element size for the SBM BC imposition
            const double h = p_element_size_func(r_geom);

            // Get vector defining positive and negative side of the boundary
            const auto& sides_vector = mSidesVectorMap[p_element];

            // Iterate over the element's skin points adding a positive side and a negative side condition for each skin point
            for (std::size_t i_skin_pt = 0; i_skin_pt < skin_points_data_vector.size(); ++i_skin_pt) {
                // Get the skin point's position and area normal (integration point of the boundary)
                const auto skin_pt_data = skin_points_data_vector[i_skin_pt];
                const array_1d<double,3> skin_pt_position = std::get<0>(skin_pt_data);
                const array_1d<double,3> skin_pt_area_normal = std::get<1>(skin_pt_data);

                // Get the element's shape function values and derivatives at the skin/ integration point
                Vector skin_pt_N(n_nodes);
                Matrix skin_pt_DN_DX = ZeroMatrix(n_nodes, n_dim);
                GetDataForSkinPointInElement(*p_element, skin_pt_position, skin_pt_N, skin_pt_DN_DX);

                // Add skin pt. condition for positive side of boundary - using support cloud data for negative nodes
                // NOTE that the boundary normal is negative in order to point outward (from positive to negative side),
                // because positive side is where dot product of vector to node with average normal is positive
                n_skin_pt_conditions_added_pos += AddSkinPointCondition(*p_element, sides_vector, h, skin_pt_position, -skin_pt_area_normal,
                mExtensionOperatorMap, cloud_nodes_vector_pos, skin_pt_N, skin_pt_DN_DX, max_cond_id, /*ConsiderPositiveSide=*/true);

                // Add skin pt. condition for negative side of boundary - using support cloud data for positive nodes
                // NOTE that boundary normal is pointing outward (from negative to positive side)
                n_skin_pt_conditions_added_neg += AddSkinPointCondition(*p_element, sides_vector, h, skin_pt_position, skin_pt_area_normal,
                mExtensionOperatorMap, cloud_nodes_vector_neg, skin_pt_N, skin_pt_DN_DX, max_cond_id, /*ConsiderPositiveSide=*/false);
            }
        }
        if (n_skin_pt_conditions_added_pos != n_skin_points) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedUtility") << "Integration point conditions were NOT successfully added for the positive side of "
                << n_skin_points-n_skin_pt_conditions_added_pos << " skin points." << std::endl;
        }
        if (n_skin_pt_conditions_added_neg != n_skin_points) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedUtility") << "Integration point conditions were NOT successfully added for the negative side of "
                << n_skin_points-n_skin_pt_conditions_added_neg << " skin points." << std::endl;
        }
        KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "'" << mSkinModelPartName << "' skin point conditions were added." << std::endl;
    }

    void ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPoints()
    {
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        switch (n_dim) {
            case 2:
                CalculateVariablesAtSkinPointsTemplated<2>();
                break;
            case 3:
                CalculateVariablesAtSkinPointsTemplated<3>();
                break;
            default:
                KRATOS_ERROR << "Wrong domain size.";
        }
    }

    void ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPointsAndNodes()
    {
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        switch (n_dim) {
            case 2:
                CalculateVariablesAtSkinPointsAndNodesTemplated<2>();
                break;
            case 3:
                CalculateVariablesAtSkinPointsAndNodesTemplated<3>();
                break;
            default:
                KRATOS_ERROR << "Wrong domain size.";
        }
    }

    //----------------------------------------------------------------
    //     PROTECTED METHODS
    //----------------------------------------------------------------

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedUtility::MapSkinPointsToElementsTemplated(
        SkinPointsToElementsMapType& rSkinPointsMap)
    {
        // Check that the skin model part has elements
        KRATOS_ERROR_IF_NOT(mpSkinPointsSubModelPart->NumberOfNodes())
            << "There are no nodes in skin points model part '" << mpSkinPointsSubModelPart->FullName() << "'." << std::endl;

        // Check that BOUNDARY elements were already founnd by FindElementsAtTessellatedBoundary()
        const std::size_t n_boundary_elements_found = mBoundaryElementsSet.size();
        KRATOS_ERROR_IF_NOT(n_boundary_elements_found)
            << "There are no elements in mBoundaryElementsSet. Maybe FindElementsAtTessellatedBoundary() was not called." << std::endl;

        // Create a set of background mesh candidates for the skin point search
        ElementsSetType bg_candidates;
        bg_candidates.reserve(n_boundary_elements_found * 2 * TDim);
        // Add elements, in which the skin boundary was already found and their immediate neighbors to the search candidates
        const std::size_t n_faces = mBoundaryElementsSet.begin()->get()->GetGeometry().FacesNumber();
        for (auto p_elem : mBoundaryElementsSet) {
            bg_candidates.emplace(p_elem);
            auto& r_neigh_elems = p_elem->GetValue(NEIGHBOUR_ELEMENTS);
            for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
                auto p_neigh_elem = r_neigh_elems(i_face).get();
                if (p_neigh_elem != nullptr) {
                    bg_candidates.emplace(p_neigh_elem);
                }
            }
        }

        const std::size_t n_candidates = bg_candidates.size();
        std::vector<ShiftedBoundaryUtilityInternals::BBox> candidate_boxes(n_candidates);
        std::vector<ElementType::Pointer> idx_to_candidate_element_pointer(n_candidates);
        idx_to_candidate_element_pointer.reserve(n_candidates);
        std::size_t idx = 0;
        for (auto p_elem : bg_candidates) {
            candidate_boxes[idx] = ShiftedBoundaryUtilityInternals::BBox(p_elem->GetGeometry(), 1e-10);
            idx_to_candidate_element_pointer[idx] = p_elem;
            idx += 1;
        }
        ShiftedBoundaryUtilityInternals::ElementBVH candidate_bvh;
        candidate_bvh.build(candidate_boxes);

        // Create vectors for skin point data and sums for parallelization
        const std::size_t num_threads = ParallelUtilities::GetNumThreads();
        const std::size_t n_skin_pts_expected = mpSkinPointsSubModelPart->NumberOfNodes();
        const std::size_t n_local_skin_pts_expected = n_skin_pts_expected/num_threads;
        std::vector< std::vector<Element::Pointer> > local_skin_point_located_elements(num_threads);
        for (auto& vec : local_skin_point_located_elements) vec.reserve(n_local_skin_pts_expected);
        std::vector< std::vector<array_1d<double,3>> > local_skin_point_positions(num_threads);
        for (auto& vec : local_skin_point_positions) vec.reserve(n_local_skin_pts_expected);
        std::vector< std::vector<array_1d<double,3>> > local_skin_point_normals(num_threads);
        for (auto& vec : local_skin_point_normals) vec.reserve(n_local_skin_pts_expected);
        std::vector< std::vector<std::size_t> > local_skin_point_ids(num_threads);
        for (auto& vec : local_skin_point_ids) vec.reserve(n_local_skin_pts_expected);
        std::vector< std::size_t > local_n_skin_points_not_found(num_threads, 0);
        std::vector< std::size_t > local_n_skin_points_found(num_threads, 0);

        // Search the skin points in the volume mesh elements
        block_for_each(mpSkinPointsSubModelPart->Nodes(), [&](NodeType& rSkinPoint){
            // Get position of skin point
            array_1d<double,3> skin_pt_position = rSkinPoint.Coordinates();
            // Get the normal at the skin point
            // NOTE that we assume here that the norm/ length of the normal is a measure of the area/ integration point weight
            array_1d<double,3> skin_pt_area_normal = rSkinPoint.GetValue(NORMAL);
            // Search for the skin point in the volume mesh candidates to get the element containing the point
            Element::Pointer p_element = nullptr;
            const bool is_found = LocatePoint(candidate_bvh, idx_to_candidate_element_pointer, skin_pt_position, p_element);

            // Add data to local vectors
            const std::size_t thread_num = omp_get_thread_num();
            if (is_found) {
                local_skin_point_located_elements[thread_num].emplace_back(p_element);
                local_skin_point_positions[thread_num].emplace_back(skin_pt_position);
                local_skin_point_normals[thread_num].emplace_back(skin_pt_area_normal);
                local_skin_point_ids[thread_num].emplace_back(rSkinPoint.Id());
                local_n_skin_points_found[thread_num]++;
            } else {
                local_n_skin_points_not_found[thread_num]++;
            }
        });

        // Create and reserve skin point data vectors
        std::vector<Element::Pointer> skin_point_located_elements;
        skin_point_located_elements.reserve(n_skin_pts_expected);
        std::vector<array_1d<double,3>> skin_point_positions;
        skin_point_positions.reserve(n_skin_pts_expected);
        std::vector<array_1d<double,3>> skin_point_normals;
        skin_point_normals.reserve(n_skin_pts_expected);
        std::vector<std::size_t> skin_point_ids;
        skin_point_ids.reserve(n_skin_pts_expected);

        // Merge skin point data and sums
        std::size_t n_skin_points_found = 0;
        std::size_t n_skin_points_not_found = 0;
        for (std::size_t i = 0; i < num_threads; ++i) {
            skin_point_located_elements.insert(skin_point_located_elements.end(),
                local_skin_point_located_elements[i].begin(), local_skin_point_located_elements[i].end());
            skin_point_positions.insert(skin_point_positions.end(),
                local_skin_point_positions[i].begin(),  local_skin_point_positions[i].end());
            skin_point_normals.insert(skin_point_normals.end(),
                local_skin_point_normals[i].begin(), local_skin_point_normals[i].end());
            skin_point_ids.insert(skin_point_ids.end(),
                local_skin_point_ids[i].begin(), local_skin_point_ids[i].end());
            n_skin_points_found += local_n_skin_points_found[i];
            n_skin_points_not_found += local_n_skin_points_not_found[i];
        }

        if (n_skin_points_not_found > 0) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedUtility")
                << n_skin_points_not_found << " skin points have not been found in any volume model part element." << std::endl;
        }

        // Store skin points and their data for each model part element with skin points in the skin points map
        for (std::size_t i_skin_pt = 0; i_skin_pt < skin_point_located_elements.size(); ++i_skin_pt) {
            auto p_element =  skin_point_located_elements[i_skin_pt];
            auto& skin_pt_position = skin_point_positions[i_skin_pt];
            auto& skin_pt_normal = skin_point_normals[i_skin_pt];
            auto& skin_pt_node_id = skin_point_ids[i_skin_pt];

            // Add skin point data to the element it was found in
            auto& skin_points_data_vector = rSkinPointsMap[p_element]; // automatically inserts if not present
            skin_points_data_vector.emplace_back(skin_pt_position, skin_pt_normal, skin_pt_node_id);
        }

        // Add elements in which skin points were found to the set of boundary elements
        for (auto& [p_elem, skin_pt_data] : mSkinPointsMap) {
            mBoundaryElementsSet.insert(p_elem);
        }

        KRATOS_INFO("ShiftedBoundaryPointBasedUtility") << "'" << mSkinModelPartName << "' skin points ("
            << n_skin_points_found << ") were mapped to volume mesh elements (" << rSkinPointsMap.size() << ")." << std::endl;
    }

    bool ShiftedBoundaryPointBasedUtility::LocatePoint(
        ShiftedBoundaryUtilityInternals::ElementBVH& rCandidatesBvh,
        std::vector<ElementType::Pointer>& rIdxToPointer,
        array_1d<double,3>& rPointCoords,
        ElementType::Pointer& pElement)
    {
        ShiftedBoundaryUtilityInternals::BBox pt_box;

        pt_box.min[0] = pt_box.max[0] = rPointCoords[0];
        pt_box.min[1] = pt_box.max[1] = rPointCoords[1];
        pt_box.min[2] = pt_box.max[2] = rPointCoords[2];

        std::vector<std::size_t> bbox_candidates;
        bbox_candidates.reserve(20);
        rCandidatesBvh.query(pt_box, bbox_candidates);  // read-only query of the BVH - safe for parallel execution

        Element::GeometryType::CoordinatesArrayType coords;
        bool is_inside = false;
        for (auto elem_idx : bbox_candidates) {
            auto p_elem = rIdxToPointer[elem_idx];
            is_inside = p_elem->GetGeometry().IsInside(rPointCoords, coords);
            if(is_inside) {
                pElement = p_elem;
                break;
            }
        }

        return is_inside;
    }

    void ShiftedBoundaryPointBasedUtility::FindAndDeactivateUnstableClusters()
    {
        // Get clusters of connected elements (clusters of ACTIVE elements that are connected via their edges) in the domain
        std::vector<std::vector<ElementType::Pointer>> clusters = FindClusters();

        //Find biggest cluster
        std::size_t biggest_cluster_id = 0;
        std::size_t biggest_cluster_n = 0;
        for (std::size_t i = 0; i < clusters.size(); ++i) {
            if (clusters[i].size() > biggest_cluster_n) {
                biggest_cluster_id = i;
                biggest_cluster_n = clusters[i].size();
            }
        }

        // Deactivate the elements of all clusters except for the elements of the biggest cluster
        //TODO deactivate only clusters that do not have a fixed DOF?
        for (std::size_t i = 0; i < clusters.size(); ++i) {
            if (i != biggest_cluster_id) {
                for (auto& p_elem : clusters[i]) {
                    p_elem->Set(ACTIVE, false);
                }
            }
        }
    }

    std::vector<std::vector<ShiftedBoundaryPointBasedUtility::ElementType::Pointer>> ShiftedBoundaryPointBasedUtility::FindClusters()
    {
        struct UnionFind {
            std::vector<std::size_t> parent;
            std::vector<std::size_t> rank;

            UnionFind(std::size_t max_id) : parent(max_id + 1), rank(max_id + 1, 0) {
                // Initialize the parent as itself.
                for (std::size_t i = 0; i <= max_id; ++i) {
                    parent[i] = i;
                }
            }

            std::size_t find_parent(std::size_t id) {
                // If the parent points to another Id assign that IDs parent recursively until the ID is its own parent (path compression).
                if (parent[id] != id) {
                    parent[id] = find_parent(parent[id]);
                }
                return parent[id];
            }

            void unite_parents(std::size_t id_A, std::size_t id_B) {
                std::size_t root_A = find_parent(id_A);
                std::size_t root_B = find_parent(id_B);

                // Do nothing if parents are already the same.
                if (root_A == root_B) return;

                // Assign higher ranked parent as lower ranked one as well.
                if (rank[root_A] < rank[root_B]) {
                    parent[root_A] = root_B;
                } else if (rank[root_A] > rank[root_B]) {
                    parent[root_B] = root_A;
                } else {
                    parent[root_A] = root_B;
                    ++rank[root_B];
                }
            }
        };

        // Find connected elements based on the adjacency graph of ACTIVE elements
        EdgesVectorType all_active_edges = GetActiveAdjacencyGraph();

        // Create a mapping from element IDs to consecutive indices for the Union-Find data structure
        std::unordered_map<std::size_t, std::size_t> element_id_to_index;
        std::size_t index = 0;
        for (const auto& rElement : mpModelPart->Elements()) {
            element_id_to_index[rElement.Id()] = index++;
        }

        // Create Union-Find data structure to unite elements that are connected via an edge in the adjacency graph
        const std::size_t num_elem = element_id_to_index.size();
        UnionFind union_find(num_elem);

        // Find roots of elements in chunks
        // NOTE that it is faster sequentially, instead of dividing domain into chunks for parallelization
        for (const auto& edge : all_active_edges) {
            const std::size_t idx_A = element_id_to_index[edge.first];
            const std::size_t idx_B = element_id_to_index[edge.second];
            union_find.unite_parents(idx_A, idx_B);
        }

        // Define clusters based on the found parents
        std::vector<std::vector<ElementType::Pointer>> clusters;
        std::unordered_map<std::size_t, std::size_t> root_to_cluster;
        for (auto& rElement : mpModelPart->Elements()) {
            if (rElement.Is(ACTIVE)) {
                std::size_t root = union_find.find_parent(element_id_to_index[rElement.Id()]);

                // Add the element to the corresponding cluster
                if (root_to_cluster.count(root) == 0) {
                    root_to_cluster[root] = clusters.size();
                    clusters.emplace_back();
                }
                clusters[root_to_cluster[root]].push_back(&rElement);
            }
        }

        return clusters;
    }

    ShiftedBoundaryPointBasedUtility::EdgesVectorType ShiftedBoundaryPointBasedUtility::GetActiveAdjacencyGraph()
    {
        const std::size_t num_threads = ParallelUtilities::GetNumThreads();
        const std::size_t num_elem = mpModelPart->NumberOfElements();

        // Create a vector for all edges found by all threads
        std::vector<EdgesVectorType> thread_edges(num_threads);
        for (auto& edges : thread_edges) {
            edges.reserve(num_elem * 3 / num_threads); // Rough estimate of edges per thread
        }

        // Loop over all elements and add an edge for each neighboring element if both elements are ACTIVE and if current element ID is smaller than neighbor element ID (undirected graph)
        block_for_each(mpModelPart->Elements(), [&](ElementType& rElement){
            auto& edges = thread_edges[omp_get_thread_num()];
            if (rElement.Is(ACTIVE)) {
                const std::size_t n_faces = rElement.GetGeometry().FacesNumber();
                auto& r_neigh_elems = rElement.GetValue(NEIGHBOUR_ELEMENTS);
                for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
                    // If neighbour corresponding to the current face is ACTIVE, add edge between the current element and the neighbour element to the graph.
                    auto p_neigh_elem = r_neigh_elems(i_face).get();
                    if (p_neigh_elem != nullptr) {
                        if (p_neigh_elem->Is(ACTIVE) && rElement.Id() < p_neigh_elem->Id()) {
                            edges.emplace_back(rElement.Id(), p_neigh_elem->Id());
                        }
                    }
                }
            }
        });

        EdgesVectorType all_edges;
        all_edges.reserve(num_elem * 3);
        for (auto& edges : thread_edges) {
            all_edges.insert(all_edges.end(), edges.begin(), edges.end());
        }

        return all_edges;
    }

    void ShiftedBoundaryPointBasedUtility::SetSidesForSkinPointElements(
        const SkinPointsToElementsMapType& rSkinPointsMap,
        SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap)
    {
        //TODO SPEED UP

        // Set DISTANCE values for all nodes to zero as variable will be used to have a majority vote on the definition of the positive and negative side
        //TODO faster than looping through the nodes of elements with skin points without parallelization?
        VariableUtils().SetVariable(DISTANCE, 0.0, mpModelPart->Nodes());

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        const auto p_element_size_func = GetElementSizeFunction(mpModelPart->ElementsBegin()->GetGeometry());

        LockObject mutex_1;
        LockObject mutex_2;
        //TODO is parallelization really faster here?
        std::for_each(rSkinPointsMap.begin(), rSkinPointsMap.end(), [&rAvgSkinMap, &mutex_1, &mutex_2, &p_element_size_func](const std::pair<ElementType::Pointer, SkinPointsDataVectorType>& rKeyData){
            const auto p_element = rKeyData.first;
            const auto& skin_points_data_vector = rKeyData.second;

            // Get access to the nodes of the element and an estimate of the element's size
            auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            const double h = p_element_size_func(r_geom);

            // Initialize average position and average normal of the element's skin points
            std::size_t n_skin_points = 0;
            array_1d<double,3> avg_position(3, 0.0);
            array_1d<double,3> avg_normal(3, 0.0);

            // Loop over the skin points located inside the element
            for (const auto& skin_pt_data: skin_points_data_vector) {
                const auto& skin_pt_position = std::get<0>(skin_pt_data);
                const auto& skin_pt_area_normal = std::get<1>(skin_pt_data);
                const double skin_pt_area = norm_2(skin_pt_area_normal);

                // Compute the dot product for each skin point and node of the element between a vector from the skin position to the node and the skin point's normal
                // NOTE that for a positive dot product the node is saved as being on the positive side of the boundary, negative dot product equals negative side
                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                    auto& r_node = r_geom[i_node];
                    const array_1d<double,3> skin_pt_to_node = r_node.Coordinates() - skin_pt_position;  // Should never be zero because of node modifications
                    const double dot_product = inner_prod(skin_pt_to_node, skin_pt_area_normal);
                    double& side_voting = r_node.FastGetSolutionStepValue(DISTANCE);
                    {
                        const double vote_weighting = skin_pt_area * h / std::pow(norm_2(skin_pt_to_node),3);
                        std::scoped_lock<LockObject> lock(mutex_1);
                        if (dot_product > 0.0) {
                            side_voting += vote_weighting;
                        } else {
                            side_voting -= vote_weighting;
                        }
                    }
                }

                // Update average position and average normal
                n_skin_points++;
                avg_position += skin_pt_position;
                avg_normal += skin_pt_area_normal;
            }

            // Calculate and store the average position and average normal of the skin points located in the element
            avg_normal /= norm_2(avg_normal);
            avg_position /= n_skin_points;
            const auto avg_position_and_normal = std::make_pair(avg_position, avg_normal);
            {
                std::scoped_lock<LockObject> lock(mutex_2);
                rAvgSkinMap.insert(std::make_pair(p_element, avg_position_and_normal));
            }
        });

        // Decide on positive and negative side of an element based on the voting of all skin points in the surrounding elements
        for (const auto& [p_element, skin_points_data_vector]: rSkinPointsMap) {
            auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();

            // Store a vector deciding on the positive and negative side of the element's nodes
            // NOTE that the positive side of the boundary equals a positive inwards skin normal, negative dot product equals a negative inward skin normal
            // NOTE that it is necessary to define the side of points of gamma_tilde here for the search of the support clouds afterwards (SetLateralSupportCloud)
            // NOTE that this will cause troubles if inverted elements exist as opposed to having a local definition of sides
            // NOTE that it is necessary here to set the other side's flag to false because it might have been set to true by another skin geometry embedded previously.
            Vector sides_vector(n_nodes);
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                auto& r_node = r_geom[i_node];
                const double& side_voting = r_node.FastGetSolutionStepValue(DISTANCE);
                if (side_voting > 0.0) {
                    sides_vector[i_node] =  1.0;
                    r_node.Set(SBM_BOUNDARY, true);
                } else {
                    sides_vector[i_node] = -1.0;
                    r_node.Set(SBM_INTERFACE, true);
                }
            }
            rSidesVectorMap.insert(std::make_pair(p_element, sides_vector));
        }
    }

    void ShiftedBoundaryPointBasedUtility::SetExtensionForSkinPointElements(
        const SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap,
        NodesCloudMapType& rExtensionOperatorMap)
    {
        // Get the extension operator shape functions function
        auto p_meshless_sh_func = mExtensionOperator == ExtensionOperator::MLS ? GetMLSShapeFunctionsFunction() : GetRBFShapeFunctionsFunction();

        // Get support node clouds for all nodes of all elements which contain skin integration points and calculate their extension operators
        // NOTE that only extension operators are calculated and added to the map if a sufficient number of support nodes was found
        //TODO make parallel
        for (const auto& [p_element, sides_vector]: rSidesVectorMap) {
            const auto& r_geom = p_element->GetGeometry();

            // Get averaged position and normal of the skin points located inside the element
            auto avg_position_and_normal = rAvgSkinMap[p_element];
            const array_1d<double,3> avg_position = std::get<0>(avg_position_and_normal);
            const array_1d<double,3> avg_normal = std::get<1>(avg_position_and_normal);

            // Get support node cloud and calculate the extension operator for all nodes of the element for which it has not been calculated yet
            const std::size_t n_nodes = r_geom.PointsNumber();
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                const auto& p_node = r_geom(i_node);

                // Check if extension operator already has been calculated for the current node
                const std::size_t found_in_map = rExtensionOperatorMap.count(p_node);
                if (!found_in_map) {

                    // Initialize the storage for the support/ cloud nodes and their coordinates
                    Matrix cloud_nodes_coordinates;
                    PointerVector<NodeType> cloud_nodes;

                    // Get support cloud for given node
                    if (sides_vector[i_node] < 0.0) {
                        // Use and declare SBM_BOUNDARY nodes on the positive side for the support cloud of a node on the negative side
                        SetLateralSupportCloud(p_node, avg_position,  avg_normal, cloud_nodes, cloud_nodes_coordinates, SBM_BOUNDARY);
                    } else {
                        // Use and declare SBM_INTERFACE nodes on the negative side for the support cloud of a node on the positive side
                        SetLateralSupportCloud(p_node, avg_position, -avg_normal, cloud_nodes, cloud_nodes_coordinates, SBM_INTERFACE);
                    }

                    // Continue if the number of support nodes is sufficient for the calculation of the extension operator
                    const std::size_t n_cloud_nodes = cloud_nodes.size();
                    if (n_cloud_nodes >= GetRequiredNumberOfPoints()) {

                        // Calculate the extension basis in the current node (MLS shape functions)
                        Vector N_container;
                        const array_1d<double,3> r_coords = p_node->Coordinates();
                        const double kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
                        p_meshless_sh_func(cloud_nodes_coordinates, r_coords, kernel_rad, N_container);

                        // if (sides_vector[i_node] < 0.0) {
                        //     Vector N_container_aux;
                        //     Matrix DNDX_container;
                        //     MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,1>(cloud_nodes_coordinates, r_coords, kernel_rad, N_container_aux, DNDX_container);
                        //     KRATOS_WATCH(DNDX_container);
                        // }

                        // Save the extension operator nodal data to the extension operator map
                        CloudDataVectorType cloud_data_vector(n_cloud_nodes);
                        for (std::size_t i_cl_nod = 0; i_cl_nod < n_cloud_nodes; ++i_cl_nod) {
                            auto p_cl_node = cloud_nodes(i_cl_nod);
                            auto i_data = std::make_pair(p_cl_node, N_container[i_cl_nod]);
                            cloud_data_vector[i_cl_nod] = i_data;
                        }
                        const auto ext_op_key_data = std::make_pair(p_node, cloud_data_vector);
                        //TODO make this threadsafe for parallelization
                        rExtensionOperatorMap.insert(ext_op_key_data);
                    // } else {
                    //     KRATOS_WARNING("ShiftedBoundaryPointBasedUtility")
                    //     << "No enough support nodes were found for node " << p_node->Id() << ". Extension basis can not be calculated." << std::endl;
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedUtility::SetLateralSupportCloud(
        const NodeType::Pointer pOtherSideNode,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates,
        const Kratos::Flags& rSearchSideFlag)
    {
        //TODO SPEED UP

        // Find the support cloud of nodes on the search side (other side as the given node)
        // NOTE that we use an unordered_set to ensure that these are unique
        // NOTE that we check the order of the MLS interpolation to add nodes from enough layers
        NodesSetType aux_set;
        std::vector<NodeType::Pointer> cur_layer_nodes;
        std::vector<NodeType::Pointer> prev_layer_nodes;
        const std::size_t n_layers = mMLSExtensionOperatorOrder + 1;

        // Find elemental neighbors of the given node and add their nodes to the cloud nodes set if they are located on the search side
        // This is the first layer of sampling/ support points
        // NOTE that the sides of the first layer of nodes at gamma_tilde already need to be defined in their flags (done by SetSidesForSkinPointElements)
        // NOTE that taking the nodes of neighboring elements is the same as adding the nodal neighbors directly for triangles and tetrahedra
        //TODO add neighboring nodes directly? for tetra and hex elements?
        auto& r_elem_neigh_vect = pOtherSideNode->GetValue(NEIGHBOUR_ELEMENTS);
        for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
            auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
            if (p_elem_neigh != nullptr) {
                const auto& r_geom = p_elem_neigh->GetGeometry();
                for (std::size_t i_neigh_node = 0; i_neigh_node < r_geom.PointsNumber(); ++i_neigh_node) {
                    NodeType::Pointer p_neigh = r_geom(i_neigh_node);

                    // Add node of neighboring element to the support node set if it is active and located on the search side
                    if (p_neigh->Is(ACTIVE) && p_neigh->Is(rSearchSideFlag)) {
                        aux_set.insert(p_neigh);
                        prev_layer_nodes.push_back(p_neigh);
                    }
                }
            }
        }
        // Check number of first layer points
        if (aux_set.size() == 0) {
            // KRATOS_WARNING("ShiftedBoundaryPointBasedUtility")
            //     << "No nodal neighbors on the other side were found for node " << pOtherSideNode->Id() << ". Extension basis can not be calculated." << std::endl;
            return;
        }

        // Add more layers of nodal neighbors of the current nodes to the cloud of nodes
        // Add those layers in normal direction of the boundary, so that both sides are more clearly separated at edges (3D) and tips of the skin geometry
        // NOTE that we start from 1 here as the first layer already has been added, so only one more layer will be added for a linear MLS extension
        for (std::size_t i_layer = 1; i_layer < n_layers; ++i_layer) {
            AddLateralSupportLayer(rAvgSkinPosition, rAvgSkinNormal, prev_layer_nodes, cur_layer_nodes, aux_set);
            prev_layer_nodes = cur_layer_nodes;
            cur_layer_nodes.clear();
        }

        // If there are not enough active support nodes to perform the MLS calculation add another layer of neighboring nodes
        // Add maximal three extra layers, these do not have to be in normal direction away from the averaged skin geometry anymore
        std::size_t n_cloud_nodes = aux_set.size();
        std::size_t n_extra_layers = 0;
        while (n_cloud_nodes <= GetRequiredNumberOfPoints() && n_extra_layers < 3) {
            AddLateralSupportLayer(prev_layer_nodes, cur_layer_nodes, aux_set);
            n_extra_layers++;
            n_cloud_nodes = aux_set.size();
            prev_layer_nodes = cur_layer_nodes;
            cur_layer_nodes.clear();
        }
        // if (n_extra_layers > 0) {
        //    KRATOS_WARNING("ShiftedBoundaryPointBasedUtility") << n_extra_layers << " extra layers of points needed for MLS calculation." << std::endl;
        // }

        // Add obtained cloud nodes to the cloud node vector and sort them by id
        //TODO sorting really necessary or helpful??
        rCloudNodes.resize(n_cloud_nodes);
        std::size_t aux_i = 0;
        for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
            rCloudNodes(aux_i++) = *it_set;
            // Mark support node for visualization - TODO make threadsafe/ put somewhere else for parallelization
            //(*it_set)->Set(rSearchSideFlag, true);
        }
        std::sort(rCloudNodes.ptr_begin(), rCloudNodes.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});

        // Fill the coordinates matrix
        rCloudCoordinates.resize(n_cloud_nodes, 3);
        IndexPartition<std::size_t>(n_cloud_nodes).for_each(array_1d<double,3>(), [&rCloudNodes, &rCloudCoordinates](std::size_t iNode, array_1d<double,3>& rAuxCoordTLS){
            noalias(rAuxCoordTLS) = rCloudNodes[iNode].Coordinates();
            rCloudCoordinates(iNode, 0) = rAuxCoordTLS[0];
            rCloudCoordinates(iNode, 1) = rAuxCoordTLS[1];
            rCloudCoordinates(iNode, 2) = rAuxCoordTLS[2];
        });
    }

    void ShiftedBoundaryPointBasedUtility::AddLateralSupportLayer(
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet)
    {
        //TODO SPEED UP

        // Find elemental neighbors of the nodes of the previous layer and add their nodes
        // NOTE that taking the nodes of neighboring elements is the same as adding the nodal neighbors directly for triangles and tetrahedra
        for (auto& p_node : PreviousLayerNodes) {
            const auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);

            // Add all nodes of neighboring elements to cloud nodes set if element is not (!) SBM_BOUNDARY
            // This way the boundary cannot be crossed (note that 'ACTIVE' might be used instead here)
            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                const auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                if (p_elem_neigh != nullptr) {
                    if (!p_elem_neigh->Is(SBM_BOUNDARY)) {
                        const auto& r_geom = p_elem_neigh->GetGeometry();
                        for (std::size_t i_neigh_node = 0; i_neigh_node < r_geom.PointsNumber(); ++i_neigh_node) {
                            // Add node of neighboring element to the support node set if it has not been added yet and is active
                            NodeType::Pointer p_neigh = r_geom(i_neigh_node);
                            if (p_neigh->Is(ACTIVE)) {
                                auto set_return = SupportNodesSet.insert(p_neigh);
                                // If the node was inserted into the set as a new element, then add it to the current layer (otherwise already visited nodes are visited again)
                                if (set_return.second) {
                                    CurrentLayerNodes.push_back(p_neigh);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedUtility::AddLateralSupportLayer(
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet)
    {
        //TODO SPEED UP

        // Find elemental neighbors of the nodes of the previous layer
        // NOTE that taking the nodes of neighboring elements is the same as adding the nodal neighbors directly for triangles and tetrahedra
        for (auto& p_node : PreviousLayerNodes) {
            const auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);

            // Add all nodes of neighboring elements to cloud nodes set if element is not (!) SBM_BOUNDARY and in normal direction
            // This way the boundary cannot be crossed (not 'ACTIVE' might be used instead here)
            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                const auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                if (p_elem_neigh != nullptr) {
                    if (!p_elem_neigh->Is(SBM_BOUNDARY)) {
                        const auto& r_geom = p_elem_neigh->GetGeometry();
                        for (std::size_t i_neigh_node = 0; i_neigh_node < r_geom.PointsNumber(); ++i_neigh_node) {
                            // Add node of neighboring element only if it is in the inwards normal direction of the element's averaged skin points
                            // NOTE this is done for a more robust separation of both sides of the boundary
                            NodeType::Pointer p_neigh = r_geom(i_neigh_node);
                            // Calculate dot product of average skin normal of the element and the normalized vector between averaged skin point and the element's node
                            array_1d<double,3> avg_skin_pt_to_node = p_neigh->Coordinates() - rAvgSkinPosition;
                            //avg_skin_pt_to_node /= norm_2(avg_skin_pt_to_node);  // normalization recommended for dot_product check another value than zero
                            const double dot_product = inner_prod(avg_skin_pt_to_node, rAvgSkinNormal);
                            if (p_neigh->Is(ACTIVE) && dot_product > 0.0) {
                                auto set_return = SupportNodesSet.insert(p_neigh);
                                // If the node was inserted into the set as a new element, then add it to the current layer (otherwise already visited nodes are visited again)
                                if (set_return.second) {
                                    CurrentLayerNodes.push_back(p_neigh);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedUtility::CreateCloudVectorsForSkinPointElement(
        const ElementType& rElement,
        const Vector& rSidesVector,
        NodesCloudMapType& rExtensionOperatorMap,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide)
    {
        //TODO SPEED UP

        // Create an auxiliary set with all the cloud nodes that affect the current element for each side separately
        // NOTE that a node can only be found if sufficient cloud nodes were found for the creation of the extension basis
        // NOTE that only active nodes are part of the extension operator support nodes
        NodesSetType cloud_nodes_set_pos;
        NodesSetType cloud_nodes_set_neg;
        const auto& r_geom = rElement.GetGeometry();
        for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            NodeType::Pointer p_node = r_geom(i_node);
            //TODO edge node treatment
            //if (p_node->Id() == 908 || p_node->Id() == 924 || p_node->Id() == 347 || p_node->Id() == 372) {  //p_node->Id() == 361 || p_node->Id() == 878 ||
            //    cloud_nodes_set_pos.insert(p_node);
            //    cloud_nodes_set_neg.insert(p_node);
            //}
            if (rSidesVector(i_node) > 0.0) {
                // Add positive side node to cloud nodes set of positive side of the boundary
                // NOTE they might not be part of the negative node's support because they are too close to the other side or not active
                cloud_nodes_set_pos.insert(p_node);
                // Add positive side's node's cloud nodes to cloud nodes set of negative side of the boundary
                const std::size_t found = rExtensionOperatorMap.count(p_node);
                if (found) {
                    auto& r_ext_op_data = rExtensionOperatorMap[p_node];
                    for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                        auto& p_cl_node = std::get<0>(*it_data);
                        cloud_nodes_set_neg.insert(p_cl_node);
                    }
                }
            } else {
                // Add negative side node to cloud nodes set of negative side of the boundary
                // NOTE they might not be part of the positive node's support because they are too close to the other side or not active
                cloud_nodes_set_neg.insert(p_node);
                // Add negative side's node's cloud nodes to cloud nodes set of positive side of the boundary
                const std::size_t found = rExtensionOperatorMap.count(p_node);
                if (found) {
                    auto& r_ext_op_data = rExtensionOperatorMap[p_node];
                    for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                        auto& p_cl_node = std::get<0>(*it_data);
                        cloud_nodes_set_pos.insert(p_cl_node);
                    }
                }
            }
        }

        // Save node clouds in pointer vectors to be used in the creation of the condition
        const std::size_t n_cloud_nodes_pos = cloud_nodes_set_pos.size();
        const std::size_t n_cloud_nodes_neg = cloud_nodes_set_neg.size();
        rCloudNodeVectorPositiveSide.resize(n_cloud_nodes_pos);
        rCloudNodeVectorNegativeSide.resize(n_cloud_nodes_neg);
        std::size_t aux_i = 0;
        for (auto it_set = cloud_nodes_set_pos.begin(); it_set != cloud_nodes_set_pos.end(); ++it_set) {
            rCloudNodeVectorPositiveSide(aux_i++) = *it_set;
        }
        aux_i = 0;
        for (auto it_set = cloud_nodes_set_neg.begin(); it_set != cloud_nodes_set_neg.end(); ++it_set) {
            rCloudNodeVectorNegativeSide(aux_i++) = *it_set;
        }

        // Sort obtained cloud node vectors by ID to properly get the extension operator data  //TODO really necessary or faster??
        std::sort(rCloudNodeVectorPositiveSide.ptr_begin(), rCloudNodeVectorPositiveSide.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});
        std::sort(rCloudNodeVectorNegativeSide.ptr_begin(), rCloudNodeVectorNegativeSide.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});
    }

    void ShiftedBoundaryPointBasedUtility::GetDataForSkinPointInElement(
        const ElementType& rElement,
        const array_1d<double,3>& rSkinPtCoordinates,
        Vector& rSkinPtShapeFunctionValues,
        Matrix& rSkinPtShapeFunctionDerivatives)
    {
        //TODO get this when searching for skin point for the first time and store it in map

        const auto& r_geom = rElement.GetGeometry();

        // Compute the local coordinates of the integration point in the element's geometry
        array_1d<double,3> int_pt_local_coords = ZeroVector(3);
        r_geom.PointLocalCoordinates(int_pt_local_coords, rSkinPtCoordinates);

        // Get N of the element at the integration point
        r_geom.ShapeFunctionsValues(rSkinPtShapeFunctionValues, int_pt_local_coords);

        // Get DN_DX of the element at the integration point
        Matrix aux_DN_DXi_parent, aux_J_parent, aux_J_inv_parent;
        double aux_detJ_parent;
        r_geom.ShapeFunctionsLocalGradients(aux_DN_DXi_parent, int_pt_local_coords);
        r_geom.Jacobian(aux_J_parent, int_pt_local_coords);
        MathUtils<double>::InvertMatrix(aux_J_parent, aux_J_inv_parent, aux_detJ_parent);
        rSkinPtShapeFunctionDerivatives = prod(aux_DN_DXi_parent, aux_J_inv_parent);
    }

    bool ShiftedBoundaryPointBasedUtility::AddSkinPointCondition(
        const ElementType& rElement,
        const Vector& rSidesVector,
        const double ElementSize,
        const array_1d<double,3>& rSkinPtCoordinates,
        const array_1d<double,3>& rSkinPtAreaNormal,
        NodesCloudMapType& rExtensionOperatorMap,
        const PointerVector<NodeType>& rCloudNodeVector,
        const Vector& rSkinPtShapeFunctionValues,
        const Matrix& rSkinPtShapeFunctionDerivatives,
        std::size_t& r_ConditionId,
        const bool ConsiderPositiveSide)
    {
        //TODO SPEED UP

        const auto& r_geom = rElement.GetGeometry();

        // Initialize the extension operator containers
        const std::size_t n_cl_nodes = rCloudNodeVector.size();
        const std::size_t n_dim = r_geom.WorkingSpaceDimension();
        Vector N_container = ZeroVector(n_cl_nodes);
        Matrix DN_DX_container = ZeroMatrix(n_cl_nodes, n_dim);

        array_1d<double,3> area_normal = rSkinPtAreaNormal;
        double skin_pt_weight = norm_2(area_normal);

        // Do not add wall condition for 0D skin elements, this also prevents problems with zero normal
        if (skin_pt_weight < 1e-10) {
            return false;
        }

        // Loop the nodes that are involved in the current element
        for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            const auto p_node = r_geom(i_node);
            // If node is on the side that is being considered, then add the standard shape function contribution of the node at the position of the skin point
            if (ConsiderPositiveSide != (rSidesVector[i_node] <= 0.0)) {
                // If a node on the side that is being considered is not active, then no wall condition is created
                if (!p_node->Is(ACTIVE)) {
                    // if (ConsiderPositiveSide) {
                    //     KRATOS_WARNING("ShiftedBoundaryPointBasedUtility") << "No wall condition will be created for positive side of the skin point because Node No." << p_node->Id() << " is not active." << std::endl;
                    // } else {
                    //     KRATOS_WARNING("ShiftedBoundaryPointBasedUtility") << "No wall condition will be created for negative side of the skin point because Node No." << p_node->Id() << " is not active." << std::endl;
                    // }
                    return false;
                }
                // Note that we need to check for the ids to match as we do not know the node's position in the node vector
                for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                    auto& p_cl_node = rCloudNodeVector(i_cl);
                    if (p_node->Id() == p_cl_node->Id()) {
                        N_container(i_cl) += rSkinPtShapeFunctionValues(i_node);
                        for (std::size_t d = 0; d < n_dim; ++d) {
                            DN_DX_container(i_cl, d) += rSkinPtShapeFunctionDerivatives(i_node, d);
                        }
                        break;
                    }
                }
            // If node is on the other side of the boundary, then get its shape function values and derivatives for the skin point and its extension operator data
            } else {
                // Get the weight as the corresponding nodal shape function value of the node at the position of the skin point
                const double i_node_N = rSkinPtShapeFunctionValues(i_node);
                const auto i_node_grad_N = row(rSkinPtShapeFunctionDerivatives, i_node);

                // If node on the other side does not have an extension basis, then no wall condition is created
                const std::size_t found = rExtensionOperatorMap.count(p_node);
                if (!found) {
                    //KRATOS_WARNING("ShiftedBoundaryPointBasedUtility") << "No wall condition will be created for one side of the skin point because no extension operator was available for Node No." << p_node->Id() << std::endl;
                    return false;
                }

                // Get the node's extension operator data
                const auto& ext_op_data = rExtensionOperatorMap[p_node];

                // Iterate over the node's extension operator data and apply the support node weight (i_cl_node_N) to make the basis conformant
                // Note that we need to check for the ids to match as we do not know the node's position in the node vector
                for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
                    const auto& r_node_data = *it_data;
                    const std::size_t data_node_id = (std::get<0>(r_node_data))->Id();
                    for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                        const auto& p_cl_node = rCloudNodeVector(i_cl);
                        if (p_cl_node->Id() == data_node_id) {
                            const double i_cl_node_N = std::get<1>(r_node_data);
                            N_container(i_cl) += i_node_N * i_cl_node_N;
                            for (std::size_t d = 0; d < n_dim; ++d) {
                                DN_DX_container(i_cl,d) += i_node_grad_N(d) * i_cl_node_N;
                            }
                            break;
                        }
                    }
                }
            }
        }

        /*const std::size_t n_nodes = r_geom.PointsNumber();
        PointerVector<NodeType> cloud_node_vector;
        cloud_node_vector.resize(n_nodes);
        Vector N_container_2 = ZeroVector(n_nodes);
        Matrix DN_DX_container_2 = ZeroMatrix(n_nodes, n_dim);
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            cloud_node_vector(i_node) = r_geom(i_node);
            N_container_2(i_node) = rSkinPtShapeFunctionValues(i_node);
            for (std::size_t d = 0; d < n_dim; ++d) {
                DN_DX_container_2(i_node, d) += rSkinPtShapeFunctionDerivatives(i_node, d);
            }
        }*/

        // Create a new condition with a geometry made up with the basis nodes
        auto p_prop = rElement.pGetProperties();
        auto p_cond = mpConditionPrototype->Create(++r_ConditionId, rCloudNodeVector, p_prop);
        p_cond->Set(ACTIVE, true);
        mpBoundarySubModelPart->AddCondition(p_cond);

        //TODO for laplacian static heat testing:
        // Set Dirichlet boundary condition
        //const double dirichlet_value = std::pow(rSkinPtCoordinates(0),2) + std::pow(rSkinPtCoordinates(1),2);
        //p_cond->SetValue(TEMPERATURE, dirichlet_value);

        // Store the SBM BC data in the condition database
        p_cond->SetValue(ELEMENT_H, ElementSize);
        p_cond->SetValue(INTEGRATION_COORDINATES, rSkinPtCoordinates);
        p_cond->SetValue(NORMAL, area_normal/ skin_pt_weight);
        p_cond->SetValue(INTEGRATION_WEIGHT, skin_pt_weight);
        p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
        p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);

        return true;
    }

    bool ShiftedBoundaryPointBasedUtility::FixPressureOfEnclosedNode(
        ElementType& rElement,
        const Vector& rSidesVector,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal)
    {
        auto& r_geom = rElement.GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {

            if ((mPositiveSideIsEnclosed && rSidesVector[i_node] > 0) or (mNegativeSideIsEnclosed && rSidesVector[i_node] < 0)) {
                auto& r_node = r_geom[i_node];
                if (r_node.Is(ACTIVE)) {
                    r_node.Fix(PRESSURE);
                    r_node.FastGetSolutionStepValue(PRESSURE) = 0.0;
                    return true;
                }
            }
        }
        return false;
    }

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPointsTemplated()
    {
        constexpr std::size_t voigt_size = 3 * (TDim-1);
        constexpr std::size_t block_size = TDim +1;

        std::size_t n_elements_without_correct_extension = 0;

        // Loop over all elements containing skin points to integrate traction=sigma*n over the interface
        // NOTE that both interface sides need to be integrated
        LockObject mutex_valid_ex;
        std::for_each(mSkinPointsMap.begin(), mSkinPointsMap.end(), [&](const std::pair<ElementType::Pointer, SkinPointsDataVectorType>& rKeyData){
            const auto p_element = rKeyData.first;
            const auto& skin_points_data_vector = rKeyData.second;
            const auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            const std::size_t local_size = n_nodes * block_size;

            // Calculate unknowns at the nodes of the element for the positive and the negative side of Gamma
            Vector unknowns_pos = ZeroVector(local_size);
            Vector unknowns_neg = ZeroVector(local_size);
            const bool unknowns_calculated_successfully = CalculateUnknownsForSkinPointElement<TDim>(p_element, unknowns_pos, unknowns_neg);
            if (!unknowns_calculated_successfully) {
                std::scoped_lock<LockObject> lock(mutex_valid_ex);
                n_elements_without_correct_extension++;
            }

            // Get constitutive law, see FluidElementData and FluidElement::CalculateMaterialResponse
            const auto p_constitutive_law =  p_element->GetProperties()[CONSTITUTIVE_LAW];
            ConstitutiveLaw::Parameters constitutive_law_values(r_geom, p_element->GetProperties(), mpModelPart->GetProcessInfo());
            Flags& r_cl_options = constitutive_law_values.GetOptions();
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            // Iterate over the element's skin points adding a positive side and a negative side drag contribution
            for (std::size_t i_skin_pt = 0; i_skin_pt < skin_points_data_vector.size(); ++i_skin_pt) {
                // Get the skin point's position and area normal
                const auto skin_pt_data = skin_points_data_vector[i_skin_pt];
                const array_1d<double,3> skin_pt_position = std::get<0>(skin_pt_data);
                const array_1d<double,3> skin_pt_area_normal = std::get<1>(skin_pt_data);
                const std::size_t skin_pt_node_id = std::get<2>(skin_pt_data);

                // Get normal and represented area of skin point
                // NOTE that the normal is defined to point in fluid inward direction for the positive side of the skin and fluid outward for the negative side of the skin
                const double skin_pt_area = norm_2(skin_pt_area_normal);
                const array_1d<double,3> aux_unit_normal = skin_pt_area_normal / skin_pt_area;

                // Get shape function values and derivatives of the element at the skin point
                Vector skin_pt_N(n_nodes);
                Matrix skin_pt_DN_DX = ZeroMatrix(n_nodes, TDim);
                GetDataForSkinPointInElement(*p_element, skin_pt_position, skin_pt_N, skin_pt_DN_DX);

                // Calculate velocity and pressure at skin point for positive and negative side of Gamma
                array_1d<double, 3> u_pos = ZeroVector(3);
                array_1d<double, 3> u_neg = ZeroVector(3);
                double p_pos = 0.0;
                double p_neg = 0.0;
                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                    const std::size_t row = i_node*block_size;
                    for (std::size_t d = 0; d < TDim; ++d) {
                        u_pos[d] += skin_pt_N(i_node) * unknowns_pos(row+d);
                        u_neg[d] += skin_pt_N(i_node) * unknowns_neg(row+d);
                    }
                    p_pos += skin_pt_N(i_node) * unknowns_pos(row+TDim);
                    p_neg += skin_pt_N(i_node) * unknowns_neg(row+TDim);
                }
                // NOTE that the pressure on the positive side of Gamma needs to be multiplied by the negative normal of the skin point,
                // so it is pointing outwards of the fluid. The pressure on the negative side points outwards multiplied by the positive normal.
                const array_1d<double, 3> traction_p = (p_neg-p_pos) * aux_unit_normal;

                // Get the normal projection matrix in Voigt notation
                BoundedMatrix<double, TDim, voigt_size> voigt_normal_proj_matrix = ZeroMatrix(TDim, voigt_size);
                ShiftedBoundaryUtilityInternals::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

                // Calculate strain rate for positive and negative side of Gamma
                Matrix B_matrix = ZeroMatrix(voigt_size, local_size);
                ShiftedBoundaryUtilityInternals::CalculateStrainMatrix<TDim>(skin_pt_DN_DX, n_nodes, B_matrix);
                Vector strain_rate_pos = prod(B_matrix, unknowns_pos);
                Vector strain_rate_neg = prod(B_matrix, unknowns_neg);

                // Calculate the stress for the positive side
                Vector shear_stress_pos = ZeroVector(voigt_size);
                constitutive_law_values.SetStrainVector(strain_rate_pos);           //input
                constitutive_law_values.SetStressVector(shear_stress_pos);          //output
                p_constitutive_law->CalculateMaterialResponseCauchy(constitutive_law_values);

                // Calculate the stress for the negative side
                Vector shear_stress_neg = ZeroVector(voigt_size);
                constitutive_law_values.SetStrainVector(strain_rate_neg);           //input
                constitutive_law_values.SetStressVector(shear_stress_neg);          //output
                p_constitutive_law->CalculateMaterialResponseCauchy(constitutive_law_values);

                // Calculate shear stress at the skin point
                // NOTE that the negative side shear stress needs to be multiplied by the negative normal, so that the normal is pointing inwards.
                const array_1d<double, TDim> shear_proj_pos =  prod(voigt_normal_proj_matrix, shear_stress_pos);
                const array_1d<double, TDim> shear_proj_neg = -prod(voigt_normal_proj_matrix, shear_stress_neg);
                array_1d<double, 3> traction_tau = ZeroVector(3);
                for (std::size_t d = 0; d < TDim ; ++d) {
                    traction_tau(d) += shear_proj_pos(d) + shear_proj_neg(d);
                }
                const array_1d<double, 3> skin_pt_force = skin_pt_area * (traction_p + traction_tau);

                // Store velocity, pressure and traction in a skin point model part
                auto& skin_pt_in_model_part = mpSkinPointsSubModelPart->GetNode(skin_pt_node_id);
                skin_pt_in_model_part.FastGetSolutionStepValue(POSITIVE_FACE_FLUID_VELOCITY) = u_pos;
                skin_pt_in_model_part.FastGetSolutionStepValue(NEGATIVE_FACE_FLUID_VELOCITY) = u_neg;
                skin_pt_in_model_part.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = p_pos;
                skin_pt_in_model_part.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) = p_neg;
                skin_pt_in_model_part.FastGetSolutionStepValue(TRACTION_FROM_FLUID_PRESSURE) = traction_p;
                skin_pt_in_model_part.FastGetSolutionStepValue(TRACTION_FROM_FLUID_STRESS) = traction_tau;
                skin_pt_in_model_part.FastGetSolutionStepValue(DRAG_FORCE) = skin_pt_force;
            }
        });
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedUtility", n_elements_without_correct_extension > 0)
        << "The unknowns inside " << n_elements_without_correct_extension << " elements which contain skin integration points were calculated without valid extension." << std::endl;
    }

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPointsAndNodesTemplated()
    {
        // Calculate variables at all skin points (integration points) and store them in the skin points model part
        CalculateVariablesAtSkinPointsTemplated<TDim>();

        // Create dictionaries to accumulate the weighted values and weights for each skin node for interpolation
        std::unordered_map<std::size_t, array_1d<double,2>> dict_sum_weighted_p; //positive and negative side pressure
        std::unordered_map<std::size_t, array_1d<array_1d<double,3>,2>> dict_sum_weighted_u; //positive and negative side velocity
        std::unordered_map<std::size_t, double> dict_sum_weights;
        for (auto& r_skin_node : mpSkinDiscSubModelPart->Nodes()) {
            dict_sum_weighted_p[r_skin_node.Id()] = ZeroVector(2);
            dict_sum_weighted_u[r_skin_node.Id()][0] = ZeroVector(3);
            dict_sum_weighted_u[r_skin_node.Id()][1] = ZeroVector(3);
            dict_sum_weights[r_skin_node.Id()] = 0.0;
        }

        //TODO faster to store the skin points element and shape function values once?!?
        // Loop over all skin points and locate them in the skin model part to add their contribution to the skin nodes
        //TODO use local vectors and index instead of mutex!
        LockObject mutex;
        block_for_each(mpSkinPointsSubModelPart->Nodes(), [&](NodeType& rSkinPoint){
            // Search for the skin point in the skin mesh to get the element containing the point and the shape function values
            array_1d<double,3> skin_pt_position = rSkinPoint.Coordinates();
            Vector skin_pt_N(TDim+1); //Note that this restricts the use to tri and tetra
            Element::Pointer p_element = nullptr;

            // Create bounding box for the skin point
            ShiftedBoundaryUtilityInternals::BBox pt_box;
            pt_box.min[0] = pt_box.max[0] = skin_pt_position[0];
            pt_box.min[1] = pt_box.max[1] = skin_pt_position[1];
            pt_box.min[2] = pt_box.max[2] = skin_pt_position[2];

            // Get skin element candidates from the skin AABB tree
            std::vector<std::size_t> candidates;
            candidates.reserve(20);
            mSkinBVH.query(pt_box, candidates);  // read-only query of the BVH - safe for parallel execution

            // Check candidates until the skin element containing the skin point is found and get its shape function values
            Element::GeometryType::CoordinatesArrayType local_coords;
            bool is_inside = false;
            for (auto elem_idx : candidates) {
                auto p_elem = mBVHIdxToSkinElementVector[elem_idx];
                auto& r_geom = p_elem->GetGeometry();
                is_inside = r_geom.IsInside(skin_pt_position, local_coords);
                if(is_inside) {
                    p_element = p_elem;
                    r_geom.ShapeFunctionsValues(skin_pt_N, local_coords);
                    break;
                }
            }

            if (is_inside) {
                // Get skin point variables
                const double& r_p_pos = rSkinPoint.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                const double& r_p_neg = rSkinPoint.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
                const array_1d<double,3>& r_u_pos = rSkinPoint.FastGetSolutionStepValue(POSITIVE_FACE_FLUID_VELOCITY);
                const array_1d<double,3>& r_u_neg = rSkinPoint.FastGetSolutionStepValue(NEGATIVE_FACE_FLUID_VELOCITY);

                // Loop over the nodes of the element to add the skin point's contribution to each node
                const auto& r_geom = p_element->GetGeometry();
                for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                    const std::size_t node_id = r_geom[i_node].Id();
                    const double N_i = skin_pt_N(i_node);

                    std::scoped_lock<LockObject> lock(mutex);
                    dict_sum_weighted_p[node_id][0] += N_i * r_p_pos;
                    dict_sum_weighted_p[node_id][1] += N_i * r_p_neg;
                    dict_sum_weighted_u[node_id][0] += N_i * r_u_pos;
                    dict_sum_weighted_u[node_id][1] += N_i * r_u_neg;
                    dict_sum_weights[node_id] += N_i;
                }
            } else {
                KRATOS_WARNING("ShiftedBoundaryPointBasedUtility") << "Skin point " << rSkinPoint.Id() << " could not be located in the skin mesh." << std::endl;
            }
        });

        // Set the skin nodes variables based on the skin points contributions
        block_for_each(mpSkinDiscSubModelPart->Nodes(), [&](NodeType& rSkinNode){
            const auto& sum_weighted_p = dict_sum_weighted_p[rSkinNode.Id()];
            const auto& sum_weighted_u = dict_sum_weighted_u[rSkinNode.Id()];
            const double& sum_weights = dict_sum_weights[rSkinNode.Id()];

            rSkinNode.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE)       = sum_weighted_p[0] / sum_weights;
            rSkinNode.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE)       = sum_weighted_p[1] / sum_weights;
            rSkinNode.FastGetSolutionStepValue(POSITIVE_FACE_FLUID_VELOCITY) = sum_weighted_u[0] / sum_weights;
            rSkinNode.FastGetSolutionStepValue(NEGATIVE_FACE_FLUID_VELOCITY) = sum_weighted_u[1] / sum_weights;
        });
    }

    template <std::size_t TDim>
    bool ShiftedBoundaryPointBasedUtility::CalculateUnknownsForSkinPointElement(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns)
    {
        constexpr std::size_t block_size = TDim +1;
        const auto& r_geom = pElement->GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();

        Vector sides_vector(n_nodes);
        const std::size_t sides_found = mSidesVectorMap.count(pElement);
        // Typically this method is called for elements containing skin points, so a side vector should be available.
        if (sides_found) {
            // Get sides vector for the element the skin node is located in
            sides_vector = mSidesVectorMap[pElement];
        } else {
            return false;
        }

        // Calculate positive and negative side unknowns at all nodes of the element
        bool element_is_without_extension = false;
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            const auto p_node = r_geom(i_node);
            const std::size_t row = i_node * block_size;

            // Check whether extension exists for node
            const std::size_t extension_found = mExtensionOperatorMap.count(p_node);

            // Initialize positive and negative side velocity and pressure
            array_1d<double,3> u_node_pos = ZeroVector(3);
            array_1d<double,3> u_node_neg = ZeroVector(3);
            double p_node_pos = 0.0;
            double p_node_neg = 0.0;

            // Get unknowns at the node for positive and negative side
            // The extension needs to be used for values on the other side on which the node is located
            if (sides_vector[i_node] > 0.0) {
                // Get nodal velocity and pressure directly for the positive side
                u_node_pos = p_node->FastGetSolutionStepValue(VELOCITY);
                p_node_pos = p_node->FastGetSolutionStepValue(PRESSURE);
                // Calculate nodal velocity and pressure using the extension operator of the node for the negative side
                if (extension_found) {
                    const auto& ext_op_data = mExtensionOperatorMap[p_node];
                    // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the value at the support node
                    for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
                        const auto p_support_node = std::get<0>(*it_data);
                        const double weight_support_node = std::get<1>(*it_data);
                        u_node_neg += weight_support_node * p_support_node->FastGetSolutionStepValue(VELOCITY);
                        p_node_neg += weight_support_node * p_support_node->FastGetSolutionStepValue(PRESSURE);
                    }
                } else { element_is_without_extension = true; }
            } else {
                // Calculate nodal velocity and pressure using the extension operator of the node for the positive side
                if (extension_found) {
                    const auto& ext_op_data = mExtensionOperatorMap[p_node];
                    // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the value at the support node
                    for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
                        const auto p_support_node = std::get<0>(*it_data);
                        const double weight_support_node = std::get<1>(*it_data);
                        u_node_pos += weight_support_node * p_support_node->FastGetSolutionStepValue(VELOCITY);
                        p_node_pos += weight_support_node * p_support_node->FastGetSolutionStepValue(PRESSURE);
                    }
                } else { element_is_without_extension = true; }
                // Get nodal velocity and pressure directly for the negative side
                u_node_neg = p_node->FastGetSolutionStepValue(VELOCITY);
                p_node_neg = p_node->FastGetSolutionStepValue(PRESSURE);
            }

            // Store positive and negative side unknowns at the node
            for (std::size_t d = 0; d < TDim; ++d) {
                rPositiveSideUnknowns(row+d) = u_node_pos(d);
                rNegativeSideUnknowns(row+d) = u_node_neg(d);
            }
            rPositiveSideUnknowns(row+TDim) = p_node_pos;
            rNegativeSideUnknowns(row+TDim) = p_node_neg;
        }

        if (element_is_without_extension) { return false; }
        else { return true; }
    }

    ShiftedBoundaryPointBasedUtility::MeshlessShapeFunctionsFunctionType ShiftedBoundaryPointBasedUtility::GetMLSShapeFunctionsFunction() const
    {
        switch (mpModelPart->GetProcessInfo()[DOMAIN_SIZE]) {
            case 2:
                switch (mMLSExtensionOperatorOrder) {
                    case 1:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                            MLSShapeFunctionsUtility::CalculateShapeFunctions<2,1>(rPoints, rX, h, rN);};
                    case 2:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                            MLSShapeFunctionsUtility::CalculateShapeFunctions<2,2>(rPoints, rX, h, rN);};
                    default:
                        KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
                }
            case 3:
                switch (mMLSExtensionOperatorOrder) {
                    case 1:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                            MLSShapeFunctionsUtility::CalculateShapeFunctions<3,1>(rPoints, rX, h, rN);};
                    case 2:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                            MLSShapeFunctionsUtility::CalculateShapeFunctions<3,2>(rPoints, rX, h, rN);};
                    default:
                        KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
                }
            default:
                KRATOS_ERROR << "Wrong domain size. MLS shape functions utility cannot be set.";
        }
    }

    ShiftedBoundaryPointBasedUtility::MeshlessShapeFunctionsFunctionType ShiftedBoundaryPointBasedUtility::GetRBFShapeFunctionsFunction() const
    {
        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
            RBFShapeFunctionsUtility::CalculateShapeFunctions(rPoints, rX, h, rN);
        };
    }

    ShiftedBoundaryPointBasedUtility::ElementSizeFunctionType ShiftedBoundaryPointBasedUtility::GetElementSizeFunction(const GeometryType& rGeometry)
    {
        switch (rGeometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return [](const GeometryType& rGeometry)->double{return ElementSizeCalculator<2,3>::AverageElementSize(rGeometry);};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return [](const GeometryType& rGeometry)->double{return ElementSizeCalculator<3,4>::AverageElementSize(rGeometry);};
            default:
                KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
        }
    }

    double ShiftedBoundaryPointBasedUtility::CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin)
    {
        const std::size_t n_nodes = rCloudCoordinates.size1();
        const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
            return std::pow(rCloudCoordinates(I,0) - rOrigin(0),2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1),2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2),2);
        });
        return std::sqrt(squared_rad);
    }

    std::size_t ShiftedBoundaryPointBasedUtility::GetRequiredNumberOfPoints()
    {
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        switch (n_dim) {
            case 2:
                switch (mMLSExtensionOperatorOrder) {
                    case 1:
                        return 3;
                    case 2:
                        return 6;
                    default:
                        KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
                }
            case 3:
                switch (mMLSExtensionOperatorOrder) {
                    case 1:
                        return 4;
                    case 2:
                        return 10;
                    default:
                        KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
                }
            default:
                KRATOS_ERROR << "Wrong domain size.";
        }
    }

    //----------------------------------------------------------------
    //     TEMPLATE INSTANTIATIONS
    //----------------------------------------------------------------

    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPointsTemplated<2>();
    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPointsTemplated<3>();

    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPointsAndNodesTemplated<2>();
    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::CalculateVariablesAtSkinPointsAndNodesTemplated<3>();

    template bool KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::CalculateUnknownsForSkinPointElement<2>(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns);
    template bool KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::CalculateUnknownsForSkinPointElement<3>(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns);

    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::MapSkinPointsToElementsTemplated<2>(SkinPointsToElementsMapType& rSkinPointsMap);
    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedUtility::MapSkinPointsToElementsTemplated<3>(SkinPointsToElementsMapType& rSkinPointsMap);

}  // namespace Kratos.
