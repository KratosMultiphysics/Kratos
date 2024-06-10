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
#include "containers/pointer_vector.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/element_size_calculator.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/shifted_boundary_meshless_discontinuous_interface_utility.h"
#include "utilities/shifted_boundary_meshless_interface_utility.h"
#include <cstddef>
#include <ostream>
#include <vector>

namespace Kratos
{

namespace
{
    using GeometryType = ShiftedBoundaryMeshlessInterfaceUtility::GeometryType;
    using ModifiedShapeFunctionsFactoryType = ShiftedBoundaryMeshlessInterfaceUtility::ModifiedShapeFunctionsFactoryType;

    /**
     * @brief Check if current geometry is split
     * This method checks if current geometry is split from the nodal historical level set values
     * @param rElement Element to be checked
     * @param rLevelSetVariable Variable storing the level set function
     * @return true If split
     * @return false If not split
     */
    bool IsSplit(
        const ElementType& rElement,
        const Variable<Vector>& rLevelSetVariable)
    {
        std::size_t n_neg = 0;
        std::size_t n_pos = 0;

        const Vector& r_distances = rElement.GetValue(rLevelSetVariable);

        for (std::size_t i_node = 0; i_node < r_distances.size(); ++i_node){
            if (r_distances[i_node] < 0.0){
                n_neg++;
            } else {
                n_pos++;
            }
        }
        return (n_pos != 0) && (n_neg != 0);
    }

    /**
     * @brief Checks if the current element is partially intersected (incised)
     * Checks if the current element is partially intersected by checking the number of extrapolated intersected edges
     * This number will only be non-zero if user provided flag to calculate extrapolated edge distances.
     * The case in which three edges of a tetrahedra are cut and element is only incised is also considered.
     * @param rElement Element to be checked
     * @param rLevelSetVariable Variable storing the level set function
     * @return true if the element is incised
     * @return false if the element is not incised
     */
    bool IsIncised(
        const ElementType& rElement,
        const Variable<Vector>& rLevelSetVariable)
    {
        std::size_t n_intersected_edges_extrapolated = 0;  //TODO

        const Vector& r_edge_distances_extrapolated = rElement.GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);

        // Number of edges cut by extrapolated geometry, if not empty
        for (std::size_t i = 0; i < r_edge_distances_extrapolated.size(); ++i) {
            if (r_edge_distances_extrapolated[i] > 0.0) {
                n_intersected_edges_extrapolated++;
            }
        }

        return n_intersected_edges_extrapolated > 0 ? true : false;
    }

    /**
     * @brief Check if node has negative level set value
     * This method checks if given nodes has a negative nodal historical level set value.
     * @param rNode Node to be checked
     * @param rLevelSetVariable Variable storing the level set function
     * @return true If negative level set value
     * @return false If level set value is positive or zero
     */
    bool IsNegative(
        const ElementType& rElement,
        const std::size_t NodeIndex,
        const Variable<Vector>& rLevelSetVariable)
    {
        const Vector& r_distances = rElement.GetValue(rLevelSetVariable);
        if (r_distances(NodeIndex) < 0.0) {
            return true;
        }
        return false;
    }

    /**
     * @brief Set the nodal distances vector
     * This method saves the elemental historical values of level set in the provided vector
     * @param rElement Element from which the distances are retrieved
     * @param rLevelSetVariable Variable storing the level set function
     * @param rNodalDistances Vector container to store the distance values
     */
    void SetNodalDistancesVector(
        const ElementType& rElement,
        const Variable<Vector>& rLevelSetVariable,
        Vector& rNodalDistances)
    {
        const auto& r_geom = rElement.GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();
        if (rNodalDistances.size() != n_nodes) {
            rNodalDistances.resize(n_nodes);
        }
        rNodalDistances = rElement.GetValue(rLevelSetVariable);
    }
    /**
     * @brief Get the standard modified shape functions factory object
     * This function returns a prototype for the split shape functions calculation from the provided geometry
     * @param rGeometry Input geometry
     * @return ModifiedShapeFunctionsFactoryType Factory to be used for the split shape functions calculation
     */
    ModifiedShapeFunctionsFactoryType GetStandardModifiedShapeFunctionsFactory(const GeometryType& rGeometry)
    {
        switch (rGeometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
            default:
                KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
        }
    }
}

    ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility(
        Model& rModel,
        Parameters ThisParameters) : ShiftedBoundaryMeshlessInterfaceUtility()
    {

        // Validate input settings with defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Retrieve the required model parts
        const std::string model_part_name = ThisParameters["model_part_name"].GetString();
        const std::string boundary_sub_model_part_name = ThisParameters["boundary_sub_model_part_name"].GetString();
        mpModelPart = &rModel.GetModelPart(model_part_name);
        if (mpModelPart->HasSubModelPart(boundary_sub_model_part_name)) {
            mpBoundarySubModelPart = &(mpModelPart->GetSubModelPart(boundary_sub_model_part_name));
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfNodes() != 0) << "Provided SBM model part has nodes." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfElements() != 0) << "Provided SBM model part has elements." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfConditions() != 0) << "Provided SBM model part has conditions." << std::endl;
        } else {
            mpBoundarySubModelPart = &(mpModelPart->CreateSubModelPart(boundary_sub_model_part_name));
        }

        // Save a pointer to variable storing the level set
        const std::string levelset_variable_name = ThisParameters["levelset_variable_name"].GetString();
        if (KratosComponents<Variable<Vector>>::Has(levelset_variable_name)) {
            mpDiscontinuousLevelSetVariable = &KratosComponents<Variable<Vector>>::Get(levelset_variable_name);
        } else {
            KRATOS_ERROR << "Provided 'levelset_variable_name' " << levelset_variable_name << " is not registered. Be aware that discontinuous level set needs to be used for discontinuous SBM." << std::endl;
        }

        // Set the order of the MLS extension operator used in the MLS shape functions utility
        mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();

        // If true, the basis is created such that it is conforming with the linear FE space of the surrogate boundary
        mConformingBasis = ThisParameters["conforming_basis"].GetBool();

        // If true, the basis is created such that the surrogate boundary gradient is kept
        const std::string ext_op_type = ThisParameters["extension_operator_type"].GetString();
        if (ext_op_type == "MLS") {
            mExtensionOperator = ShiftedBoundaryMeshlessInterfaceUtility::ExtensionOperator::MLS;
        } else if (ext_op_type == "RBF") {
            mExtensionOperator = ShiftedBoundaryMeshlessInterfaceUtility::ExtensionOperator::RBF;
            KRATOS_WARNING("ShiftedBoundaryMeshlessInterfaceUtility") << "Only 'MLS' 'extension_operator_type' is supported by discontinuous shifted boundary interface utility." << std::endl;
        } else if (ext_op_type == "gradient_based") {
            mExtensionOperator = ShiftedBoundaryMeshlessInterfaceUtility::ExtensionOperator::GradientBased;
            KRATOS_WARNING("ShiftedBoundaryMeshlessInterfaceUtility") << "Only 'MLS' 'extension_operator_type' is supported by discontinuous shifted boundary interface utility." << std::endl;
        } else {
            KRATOS_ERROR << "Wrong 'extension_operator_type' provided. Available options are 'MLS', 'RBF' and 'gradient_based'." << std::endl;
        }

        // Check that the basis settings are correct
        KRATOS_ERROR_IF(mExtensionOperator == ShiftedBoundaryMeshlessInterfaceUtility::ExtensionOperator::RBF && !mConformingBasis)
            << "'RBF' extension operator can only be used with conforming basis." << std::endl;
        KRATOS_ERROR_IF(mExtensionOperator == ShiftedBoundaryMeshlessInterfaceUtility::ExtensionOperator::GradientBased && !mConformingBasis)
            << "'gradient_based' extension operator can only be used with conforming basis." << std::endl;

        // Set the SBD condition prototype to be used in the condition creation
        std::string interface_condition_name = ThisParameters["sbm_interface_condition_name"].GetString();
        KRATOS_ERROR_IF(interface_condition_name == "") << "SBM interface condition has not been provided." << std::endl;
        mpConditionPrototype = &KratosComponents<Condition>::Get(interface_condition_name);
    }

    void ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::CalculateMeshlessBasedConformingExtensionBasis()
    {
        // Set the required interface flags
        SetInterfaceFlags();

        // Set the modified shape functions factory
        // Note that unique geometry in the mesh is assumed
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements())
            << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;
        const auto& r_begin_geom = mpModelPart->ElementsBegin()->GetGeometry();
        auto p_mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_begin_geom);

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        auto p_element_size_func = GetElementSizeFunction(r_begin_geom);

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Loop the elements to create the an extension basis for each node of all intersected elements(MLS shape functions values for support cloud of node)
        NodesCloudMapType ext_op_map;
        SetExtensionOperatorsForSplitElementNodes(ext_op_map);

        // Create the interface conditions
        //TODO: THIS CAN BE PARALLEL (WE JUST NEED TO MAKE CRITICAL THE CONDITION ID UPDATE)
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto& r_geom = rElement.GetGeometry();
            if (IsSplit(rElement, *mpDiscontinuousLevelSetVariable)) {

                // Retrieve shape functions values, gradients, integration weights and area normals (for positive interface) of the integrations points of the boundary inside the split element
                Matrix bd_Ns;
                typename ModifiedShapeFunctions::ShapeFunctionsGradientsType bd_DN_DXs;
                Vector bd_weights;
                std::vector<array_1d<double,3>> bd_pos_area_normals;
                GetDataForSplitElementBoundary(rElement, p_mod_sh_func_factory, bd_Ns, bd_DN_DXs, bd_weights, bd_pos_area_normals);

                // Calculate parent element size for the SBM BC imposition
                const double h = p_element_size_func(r_geom);

                // For each side of the boundary separately (positive and negative side of gamma), create a pointer vector with all the cloud nodes that affect the current element
                // To be used in the creation of the condition. Positive side refers to adding the negative node node's support cloud nodes.
                // NOTE that the obtained clouds are sorted by id to properly get the extension operator data //TODO: necessary???
                PointerVector<NodeType> cloud_nodes_vector_pos;
                PointerVector<NodeType> cloud_nodes_vector_neg;
                CreateCloudNodeVectorsForSplitElement(rElement, ext_op_map, cloud_nodes_vector_pos, cloud_nodes_vector_neg);

                // Iterate the interface integration/ Gauss pts.
                const std::size_t n_nodes = r_geom.PointsNumber();
                DenseVector<double> i_g_N;
                DenseMatrix<double> i_g_DN_DX;
                array_1d<double,3> i_g_coords;
                const std::size_t n_int_pts = bd_weights.size();
                for (std::size_t i_g = 0; i_g < n_int_pts; ++i_g) {

                    // Get element's shape function values, shape function derivatives, weight and positive area normal at Gauss pt.
                    i_g_N = row(bd_Ns, i_g);
                    i_g_DN_DX = bd_DN_DXs[i_g];
                    const double i_g_w = bd_weights[i_g];
                    const array_1d<double,3> i_g_pos_n = bd_pos_area_normals[i_g] / norm_2(bd_pos_area_normals[i_g]);

                    // Calculate Gauss pt. coordinates
                    noalias(i_g_coords) = ZeroVector(3);
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        noalias(i_g_coords) += i_g_N[i_node] * r_geom[i_node].Coordinates();
                    }

                    // Add Gauss pt. condition for positive side of boundary - using support cloud data for negative nodes
                    AddIntegrationPointCondition(rElement, h, ext_op_map, cloud_nodes_vector_pos,
                    i_g_coords, i_g_N, i_g_DN_DX, i_g_w,  i_g_pos_n, ++max_cond_id, /*ConsiderPositiveSide=*/true);

                    // Add Gauss pt. condition for negative side of boundary - using support cloud data for positive nodes
                    // NOTE that boundary normal is opposite
                    AddIntegrationPointCondition(rElement, h, ext_op_map, cloud_nodes_vector_neg,
                    i_g_coords, i_g_N, i_g_DN_DX, i_g_w, -i_g_pos_n, ++max_cond_id, /*ConsiderPositiveSide=*/false);
                }
            }
        }
    }

    void ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::SetInterfaceFlags()
    {
        // Initialize flags to false
        block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
            rNode.Set(ACTIVE, true);         // Nodes that belong to the elements to be assembled
            rNode.Set(BOUNDARY, false);      // Nodes that belong to the surrogate boundary
            rNode.Set(INTERFACE, false);     // Nodes that belong to the cloud of points
        });
        block_for_each(mpModelPart->Elements(), [](Element& rElement){
            rElement.Set(ACTIVE, true);      // Elements in the positive distance region (the ones to be assembled)
            rElement.Set(BOUNDARY, false);   // Intersected elements
            rElement.Set(INTERFACE, false);  // Positive distance elements owning the surrogate boundary nodes
        });

        // Find and deactivate BOUNDARY elements (gamma)
        for (auto& rElement : mpModelPart->Elements()) {
            if (IsSplit(rElement, *mpDiscontinuousLevelSetVariable)) {
                // Mark the intersected elements as BOUNDARY
                // The intersected elements are flagged as non ACTIVE to avoid assembling them
                // Note that the split elements BC is applied by means of the extension operators
                rElement.Set(ACTIVE, false);
                rElement.Set(BOUNDARY, true);
            }
        }

        // Find the surrogate boundary elements (gamma_tilde/ INTERFACE, not gamma/ BOUNDARY)
        // Note that we rely on the fact that the neighbors are sorted according to the faces
        for (auto& rElement : mpModelPart->Elements()) {
            if (rElement.Is(BOUNDARY)) {
                const auto& r_geom = rElement.GetGeometry();
                const std::size_t n_faces = r_geom.FacesNumber();
                auto& r_neigh_elems = rElement.GetValue(NEIGHBOUR_ELEMENTS);
                for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
                    // The neighbour corresponding to the current face is ACTIVE means that the current face is surrogate boundary
                    // Flag the current neighbour owning the surrogate face as INTERFACE
                    // The nodes will be flagged if required (MLS basis) when creating the cloud
                    auto p_neigh_elem = r_neigh_elems(i_face).get();
                    if (p_neigh_elem != nullptr) {
                        if (p_neigh_elem->Is(ACTIVE)) {
                            p_neigh_elem->Set(INTERFACE, true);
                        }
                    }
                }
            }
        }

        // Mark incised elements as BOUNDARY (they already are INTERFACE and ACTIVE) using ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED //TODO
        //TODO ?? Mark neighboring elements of incised elements also as BOUNDARY, so that the support clouds of both sides are separated better
        // --> PROBLEM that nodes that need support cloud could be only surrounded by BOUNDARY elements then and could not get support
        std::size_t n_elements_incised = 0;
        for (auto& rElement : mpModelPart->Elements()) {
            if (IsIncised(rElement, *mpDiscontinuousLevelSetVariable)) {
                rElement.Set(BOUNDARY, true);
                n_elements_incised++;
            }
        }
        KRATOS_WARNING("[ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility]") << "Number of incised elements: " << n_elements_incised << std::endl;
    }

    void ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::SetExtensionOperatorsForSplitElementNodes(
        NodesCloudMapType& rExtensionOperatorMap)
    {
        // Get the MLS shape functions function
        auto p_meshless_sh_func = GetMLSShapeFunctionsFunction();

        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto r_geom = rElement.GetGeometry();
            if (IsSplit(rElement, *mpDiscontinuousLevelSetVariable)) {
                // All nodes of intersected elements are either part of one side of the boundary or the other side of the boundary
                for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                    const auto p_node = r_geom(i_node);
                    // Calculate extension basis of the node if it has not already been done
                    const std::size_t found = rExtensionOperatorMap.count(p_node);
                    if (!found) {

                        // Get the support/ cloud nodes for the respective node
                        Matrix cloud_nodes_coordinates;
                        PointerVector<NodeType> cloud_nodes;
                        if (IsNegative(rElement, i_node, *mpDiscontinuousLevelSetVariable)) {
                            // Get the element's positive nodes and that sides neighboring nodes
                            std::vector<NodeType::Pointer> elem_pos_nodes;
                            for (std::size_t j_node = 0; j_node < r_geom.PointsNumber(); ++j_node) {
                                NodeType::Pointer p_node_j = r_geom(j_node);
                                if (!IsNegative(rElement, j_node, *mpDiscontinuousLevelSetVariable)) {
                                    elem_pos_nodes.push_back(p_node_j);
                                }
                            }
                            SetLateralSupportCloud(elem_pos_nodes, cloud_nodes, cloud_nodes_coordinates);
                        } else {
                            // Get the element's negative nodes and that sides neighboring nodes
                            std::vector<NodeType::Pointer> elem_neg_nodes;
                            for (std::size_t j_node = 0; j_node < r_geom.PointsNumber(); ++j_node) {
                                NodeType::Pointer p_node_j = r_geom(j_node);
                                if (IsNegative(rElement, j_node, *mpDiscontinuousLevelSetVariable)) {
                                    elem_neg_nodes.push_back(p_node_j);
                                }
                            }
                            SetLateralSupportCloud(elem_neg_nodes, cloud_nodes, cloud_nodes_coordinates);
                        }

                        // Calculate the extension basis in the current node (MLS shape functions)
                        Vector N_container;
                        const array_1d<double,3> r_coords = p_node->Coordinates();
                        const double kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
                        p_meshless_sh_func(cloud_nodes_coordinates, r_coords, kernel_rad, N_container);

                        // Save the extension operator nodal data to the extension operator map
                        std::size_t n_cl_nod = cloud_nodes.size();
                        CloudDataVectorType cloud_data_vector(n_cl_nod);
                        for (std::size_t i_cl_nod = 0; i_cl_nod < n_cl_nod; ++i_cl_nod) {
                            auto p_cl_node = cloud_nodes(i_cl_nod);
                            auto i_data = std::make_pair(p_cl_node, N_container[i_cl_nod]);
                            cloud_data_vector(i_cl_nod) = i_data;
                        }
                        auto ext_op_key_data = std::make_pair(p_node, cloud_data_vector);
                        rExtensionOperatorMap.insert(ext_op_key_data);
                    }
                }
            }
        }
    }

    void ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::SetLateralSupportCloud(
        const std::vector<NodeType::Pointer>& rSameSideNodes,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates)
    {
        // Find the same side support cloud of nodes
        // Note that we use an unordered_set to ensure that these are unique
        NodesCloudSetType aux_set;
        std::vector<NodeType::Pointer> cur_layer_nodes;
        std::vector<NodeType::Pointer> prev_layer_nodes;
        const std::size_t n_layers = mMLSExtensionOperatorOrder + 1;

        // Add first given nodes to cloud nodes set
        for (std::size_t i_node = 0; i_node < rSameSideNodes.size(); ++i_node) {
            const auto p_node = rSameSideNodes[i_node];
            aux_set.insert(p_node);
            prev_layer_nodes.push_back(p_node);
        }

        // Add several layers of neighbors to the cloud nodes
        // Note that we check the order of the MLS interpolation to add nodes from enough interior layers
        for (std::size_t i_layer = 0; i_layer < n_layers; ++i_layer) {

            // Find elemental neighbors of the nodes of the previous layer
            for (auto& p_node : prev_layer_nodes) {
                auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);

                // Add all nodes of neighboring elements to cloud nodes set if element is not (!) BOUNDARY
                // This way only nodes of the same side as the given nodes are added
                // NOTE that 'not BOUNDARY' is used instead of 'not ACTIVE' so that incised elements, which might be active, also represent a boundary
                for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                    auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                    if (p_elem_neigh != nullptr && !p_elem_neigh->Is(BOUNDARY)) {
                        const auto& r_geom = p_elem_neigh->GetGeometry();
                        for (std::size_t i_neigh_node = 0; i_neigh_node < r_geom.PointsNumber(); ++i_neigh_node) {
                            NodeType::Pointer p_neigh = r_geom(i_neigh_node);
                            p_neigh->Set(INTERFACE, true);
                            aux_set.insert(p_neigh);
                            cur_layer_nodes.push_back(p_neigh);
                        }
                    }
                }
            }
            prev_layer_nodes = cur_layer_nodes;
            cur_layer_nodes.clear();
        }

        // Check that the current nodes are enough to perform the MLS calculation
        // If not sufficient, add the nodal neighbors of the current nodes to the cloud of points
        const std::size_t n_cloud_nodes_temp = aux_set.size();
        KRATOS_ERROR_IF(n_cloud_nodes_temp == 0) << "Degenerated case with no neighbors. Check if the node " << rSameSideNodes[0]->Id() << " belongs to an intersected and isolated element." << std::endl;
        if (n_cloud_nodes_temp < GetRequiredNumberOfPoints()) {
            NodesCloudSetType aux_extra_set(aux_set);
            for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
                auto& r_elem_neigh_vect = (*it_set)->GetValue(NEIGHBOUR_ELEMENTS);
                for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                    auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                    if (p_elem_neigh != nullptr && !p_elem_neigh->Is(BOUNDARY)) {
                        const auto& r_geom = p_elem_neigh->GetGeometry();
                        for (std::size_t i_neigh_node = 0; i_neigh_node < r_geom.PointsNumber(); ++i_neigh_node) {
                            NodeType::Pointer p_neigh = r_geom(i_neigh_node);
                            p_neigh->Set(INTERFACE, true);
                            aux_extra_set.insert(p_neigh);
                        }
                    }
                }
            }
            aux_set = aux_extra_set;
        }

        // Sort the obtained nodes by id
        const std::size_t n_cloud_nodes = aux_set.size();
        rCloudNodes.resize(n_cloud_nodes);
        std::size_t aux_i = 0;
        for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
            rCloudNodes(aux_i++) = *it_set;
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

    void ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::GetDataForSplitElementBoundary(
        const ElementType& rElement,
        ModifiedShapeFunctionsFactoryType pModifiedShapeFunctionsFactory,
        Matrix& rBoundaryShapeFunctionValues,
        ModifiedShapeFunctions::ShapeFunctionsGradientsType& rBoundaryShapeFunctionDerivatives,
        Vector& rBoundaryWeights,
        std::vector<array_1d<double,3>>& rBoundaryAreaNormals)
    {
        // Set up the distances vector
        const GeometryType::Pointer p_geom = rElement.pGetGeometry();
        const std::size_t n_nodes = p_geom->PointsNumber();
        Vector nodal_distances(n_nodes);
        SetNodalDistancesVector(rElement, *mpDiscontinuousLevelSetVariable, nodal_distances);

        // Set the modified shape functions pointer and calculate the positive interface data
        auto p_mod_sh_func = pModifiedShapeFunctionsFactory(p_geom, nodal_distances);
        //TODO: Add a method without the interface gradients
        p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(rBoundaryShapeFunctionValues, rBoundaryShapeFunctionDerivatives, rBoundaryWeights, GeometryData::IntegrationMethod::GI_GAUSS_2);
        p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(rBoundaryAreaNormals, GeometryData::IntegrationMethod::GI_GAUSS_2);
    }

    void ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::CreateCloudNodeVectorsForSplitElement(
        const ElementType& rElement,
        NodesCloudMapType& rExtensionOperatorMap,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide)
    {
        // Create an auxiliary set with all the cloud nodes that affect the current element for each side separately
        NodesCloudSetType cloud_nodes_set_pos;
        NodesCloudSetType cloud_nodes_set_neg;
        const auto& r_geom = rElement.GetGeometry();
        for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            NodeType::Pointer p_node = r_geom(i_node);
            if (IsNegative(rElement, i_node, *mpDiscontinuousLevelSetVariable)) {
                // Add negative side node to cloud nodes set of negative side of the boundary //TODO: really necessary? these nodes should be there already as other sides cloud nodes?!
                cloud_nodes_set_neg.insert(p_node);
                // Add negative side's node's cloud nodes to cloud nodes set of positive side of the boundary
                auto& r_ext_op_data = rExtensionOperatorMap[p_node];
                for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                    auto& p_node = std::get<0>(*it_data);
                    cloud_nodes_set_pos.insert(p_node);
                }
            } else {
                // Add positive side node to cloud nodes set of positive side of the boundary //TODO: really necessary? these nodes should be there already as other sides cloud nodes?!
                cloud_nodes_set_pos.insert(p_node);
                // Add positive side's node's cloud nodes to cloud nodes set of negative side of the boundary
                auto& r_ext_op_data = rExtensionOperatorMap[p_node];
                for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                    auto& p_node = std::get<0>(*it_data);
                    cloud_nodes_set_neg.insert(p_node);
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

        // Sort obtained cloud node vectors by ID to properly get the extension operator data
        std::sort(rCloudNodeVectorPositiveSide.ptr_begin(), rCloudNodeVectorPositiveSide.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});
        std::sort(rCloudNodeVectorNegativeSide.ptr_begin(), rCloudNodeVectorNegativeSide.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});
    }

    void ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility::AddIntegrationPointCondition(
        const ElementType& rElement,
        const double ElementSize,
        NodesCloudMapType& rExtensionOperatorMap,
        const PointerVector<NodeType>& rCloudNodeVector,
        const array_1d<double,3>& rIntPtCoordinates,
        const DenseVector<double>& rIntPtShapeFunctionValues,
        const DenseMatrix<double>& rIntPtShapeFunctionDerivatives,
        const double rIntPtWeight,
        const array_1d<double,3>& rIntPtNormal,
        const std::size_t ConditionId,
        bool ConsiderPositiveSide)
    {
        const auto& r_geom = rElement.GetGeometry();

        // Initialize the extension operator containers
        const std::size_t n_cl_nodes = rCloudNodeVector.size();
        const std::size_t n_dim = r_geom.WorkingSpaceDimension();
        Vector N_container = ZeroVector(n_cl_nodes);
        Matrix DN_DX_container = ZeroMatrix(n_cl_nodes, n_dim);

        // Loop the nodes that are involved in the current element
        for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            const auto& r_node = r_geom[i_node];
            if (ConsiderPositiveSide != IsNegative(rElement, i_node, *mpDiscontinuousLevelSetVariable)) {
                // If positive side is considered and node is positive OR negative side is considered and node is negative, then add the standard shape function contribution
                // Note that we need to check for the ids to match in the geometry as nodes in the map are mixed
                for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                    auto& p_cl_node = rCloudNodeVector(i_cl);
                    if (r_node.Id() == p_cl_node->Id()) {
                        N_container(i_cl) += rIntPtShapeFunctionValues(i_node);
                        for (std::size_t d = 0; d < n_dim; ++d) {
                            DN_DX_container(i_cl,d) += rIntPtShapeFunctionDerivatives(i_node,d);
                        }
                        break;
                    }
                }
            } else {
                // Get the weight as the corresponding nodal shape function value
                const double i_node_N = rIntPtShapeFunctionValues(i_node);
                const auto i_node_grad_N = row(rIntPtShapeFunctionDerivatives, i_node);

                // If it is not ACTIVE (negative side) search for the extension operator data
                auto p_node = r_geom(i_node);

                auto& ext_op_data = rExtensionOperatorMap[p_node];

                // Loop the current negative node to get its extrapolation operator data and apply the cloud node weight (i_cl_node_N) to make the basis conformant
                // Note that we need to check for the ids to match in the geometry as nodes in the map are mixed
                for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
                    auto& r_node_data = *it_data;
                    std::size_t data_node_id = (std::get<0>(r_node_data))->Id();
                    for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                        auto& p_cl_node = rCloudNodeVector(i_cl);
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

        // Create a new condition with a geometry made up with the basis nodes
        auto p_prop = rElement.pGetProperties();
        auto p_cond = mpConditionPrototype->Create(ConditionId, rCloudNodeVector, p_prop);
        p_cond->Set(ACTIVE, true);
        mpBoundarySubModelPart->AddCondition(p_cond);

        // Store the SBM BC data in the condition database
        p_cond->SetValue(ELEMENT_H, ElementSize);
        p_cond->SetValue(NORMAL, rIntPtNormal);
        p_cond->SetValue(INTEGRATION_WEIGHT, rIntPtWeight);
        p_cond->SetValue(INTEGRATION_COORDINATES, rIntPtCoordinates);
        p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
        p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);
    }

} // namespace Kratos.
