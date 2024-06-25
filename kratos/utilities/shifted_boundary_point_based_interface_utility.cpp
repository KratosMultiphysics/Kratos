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
#include "includes/define.h"
#include "includes/element.h"
#include "includes/smart_pointers.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/element_size_calculator.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/shifted_boundary_point_based_interface_utility.h"

#include "utilities/binbased_fast_point_locator.h"


namespace Kratos
{
    ShiftedBoundaryPointBasedInterfaceUtility::ShiftedBoundaryPointBasedInterfaceUtility(
        Model& rModel,
        Parameters ThisParameters)
    {
        // Validate input settings with defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Retrieve the required model parts
        const std::string model_part_name = ThisParameters["model_part_name"].GetString();
        const std::string skin_model_part_name = ThisParameters["skin_model_part_name"].GetString();
        const std::string boundary_sub_model_part_name = ThisParameters["boundary_sub_model_part_name"].GetString();
        mpModelPart = &rModel.GetModelPart(model_part_name);
        //TODO error message if skin ModelPart does not exist - automatically?
        mpSkinModelPart = &rModel.GetModelPart(skin_model_part_name);
        if (mpModelPart->HasSubModelPart(boundary_sub_model_part_name)) {
            mpBoundarySubModelPart = &(mpModelPart->GetSubModelPart(boundary_sub_model_part_name));
            KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceProcess", mpBoundarySubModelPart->NumberOfNodes() != 0) << "Provided SBM model part has nodes." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceProcess", mpBoundarySubModelPart->NumberOfElements() != 0) << "Provided SBM model part has elements." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceProcess", mpBoundarySubModelPart->NumberOfConditions() != 0) << "Provided SBM model part has conditions." << std::endl;
        } else {
            mpBoundarySubModelPart = &(mpModelPart->CreateSubModelPart(boundary_sub_model_part_name));
        }

        // Set the order of the MLS extension operator used in the MLS shape functions utility
        mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();

        // If true, the basis is created such that it is conforming with the linear FE space of the surrogate boundary
        mConformingBasis = ThisParameters["conforming_basis"].GetBool();

        // If true, the basis is created such that the surrogate boundary gradient is kept
        const std::string ext_op_type = ThisParameters["extension_operator_type"].GetString();
        if (ext_op_type == "MLS") {
            mExtensionOperator = ExtensionOperator::MLS;
        } else {
            KRATOS_ERROR << "Wrong 'extension_operator_type' provided. Only 'MLS' 'extension_operator_type' is supported by point based shifted boundary interface utility." << std::endl;
        }

        // Set the SBD condition prototype to be used in the condition creation
        std::string interface_condition_name = ThisParameters["sbm_interface_condition_name"].GetString();
        KRATOS_ERROR_IF(interface_condition_name == "") << "SBM interface condition has not been provided." << std::endl;
        mpConditionPrototype = &KratosComponents<Condition>::Get(interface_condition_name);
    };

    void ShiftedBoundaryPointBasedInterfaceUtility::AddSkinIntegrationPointConditions()
    {
        CalculateMeshlessBasedConformingExtensionBasis();
    }

    const Parameters ShiftedBoundaryPointBasedInterfaceUtility::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "skin_model_part_name" : "",
            "boundary_sub_model_part_name" : "",
            "sbm_interface_condition_name" : "",
            "conforming_basis" : true,
            "extension_operator_type" : "MLS",
            "mls_extension_operator_order" : 1
        })" );

        return default_parameters;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateMeshlessBasedConformingExtensionBasis()
    {
        // Check that the volume model part has elements
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements())
            << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;

        // Map skin points (boundary) to elements of volume model part
        SkinPointsToElementsMapType skin_points_map;
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        switch (n_dim) {
            case 2:
                MapSkinPointsToElements<2>(skin_points_map);
            case 3:
                MapSkinPointsToElements<3>(skin_points_map);
            default:
                KRATOS_ERROR << "Wrong domain size.";
        }

        // Set the required interface flags (ACTIVE, BOUNDARY, INTERFACE)
        SetInterfaceFlags(skin_points_map);

        // Iterate over the split elements to create a vector for each element defining both sides using the element's skin points normals and positions
        // The resulting vector is as long as the number of nodes of the element and a positive value stands for the positive side of the boundary, a negative one for the negative side.
        SidesVectorToElementsMapType sides_vector_map;
        SetSidesVectorsForSplitElements(skin_points_map, sides_vector_map);

        // Iterate over the split elements to create an extension basis for each node of the element (MLS shape functions values for support cloud of node)
        NodesCloudMapType ext_op_map;
        SetExtensionOperatorsForSplitElementNodes(sides_vector_map, ext_op_map);

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        auto p_element_size_func = GetElementSizeFunction(mpModelPart->ElementsBegin()->GetGeometry());

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Create the interface conditions
        //TODO: THIS CAN BE PARALLEL (WE JUST NEED TO MAKE CRITICAL THE CONDITION ID UPDATE)
        for (const auto& [p_element, skin_points_data_vector]: skin_points_map) {
            const auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            
            // For each side of the boundary separately (positive and negative side of gamma), create a pointer vector with all the cloud nodes that affect the current element
            // To be used in the creation of the condition. Positive side refers to adding the negative node's support cloud nodes.
            // NOTE that the obtained clouds are sorted by id to properly get the extension operator data //TODO: necessary?
            PointerVector<NodeType> cloud_nodes_vector_pos;
            PointerVector<NodeType> cloud_nodes_vector_neg;
            CreateCloudNodeVectorsForSplitElement(*p_element, ext_op_map, cloud_nodes_vector_pos, cloud_nodes_vector_neg, sides_vector_map[p_element]);

            // Calculate parent element size for the SBM BC imposition
            const double h = p_element_size_func(r_geom);

            // Get vector defining positive and negative side of the boundary
            auto& sides_vector = sides_vector_map[p_element];

            // Iterate over the element's skin points adding a positive side and a negative side condition for each skin point
            for (std::size_t i_skin_pt = 0; i_skin_pt < skin_points_data_vector.size(), ++i_skin_pt) {
                // Get the skin point's position and area normal (integration point of the boundary)
                auto skin_pt_data = skin_points_data_vector[i_skin_pt];
                const array_1d<double, 3> skin_pt_position = std::get<0>(skin_points_data);
                const array_1d<double, 3> skin_pt_area_normal = std::get<1>(skin_points_data);
 
                // Get the split element's shape function values and derivatives at the skin/ integration point
                Vector skin_pt_N(n_nodes);
                Matrix skin_pt_DNDX = ZeroMatrix(n_nodes, n_dim);
                GetDataForSplitElementIntegrationPoint(*p_element, skin_pt_N, skin_pt_DNDX);

                // Add skin pt. condition for positive side of boundary - using support cloud data for negative nodes
                AddIntegrationPointCondition(*p_element, sides_vector, h, skin_pt_position, skin_pt_area_normal, 
                ext_op_map, cloud_nodes_vector_pos, skin_pt_N, skin_pt_DNDX, ++max_cond_id, /*ConsiderPositiveSide=*/true);

                // Add skin pt. condition for negative side of boundary - using support cloud data for positive nodes
                // NOTE that boundary normal is opposite - TODO correct side????
                AddIntegrationPointCondition(*p_element, sides_vector, h, skin_pt_position, -skin_pt_area_normal, 
                ext_op_map, cloud_nodes_vector_neg, skin_pt_N, skin_pt_DNDX, ++max_cond_id, /*ConsiderPositiveSide=*/false);
            }
        }
    }

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedInterfaceUtility::MapSkinPointsToElements(
        SkinPointsToElementsMapType& rSkinPointsMap)
    {
        // Check that the skin model part has elements
        KRATOS_ERROR_IF_NOT(mpSkinModelPart->NumberOfElements())
            << "There are no elements in skin model part (boundary) '" << mpSkinModelPart->FullName() << "'." << std::endl;

        // Set the bin-based fast point locator utility
        BinBasedFastPointLocator<TDim> bin_based_point_locator(*mpModelPart);
        bin_based_point_locator.UpdateSearchDatabase();

        // Search the skin points (skin model part integration points) in the volume mesh elements
        // TODO get the skin model part integration points instead!!!
        block_for_each(mpSkinModelPart->Nodes(), 
            typename BinBasedFastPointLocator<TDim>::ResultContainerType(mSearchMaxResults),
            [&](auto& rNode, auto& rSearchResults){

            "search_tolerance" : 1.0e-5,
            "search_max_results" : 10000,

            Vector aux_N(TDim+1);
            Element::Pointer p_elem = nullptr;
            const bool is_found = bin_based_point_locator.FindPointOnMesh(
                rNode.Coordinates(),
                aux_N,
                p_elem,
                rSearchResults.begin(),
                mSearchMaxResults,
                mSearchTolerance);

            // Check if the node is found
            if (is_found){
                //TODO
                //TODO get integration points with normal and weight

                std::size_t n_skin_points = 1;
                SkinPointsDataVectorType skin_points_data_vector(n_skin_points);
                for (std::size_t i_skin_pt = 0; i_skin_pt < n_skin_points; ++i_skin_pt) {
                    //auto p_cl_node = cloud_nodes(i_skin_pt);
                    const array_1d<double, 3> skin_pt_position(3, 0.0);
                    const array_1d<double, 3> skin_pt_area_normal(3, 0.0);
                    auto i_data = std::make_pair(skin_pt_position, skin_pt_area_normal);
                    skin_points_data_vector(i_skin_pt) = i_data;
                }

                auto skin_points_key_data = std::make_pair(&rElement, skin_points_data_vector);
                rSkinPointsMap.insert(skin_points_key_data);

            } else {
                //TODO: coordinates of skin point!!
                KRATOS_WARNING("FixedMeshALEUtilities")
                    << "A skin point has not been found in any volume model part element. Point coordinates: (" << rNode.X() << " , " << rNode.Y() << " , " << rNode.Z() << ")" << std::endl;
            }
        });
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetInterfaceFlags(
        const SkinPointsToElementsMapType& rSkinPointsMap)
    {
        // Initialize flags to false
        block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
            //rNode.Set(ACTIVE, true);         // Nodes that belong to the elements to be assembled
            //rNode.Set(BOUNDARY, false);      // Nodes that belong to the surrogate boundary
            rNode.Set(INTERFACE, false);     // Nodes that belong to the cloud of points
        });
        block_for_each(mpModelPart->Elements(), [](Element& rElement){
            rElement.Set(ACTIVE, true);      // Elements in the positive distance region (the ones to be assembled)
            rElement.Set(BOUNDARY, false);   // Intersected elements
            rElement.Set(INTERFACE, false);  // Positive distance elements owning the surrogate boundary nodes
        });

        // Deactivate split elements and mark them as BOUNDARY (gamma)
        // Note that the split elements BC is applied by means of the extension operators
        /*for (const auto& [p_element, _]: rSkinPointsMap) {
            p_element->Set(ACTIVE, false);
            p_element->Set(BOUNDARY, true);
        }*/
        //TODO: parallel
        block_for_each(rSkinPointsMap, [](std::pair<ElementType::Pointer, SkinPointsDataVectorType> element_skin_points){
            element_skin_points.first->Set(ACTIVE, false);
            element_skin_points.first->Set(BOUNDARY, true);
        });

        // Find the surrogate boundary elements and mark them as INTERFACE (gamma_tilde)
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
    }

    void SetSidesVectorsForSplitElements(
        const SkinPointsToElementsMapType& rSkinPointsMap,
        SidesVectorToElementsMapType& rSidesVectorMap)
    {
        // TODO can be parallel if insert into rSidesVectorMap is unique?
        for (const auto& [p_element, skin_points_data_vector]: rSkinPointsMap) {
            // Calculate average position and normal of the element's skin points
            std::size_t n_skin_points = 0;
            const array_1d<double, 3> avg_position(3, 0.0);
            const array_1d<double, 3> avg_normal(3, 0.0);
            for (const auto& skin_points_data: skin_points_data_vector) {
                n_skin_points++;
                avg_position += std::get<0>(skin_points_data);
                avg_normal += std::get<1>(skin_points_data);
            }
            //avg_normal /= norm_2(avg_normal);
            avg_position /= n_skin_points;

            
            // Compute the dot product for each node between a vector from the average skin position to the node and the skin's average normal
            // Note that for a positive dot product the node is saved as being on the positive side of the boundary, negative dot product = negative side
            Vector sides_vector;
            //auto& r_geom = p_element->pGetGeometry()
            //for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            //const auto& r_node = r_geom[i_node];
            for (const auto& r_node : p_element->pGetGeometry()) {
                array_1d<double, 3> skin_pt_to_node = avg_position - r_node.Coordinates();
                if (inner_prod(skin_pt_to_node, avg_normal) > 0.0) {
                    sides_vector.push_back(1.0);
                } else {
                    sides_vector.push_back(-1.0);
                }
            }

            auto sides_vector_key_data = std::make_pair(p_element, sides_vector);
            rSidesVectorMap.insert(skin_points_key_data);
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetExtensionOperatorsForSplitElementNodes(
        const SidesVectorToElementsMapType& rSidesVectorMap,
        NodesCloudMapType& rExtensionOperatorMap)
    {
        // Get the MLS shape functions function
        auto p_meshless_sh_func = GetMLSShapeFunctionsFunction();

        for (const auto& [p_element, sides_vector]: rSidesVectorMap) {
            const auto r_geom = p_element->pGetGeometry();
            // All nodes of intersected elements are either part of one side of the boundary or the other side of the boundary
            for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                const auto p_node = r_geom(i_node);
                // Calculate extension basis of the node if it has not already been done
                const std::size_t found = rExtensionOperatorMap.count(p_node);
                if (!found) {

                    // Get the support/ cloud nodes for the respective node
                    //TODO only once per side necessary, currently done multiple times (for each node on each side)
                    Matrix cloud_nodes_coordinates;
                    PointerVector<NodeType> cloud_nodes;

                    // Node is on positive side of the boundary
                    if (sides_vector[i_node] > 0.0) {
                        // Get the element's negative nodes and that sides neighboring nodes
                        std::vector<NodeType::Pointer> elem_neg_nodes;
                        for (std::size_t j_node = 0; j_node < r_geom.PointsNumber(); ++j_node) {
                            if (sides_vector[j_node] <= 0.0) {
                                elem_neg_nodes.push_back(r_geom(j_node));
                            }
                        }
                        SetLateralSupportCloud(elem_neg_nodes, cloud_nodes, cloud_nodes_coordinates);
                
                    // Node is on negative side of the boundary
                    } else {
                        // Get the element's positive nodes and that sides neighboring nodes
                        std::vector<NodeType::Pointer> elem_pos_nodes;
                        for (std::size_t j_node = 0; j_node < r_geom.PointsNumber(); ++j_node) {
                            if (sides_vector[j_node] > 0.0) {
                                elem_pos_nodes.push_back(r_geom(j_node));
                            }
                        }
                        SetLateralSupportCloud(elem_pos_nodes, cloud_nodes, cloud_nodes_coordinates);
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

    void ShiftedBoundaryPointBasedInterfaceUtility::SetLateralSupportCloud(
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

        // Add first given nodes to cloud nodes set (1st layer of sampling points)
        for (std::size_t i_node = 0; i_node < rSameSideNodes.size(); ++i_node) {
            const auto p_node = rSameSideNodes[i_node];
            aux_set.insert(p_node);
            prev_layer_nodes.push_back(p_node);
        }

        // Add several layers of neighbors to the cloud nodes
        // NOTE that we check the order of the MLS interpolation to add nodes from enough interior 
        // NOTE that we start from 1 here as the first layer already has been added
        for (std::size_t i_layer = 1; i_layer < n_layers; ++i_layer) {

            // Find elemental neighbors of the nodes of the previous layer
            for (auto& p_node : prev_layer_nodes) {
                auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);

                // Add all nodes of neighboring elements to cloud nodes set if element is not (!) BOUNDARY
                // This way only nodes of the same side as the given nodes are added
                // NOTE that 'not BOUNDARY' is used instead of 'not ACTIVE' so that incised elements, which might be active, also represent a boundary
                for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                    auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                    // TODO add points only if they are on the correct side of the intersection plane? (for intersected)
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
            KRATOS_WARNING("[SHIFTED-BOUNDARY-POINT-BASED-DISCONTINUOUS-INTERFACE-UTILITY] Extra set of points needed for MLS calculation!!");
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

    void ShiftedBoundaryPointBasedInterfaceUtility::CreateCloudNodeVectorsForSplitElement(
        const ElementType& rElement,
        const Vector& rSidesVector,
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
            if (rSidesVector(i_node) > 0.0) {
                // Add positive side node to cloud nodes set of positive side of the boundary //TODO: really necessary? these nodes should be there already as other sides cloud nodes?!
                cloud_nodes_set_pos.insert(p_node);
                // Add positive side's node's cloud nodes to cloud nodes set of negative side of the boundary
                auto& r_ext_op_data = rExtensionOperatorMap[p_node];
                for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                    auto& p_node = std::get<0>(*it_data);
                    cloud_nodes_set_neg.insert(p_node);
                }
            } else {
                // Add negative side node to cloud nodes set of negative side of the boundary //TODO: really necessary? these nodes should be there already as other sides cloud nodes?!
                cloud_nodes_set_neg.insert(p_node);
                // Add negative side's node's cloud nodes to cloud nodes set of positive side of the boundary
                auto& r_ext_op_data = rExtensionOperatorMap[p_node];
                for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                    auto& p_node = std::get<0>(*it_data);
                    cloud_nodes_set_pos.insert(p_node);
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

    void GetDataForSplitElementIntegrationPoint(
        const ElementType& rElement,
        const array_1d<double,3>& rIntPtCoordinates,
        Vector& rIntPtShapeFunctionValues, 
        Matrix& rIntPtShapeFunctionDerivatives)
    {
        const auto& r_geom = rElement.GetGeometry();

        // Compute the local coordinates of the integration point in the element's geometry
        //Geometry<Node>::CoordinatesArrayType aux_global_coords = ZeroVector(3);
        //r_bd_geom.GlobalCoordinates(aux_global_coords, r_integration_points[i_int_pt].Coordinates());
        Geometry<Node>::CoordinatesArrayType int_pt_local_coords = ZeroVector(3);
        rElement.PointLocalCoordinates(int_pt_local_coords, rIntPtCoordinates);

        // Get N of the element at the integration point
        rElement.ShapeFunctionsValues(rIntPtShapeFunctionValues, int_pt_local_coords);

        // Get DN_DX of the element at the integration point
        Matrix aux_DN_DXi_parent, aux_J_parent, aux_J_inv_parent;
        double aux_detJ_parent;
        rElement.ShapeFunctionsLocalGradients(aux_DN_DXi_parent, int_pt_local_coords);
        rElement.Jacobian(aux_J_parent, int_pt_local_coords);
        MathUtils<double>::InvertMatrix(aux_J_parent, aux_J_inv_parent, aux_detJ_parent);
        rIntPtShapeFunctionDerivatives = prod(aux_DN_DXi_parent, aux_J_inv_parent);
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::AddIntegrationPointCondition(
        const ElementType& rElement,
        const Vector& rSidesVector,
        const double ElementSize,
        const array_1d<double,3>& rIntPtCoordinates,
        const array_1d<double,3>& rIntPtAreaNormal,
        NodesCloudMapType& rExtensionOperatorMap,
        const PointerVector<NodeType>& rCloudNodeVector,
        const Vector& rIntPtShapeFunctionValues,
        const Matrix& rIntPtShapeFunctionDerivatives,
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
            // If node is on the side that is being considered, then add the standard shape function contribution of the node at the position of the skin point
            if (ConsiderPositiveSide != (rSidesVector[i_node] <= 0.0)) {
                // Note that we need to check for the ids to match in the geometry as nodes in the map are mixed
                for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                    auto& p_cl_node = rCloudNodeVector(i_cl);
                    if (r_node.Id() == p_cl_node->Id()) {
                        N_container(i_cl) += rIntPtShapeFunctionValues(i_node);
                        for (std::size_t d = 0; d < n_dim; ++d) {
                            DN_DX_container(i_cl, d) += rIntPtShapeFunctionDerivatives(i_node, d);
                        }
                        break;
                    }
                }
            // If node is on the other side of the boundary, then get its shape function values and derivatives for the skin point and its extension operator data
            } else {
                // Get the weight as the corresponding nodal shape function value of the node at the position of the skin point
                const double i_node_N = rIntPtShapeFunctionValues(i_node);
                const auto i_node_grad_N = row(rIntPtShapeFunctionDerivatives, i_node);

                // Get the node's extension operator data
                auto& ext_op_data = rExtensionOperatorMap[r_geom(i_node)];

                // Iterate over the node's extension operator data and apply the cloud node weight (i_cl_node_N) to make the basis conformant
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
        p_cond->SetValue(INTEGRATION_COORDINATES, rIntPtCoordinates);
        const double skin_pt_weight = norm_2(rIntPtAreaNormal);  // TODO: not necessary? done in the condition anyway?!
        p_cond->SetValue(NORMAL, rIntPtAreaNormal/ skin_pt_weight);
        p_cond->SetValue(INTEGRATION_WEIGHT, skin_pt_weight);
        p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
        p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);
    }

    ShiftedBoundaryPointBasedInterfaceUtility::MeshlessShapeFunctionsFunctionType ShiftedBoundaryPointBasedInterfaceUtility::GetMLSShapeFunctionsFunction() const
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

    ShiftedBoundaryPointBasedInterfaceUtility::ElementSizeFunctionType ShiftedBoundaryPointBasedInterfaceUtility::GetElementSizeFunction(const GeometryType& rGeometry)
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

    double ShiftedBoundaryPointBasedInterfaceUtility::CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin)
    {
        const std::size_t n_nodes = rCloudCoordinates.size1();
        const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
            return std::pow(rCloudCoordinates(I,0) - rOrigin(0),2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1),2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2),2);
        });
        return std::sqrt(squared_rad);
    }

    std::size_t ShiftedBoundaryPointBasedInterfaceUtility::GetRequiredNumberOfPoints()
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

    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::MapSkinPointsToElements<2>(SkinPointsToElementsMapType& rSkinPointsMap);
    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::MapSkinPointsToElements<3>(SkinPointsToElementsMapType& rSkinPointsMap);

} // namespace Kratos.
