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
#include "includes/kratos_flags.h"
#include "includes/smart_pointers.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
#include "input_output/logger.h"
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
#include "processes/find_intersected_geometrical_objects_process.h"

#include <boost/numeric/ublas/vector_expression.hpp>
#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>


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
        mSkinModelPartName = ThisParameters["skin_model_part_name"].GetString();
        const std::string boundary_sub_model_part_name = ThisParameters["boundary_sub_model_part_name"].GetString();
        // Get skin model part
        mpSkinModelPart = &rModel.GetModelPart(mSkinModelPartName);
        // Get volume model part and check that the volume model part has elements
        mpModelPart = &rModel.GetModelPart(model_part_name);
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements()) << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;
        // Get and check or create boundary model part as sub model part of the volume model part
        if (mpModelPart->HasSubModelPart(boundary_sub_model_part_name)) {
            mpBoundarySubModelPart = &(mpModelPart->GetSubModelPart(boundary_sub_model_part_name));
            KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", mpBoundarySubModelPart->NumberOfNodes() != 0) << "Provided SBM model part has nodes." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", mpBoundarySubModelPart->NumberOfElements() != 0) << "Provided SBM model part has elements." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", mpBoundarySubModelPart->NumberOfConditions() != 0) << "Provided SBM model part has conditions." << std::endl;
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

        // Set the shifted-boundary condition prototype to be used in the condition creation
        std::string interface_condition_name = ThisParameters["sbm_interface_condition_name"].GetString();
        KRATOS_ERROR_IF(interface_condition_name == "") << "SBM interface condition has not been provided." << std::endl;
        mpConditionPrototype = &KratosComponents<Condition>::Get(interface_condition_name);

        // Set which side of the skin model part is considered
        const std::string active_side_of_skin = ThisParameters["active_side_of_skin"].GetString();
        if (active_side_of_skin != "positive" && active_side_of_skin != "negative" && active_side_of_skin != "both") {
            KRATOS_ERROR << "Unknown 'active_side_of_skin' given. 'positive', 'negative' or 'both' sides are supported by point based shifted boundary interface utility." << std::endl;
        }
        //TODO actually deactivate side that is not considered --> flood fill and majority vote like Manuel
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", active_side_of_skin != "both") << "So far always both sides of the skin are being considered." << std::endl;
        //TODO: actually being used?? - No extension operators will be calculated and no wall conditions will be created for inactive side
        if (active_side_of_skin == "positive") {
            mNegativeSideIsActive = false;
        } else if (active_side_of_skin == "negative") {
            mPositiveSideIsActive = false;
        }

        // Set which side of the skin model part is an enclosed area
        //TODO NEXT use flood fill (--> Manuel) and Octree to check and search for enclosed areas for setting the pressure!! (robustness)
        const std::string enclosed_area = ThisParameters["enclosed_area"].GetString();
        if (enclosed_area != "positive" && enclosed_area != "negative" && enclosed_area != "none") {
            KRATOS_ERROR << "Unknown 'enclosed_area' keyword given. 'positive' or 'negative' side or 'none' are supported by point based shifted boundary interface utility." << std::endl;
        }
        if (enclosed_area == "positive") {
            mPositiveSideIsEnclosed = true;
        } else if (enclosed_area == "negative") {
            mNegativeSideIsEnclosed = true;
        }

        // If true, then the support cloud for a node can be found beyond a first layer of deactivated BOUNDARY elements (direct neighbors) in normal direction
        // NOTE that this is beneficial for small cuts and overlapping boundaries
        mCrossBoundaryNeighbors = ThisParameters["cross_boundary_neighbors"].GetBool();

        // If true, then elements will be declared BOUNDARY which are intersected by the tessellated skin geometry
        mUseTessellatedBoundary = ThisParameters["use_tessellated_boundary"].GetBool();
    };

    void ShiftedBoundaryPointBasedInterfaceUtility::ResetFlags()
    {
        // Check that the volume model part has elements
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements())
            << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;

        // Activate all elements and initialize flags to false
        // NOTE Resetting the interface flags will eliminate previously embedded model parts except for their wall conditions!
        block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
            rNode.Set(ACTIVE, true);         // Nodes that belong to the elements to be assembled
            rNode.Set(BOUNDARY, false);      // TODO Nodes that belong to the surrogate boundary
            rNode.Set(INTERFACE, false);     // TODO Nodes that belong to the cloud of points
        });
        block_for_each(mpModelPart->Elements(), [](Element& rElement){
            rElement.Set(ACTIVE, true);      // Elements in the positive distance region (the ones to be assembled)
            rElement.Set(BOUNDARY, false);   // Elements in which the skin geometry is located (true boundary gamma)
            rElement.Set(INTERFACE, false);  // Positive distance elements owning the surrogate boundary nodes
        });

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Boundary and interface flags were reset and all elements and nodes re-activated." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetTessellatedBoundaryFlags()
    {
        // Reset SELECTED flag
        block_for_each(mpModelPart->Elements(), [](Element& rElement){
            rElement.Set(SELECTED, false);
        });

        // Create and initialize find_intersected_objects_process
        Flags find_intersected_options;
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_CONDITIONS, false);
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_ELEMENTS, true);
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS, false);
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS, true);
        FindIntersectedGeometricalObjectsProcess find_intersected_objects_process = FindIntersectedGeometricalObjectsProcess(*mpModelPart, *mpSkinModelPart, find_intersected_options);
        find_intersected_objects_process.ExecuteInitialize();

        // Find intersections
        find_intersected_objects_process.FindIntersections();  // FindIntersectedSkinObjects - NOTE that this marks intersected elements as SELECTED
        std::vector<PointerVector<GeometricalObject>>&  r_intersected_objects = find_intersected_objects_process.GetIntersections();

        // Get elements in the same order as find_intersected_objects_process
        const std::size_t n_elements = (find_intersected_objects_process.GetModelPart1()).NumberOfElements();
        auto& r_elements = (find_intersected_objects_process.GetModelPart1()).ElementsArray();

        // Check for nodes that would be surrounded by BOUNDARY elements only, as these nodes are probably intersected by the skin
        //TODO NEXT instead check for nodes that would be surrounded by BOUNDARY elements mostly to detect small cuts
        //TODO make parallel and check only nodes of SELECTED elements?
        const double relocation_distance = 1e-5;  //1e-10 OR 1e-3*length
        std::size_t n_relocated_nodes = 0;
        std::size_t n_majority_selected = 0;
        for (auto& rNode : mpModelPart->Nodes()) {
            // Check whether any of the elements surrounding the node is not SELECTED (marked by process as intersected)
            auto& r_neigh_elems = rNode.GetValue(NEIGHBOUR_ELEMENTS);
            std::size_t n_selected = 0;
            std::size_t n_neighbors = 0;
            for (std::size_t i_n = 0; i_n < r_neigh_elems.size(); ++i_n) {
                auto p_neigh = r_neigh_elems(i_n).get();
                // Continue with next node if the current node is part of an element that is not intersected
                n_neighbors++;
                if (p_neigh->Is(SELECTED)) {
                    n_selected++;
                }
            }
            const bool majority_is_selected = n_selected > n_neighbors/2;

            // If there are only intersected elements around the current node, we will check whether the node is intersected and get the intersecting object's normal
            bool relocate_node = false;
            array_1d<double,3> unit_normal;
            //double obj_length = 1e-6;
            if (majority_is_selected) {
                for (std::size_t i_n = 0; i_n < r_neigh_elems.size(); ++i_n) {
                    auto p_neigh = r_neigh_elems(i_n).get();
                    // Get intersected objects of boundary neighbor
                    for (std::size_t i = 0; i < n_elements; ++i) {
                        if (r_elements[i]->Id() == p_neigh->Id()) {
                            // Check whether node is intersected and get the intersected object's normal
                            for (auto& r_int_obj : r_intersected_objects[i]) {
                                array_1d<double, 3> local_coords;
                                // NOTE that r_int_obj would commonly be a line in 2D or a triangle in 3D
                                // NOTE that Line2D2.IsInside as well as Triangle3D3.IsInside allow a distance to the plane of <1.0e-6 * Length()>
                                relocate_node = r_int_obj.GetGeometry().IsInside(rNode.Coordinates(), local_coords);
                                if (relocate_node) {
                                    //obj_length = r_int_obj.GetGeometry().Length();
                                    unit_normal = r_int_obj.GetGeometry().Normal(local_coords);
                                    unit_normal /= norm_2(unit_normal);
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    if (relocate_node) {
                        break;
                    }
                }
                n_majority_selected++;
            }

            // Relocate node (initial position) into the direction of the intersected object's normal
            // TODO NEXT Set node position to initial position?? Move initial position?
            // TODO NEXT Change node positions only for calculation of intersections? -> no, because skin points should be found in intersected elements!
            // TODO NEXT remove crossing of boundary? -> useful if skin points are not exactly on tessellated geometry?
            // NOTE that initial position is mesh output for visualization (ParaView)
            if (relocate_node) {
                rNode.X0() += relocation_distance * unit_normal[0];
                rNode.Y0() += relocation_distance * unit_normal[1];
                rNode.Z0() += relocation_distance * unit_normal[2];
                rNode.X()  += relocation_distance * unit_normal[0];
                rNode.Y()  += relocation_distance * unit_normal[1];
                rNode.Z()  += relocation_distance * unit_normal[2];
                // TODO correction needed for acceleration of the node?? correction needed for FM-ALE??
                n_relocated_nodes++;
            }
            // TODO:
            // 5a. Check detJ of surrounding elements - if negative, then revert?
            // 6a. Live with small change of location of previously embedded overlapping geometries and mesh velocity

            // ALTERNATIVE:
            // 4b. Cluster elements surrounding node into two clusters separated by elements in which skin points are located
            // 5b. If there are more or less than two clusters, then the skin normal will decide whether elements are on the positive or negative side
            //     If there are two clusters a majority vote using the skin normal will decide about which cluster is on the positive and which on the negative side
            // 6b. All elements on the positive side are no BOUNDARY anymore, except if skin points are located inside them
        }
        KRATOS_WATCH(n_majority_selected);

        // Recalculate intersections if nodes where moved
        // NOTE that entirely new process is necessary for taking update node positions into consideration (TODO why??)
        if (n_relocated_nodes > 0) {
            FindIntersectedGeometricalObjectsProcess new_find_intersected_objects_process = FindIntersectedGeometricalObjectsProcess(*mpModelPart, *mpSkinModelPart, find_intersected_options);
            new_find_intersected_objects_process.ExecuteInitialize();
            new_find_intersected_objects_process.FindIntersections();
            std::vector<PointerVector<GeometricalObject>>&  r_new_intersected_objects = new_find_intersected_objects_process.GetIntersections();

            // Mark elements as BOUNDARY that are intersected by the tessellated skin geometry
            #pragma omp parallel for schedule(dynamic)
            for (std::size_t i = 0; i < n_elements; ++i) {
                if (!r_new_intersected_objects[i].empty()) {
                    r_elements[i]->Set(BOUNDARY, true);
                }
            }
        } else {
            // Mark elements as BOUNDARY that are intersected by the tessellated skin geometry
            #pragma omp parallel for schedule(dynamic)
            for (std::size_t i = 0; i < n_elements; ++i) {
                if (!r_intersected_objects[i].empty()) {
                    r_elements[i]->Set(BOUNDARY, true);
                }
            }
        }

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << mSkinModelPartName << " tessellated boundary flags were set. " << n_relocated_nodes << " nodes were relocated." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::LocateSkinPoints()
    {
        // Map skin points (boundary) to elements of volume model part
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        switch (n_dim) {
            case 2:
                MapSkinPointsToElements<2>(mSkinPointsMap);
                break;
            case 3:
                MapSkinPointsToElements<3>(mSkinPointsMap);
                break;
            default:
                KRATOS_ERROR << "Wrong domain size.";
        }

        // Mark elements as BOUNDARY (gamma) in which skin points were located
        // NOTE that BOUNDARY elements will get deactivated and that the split elements BC is applied by means of the extension operators
        //TODO: parallel
        for (const auto& [p_element, _]: mSkinPointsMap) {
            p_element->Set(BOUNDARY, true);
        }

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << mSkinModelPartName << " skin points were mapped to volume mesh elements."  << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetInterfaceFlags()
    {
        // Find the surrogate boundary elements and mark them as INTERFACE (gamma_tilde)
        // Note that we rely on the fact that the neighbors are sorted according to the faces
        // TODO faster to store BOUNDARY element pointers somewhere?
        for (auto& rElement : mpModelPart->Elements()) {
            if (rElement.Is(BOUNDARY)) {
                const std::size_t n_faces = rElement.GetGeometry().FacesNumber();
                auto& r_neigh_elems = rElement.GetValue(NEIGHBOUR_ELEMENTS);
                for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
                    // The neighbour corresponding to the current face is not also BOUNDARY, it means that the current face is surrogate boundary (INTERFACE)
                    // Flag the current neighbour owning the surrogate face as INTERFACE
                    // The nodes will be flagged if required (MLS basis) when creating the cloud
                    auto p_neigh_elem = r_neigh_elems(i_face).get();
                    if (p_neigh_elem != nullptr) {
                        if (!p_neigh_elem->Is(BOUNDARY)) {
                            p_neigh_elem->Set(INTERFACE, true);
                        }
                    }
                }
            }
        }
        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << mSkinModelPartName << " interface flags were set."  << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::DeactivateElementsAndNodes()
    {
        // Deactivate elements in which the (true) boundary is located
        //TODO: make parallel
        for (auto& rElement : mpModelPart->Elements()) {
            if (rElement.Is(BOUNDARY)) {
                rElement.Set(ACTIVE, false);
            }
        }

        // Deactivate nodes that do not belong to any active element (anymore)
        //TODO make parallel and check only nodes of BOUNDARY elements?
        for (auto& rNode : mpModelPart->Nodes()) {
            // Check whether any of the elements surrounding the node is active
            auto& r_neigh_elems = rNode.GetValue(NEIGHBOUR_ELEMENTS);
            bool active_neighbor_found = false;
            for (std::size_t i_n = 0; i_n < r_neigh_elems.size(); ++i_n) {
                auto p_neigh = r_neigh_elems(i_n).get();
                // Continue with next node if the current node is part of an active element
                if (p_neigh->Is(ACTIVE)) {
                    active_neighbor_found = true;
                    break;
                }
            }
            // If there is no active element around the current node, we will deactivate it
            if (!active_neighbor_found) {
                rNode.Set(ACTIVE, false);
            }
        }

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Boundary elements and nodes were deactivated." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::AddSkinIntegrationPointConditions()
    {
        // NOTE that its necessary to call these methods from outside for multiple interface utilities
        /*ResetFlags();
        if (mUseTessellatedBoundary) SetTessellatedBoundaryFlags();
        LocateSkinPoints();
        SetInterfaceFlags();
        DeactivateElementsAndNodes();*/

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
            "mls_extension_operator_order" : 1,
            "active_side_of_skin" : "both",
            "enclosed_area" : "none",
            "cross_boundary_neighbors" : true,
            "use_tessellated_boundary" : false
        })" );

        return default_parameters;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateMeshlessBasedConformingExtensionBasis()
    {
        // Iterate over the split elements to create a vector for each element defining both sides using the element's skin points normals and positions
        // The resulting vector is as long as the number of nodes of the element and a positive value stands for the positive side of the boundary, a negative one for the negative side.
        // Also store the average position and average normal of the skin points located in the element
        SidesVectorToElementsMapType sides_vector_map;
        AverageSkinToElementsMapType avg_skin_map;
        SetSidesVectorsAndSkinNormalsForSplitElements(mSkinPointsMap, sides_vector_map, avg_skin_map);
        //KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Sides vectors and skin normals were set." << std::endl;

        // Iterate over the split elements to create an extension basis for each node of the element (MLS shape functions values for support cloud of node)
        // NOTE that no extension bases will be calculated and added for a node for which not a sufficient number of support nodes were found
        NodesCloudMapType ext_op_map;
        SetExtensionOperatorsForSplitElementNodes(sides_vector_map, avg_skin_map, ext_op_map);
        //KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Extension operators were set." << std::endl;

        // Set the pressure of the first node of an enclosed volume to zero if one side is enclosed.
        auto skin_pt_element_iter = sides_vector_map.begin();
        bool enclosed_pressure_is_set = false;
        //TODO p=0 for all nodes bounding the enclosed area, so no wall conditions are needed
        //while (skin_pt_element_iter != sides_vector_map.end()) {
        if (mPositiveSideIsEnclosed || mNegativeSideIsEnclosed) {
            while (!enclosed_pressure_is_set && skin_pt_element_iter != sides_vector_map.end()) {
                auto p_elem = skin_pt_element_iter->first;
                const array_1d<double, 3> avg_skin_position = avg_skin_map[p_elem].first;
                const array_1d<double, 3> avg_skin_normal = avg_skin_map[p_elem].second;
                enclosed_pressure_is_set = SetEnclosedNodesPressure(*p_elem, skin_pt_element_iter->second, avg_skin_position, avg_skin_normal);
                skin_pt_element_iter++;
            }
        }

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        auto p_element_size_func = GetElementSizeFunction(mpModelPart->ElementsBegin()->GetGeometry());

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Create the interface conditions
        //TODO: THIS CAN BE PARALLEL (WE JUST NEED TO MAKE CRITICAL THE CONDITION ID UPDATE)
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        std::size_t n_skin_pt_conditions_not_added = 0;
        bool successfully_added_pos;
        bool successfully_added_neg;
        for (const auto& [p_element, skin_points_data_vector]: mSkinPointsMap) {
            const auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();

            // For each side of the boundary separately (positive and negative side of gamma), create a pointer vector with all the nodes that affect that side of the current element
            // To be used in the creation of the condition. Positive side refers to adding the positive side's nodes of the element and the negative node's support cloud nodes.
            // NOTE that the obtained clouds are sorted by id to properly get the extension operator data //TODO: necessary? ids are compared anyway?!
            PointerVector<NodeType> cloud_nodes_vector_pos;
            PointerVector<NodeType> cloud_nodes_vector_neg;
            CreateCloudNodeVectorsForSplitElement(*p_element, sides_vector_map[p_element], ext_op_map, cloud_nodes_vector_pos, cloud_nodes_vector_neg);

            // Calculate parent element size for the SBM BC imposition
            const double h = p_element_size_func(r_geom);

            // Get vector defining positive and negative side of the boundary
            auto& sides_vector = sides_vector_map[p_element];

            // Iterate over the element's skin points adding a positive side and a negative side condition for each skin point
            for (std::size_t i_skin_pt = 0; i_skin_pt < skin_points_data_vector.size(); ++i_skin_pt) {
                // Get the skin point's position and area normal (integration point of the boundary)
                auto skin_pt_data = skin_points_data_vector[i_skin_pt];
                const array_1d<double, 3> skin_pt_position = std::get<0>(skin_pt_data);
                const array_1d<double, 3> skin_pt_area_normal = std::get<1>(skin_pt_data);

                // Get the split element's shape function values and derivatives at the skin/ integration point
                Vector skin_pt_N(n_nodes);
                Matrix skin_pt_DNDX = ZeroMatrix(n_nodes, n_dim);
                GetDataForSplitElementIntegrationPoint(*p_element, skin_pt_position, skin_pt_N, skin_pt_DNDX);

                // Add skin pt. condition for positive side of boundary - using support cloud data for negative nodes
                // NOTE that the boundary normal is negative in order to point outwards (from positive to negative side),
                // because positive side is where dot product of vector to node with average normal is positive
                successfully_added_pos = AddIntegrationPointCondition(*p_element, sides_vector, h, skin_pt_position, -skin_pt_area_normal,
                ext_op_map, cloud_nodes_vector_pos, skin_pt_N, skin_pt_DNDX, ++max_cond_id, /*ConsiderPositiveSide=*/true);

                // Add skin pt. condition for negative side of boundary - using support cloud data for positive nodes
                // NOTE that boundary normal is opposite (from negative to positive side)
                successfully_added_neg = AddIntegrationPointCondition(*p_element, sides_vector, h, skin_pt_position, skin_pt_area_normal,
                ext_op_map, cloud_nodes_vector_neg, skin_pt_N, skin_pt_DNDX, ++max_cond_id, /*ConsiderPositiveSide=*/false);

                if (!successfully_added_pos || !successfully_added_neg) {
                    n_skin_pt_conditions_not_added++;
                }
            }
        }
        if (n_skin_pt_conditions_not_added > 0) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "Integration point conditions were NOT successfully added for at least one side of "
                << n_skin_pt_conditions_not_added << " skin points." << std::endl;
        }

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << mSkinModelPartName << " skin point conditions were created." << std::endl;
    }

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedInterfaceUtility::MapSkinPointsToElements(
        SkinPointsToElementsMapType& rSkinPointsMap)
    {
        // Check that the skin model part has elements
        KRATOS_ERROR_IF_NOT(mpSkinModelPart->NumberOfElements())
            << "There are no elements in skin model part (boundary) '" << mpSkinModelPart->FullName() << "'." << std::endl;

        // Set the bin-based fast point locator utility
        const std::size_t point_locator_max_results = 10000;
        const double point_locator_tolerance = 1.0e-5;
        BinBasedFastPointLocator<TDim> point_locator(*mpModelPart);
        point_locator.UpdateSearchDatabase();
        typename BinBasedFastPointLocator<TDim>::ResultContainerType search_results(point_locator_max_results);

        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;

        //const std::size_t n_gp_per_element = mpSkinModelPart->ElementsBegin()->GetGeometry().IntegrationPointsNumber(integration_method);
        //const std::size_t n_skin_points = mpSkinModelPart->NumberOfElements() * n_gp_per_element;
        std::vector<Element::Pointer> skin_point_located_elements;
        std::vector<array_1d<double, 3>> skin_point_positions;
        std::vector<array_1d<double, 3>> skin_point_normals;

        // Search the skin points (skin model part integration points) in the volume mesh elements
        //TODO make parallel and push back into vectors unique?
        std::size_t n_skin_points_not_found = 0;
        for (Element& rSkinElement : mpSkinModelPart->Elements()) {
            const auto& r_skin_geom = rSkinElement.GetGeometry();
            const std::size_t n_gp = r_skin_geom.IntegrationPointsNumber(integration_method);
            const GeometryType::IntegrationPointsArrayType& integration_points = r_skin_geom.IntegrationPoints(integration_method);
            // Get detJ for all integration points of the skin element
            Vector integration_point_jacobians;
            r_skin_geom.DeterminantOfJacobian(integration_point_jacobians, integration_method);

            for (std::size_t i_gp = 0; i_gp < n_gp; ++i_gp) {
                // Get position of skin point
                const array_1d<double, 3> skin_pt_local_coords = integration_points[i_gp].Coordinates();
                array_1d<double, 3> skin_pt_position = ZeroVector(3);
                r_skin_geom.GlobalCoordinates(skin_pt_position, skin_pt_local_coords);

                // Get normal at the skin point and make its length a measure of the area/ integration weight
                array_1d<double,3> skin_pt_area_normal = r_skin_geom.Normal(skin_pt_local_coords);
                // Normalize normal
                skin_pt_area_normal /= std::max(norm_2(skin_pt_area_normal), 1e-10);  // tolerance = std::pow(1e-3 * h, Dim-1)
                // Scale normal with integration weight
                const double integration_weight = integration_point_jacobians(i_gp) * integration_points[i_gp].Weight();
                skin_pt_area_normal *= integration_weight;

                // Search for the skin point in the volume mesh to get the element containing the point
                Vector aux_N(TDim+1);
                Element::Pointer p_element = nullptr;
                const bool is_found = point_locator.FindPointOnMesh(
                    skin_pt_position, aux_N, p_element,
                    search_results.begin(), point_locator_max_results, point_locator_tolerance);

                // Add data to vectors
                if (is_found){
                    skin_point_located_elements.push_back(p_element);
                    skin_point_positions.push_back(skin_pt_position);
                    skin_point_normals.push_back(skin_pt_area_normal);
                } else {
                    n_skin_points_not_found++;
                    //KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility")
                    //    << "A skin point has not been found in any volume model part element. Point coordinates: ("
                    //    << skin_pt_position[0] << " , " << skin_pt_position[1] << " , " << skin_pt_position[2] << ")" << std::endl;
                }
            }
        }
        if (n_skin_points_not_found > 0) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility")
                << n_skin_points_not_found << " skin points have not been found in any volume model part element." << std::endl;
        }

        for (Element& rElement : mpModelPart->Elements()) {
            // Find the indices of all skin points that are located in the volume mesh element
            std::vector<std::size_t> skin_point_indices;
            for (std::size_t i_skin_pt = 0; i_skin_pt < skin_point_located_elements.size(); ++i_skin_pt) {
                Element::Pointer p_elem = skin_point_located_elements[i_skin_pt];
                if (p_elem->Id() == rElement.Id()) {
                    skin_point_indices.push_back(i_skin_pt);
                }
            }
            // Add the position and area normal of all skin points located inside the element to the data vector of the element
            std::size_t n_skin_pts_in_element = skin_point_indices.size();
            if (n_skin_pts_in_element > 0) {
                SkinPointsDataVectorType skin_points_data_vector(n_skin_pts_in_element);
                for (std::size_t i_pt = 0; i_pt < n_skin_pts_in_element; ++i_pt) {
                    const std::size_t skin_pt_index = skin_point_indices[i_pt];
                    skin_points_data_vector[i_pt] = std::make_pair(skin_point_positions[skin_pt_index], skin_point_normals[skin_pt_index]);
                }
                // Add data vector of the element to the skin points map
                auto skin_points_key_data = std::make_pair(&rElement, skin_points_data_vector);
                rSkinPointsMap.insert(skin_points_key_data);
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetSidesVectorsAndSkinNormalsForSplitElements(
        const SkinPointsToElementsMapType& rSkinPointsMap,
        SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap)
    {
        // TODO can be parallel if insert into rSidesVectorMap is unique?
        // TODO NEXT decide on sides based on all skin points individually using majority vote?
        for (const auto& [p_element, skin_points_data_vector]: rSkinPointsMap) {
            // Calculate average position and normal of the element's skin points
            std::size_t n_skin_points = 0;
            array_1d<double, 3> avg_position(3, 0.0);
            array_1d<double, 3> avg_normal(3, 0.0);
            for (const auto& skin_pt_data: skin_points_data_vector) {
                n_skin_points++;
                avg_position += std::get<0>(skin_pt_data);
                avg_normal += std::get<1>(skin_pt_data);
            }
            avg_normal /= norm_2(avg_normal);
            avg_position /= n_skin_points;

            // Store the average normal of the skin points located in the element
            auto avg_position_and_normal = std::make_pair(avg_position, avg_normal);
            rAvgSkinMap.insert(std::make_pair(p_element, avg_position_and_normal));

            // Compute the dot product for each node between a vector from the average skin position to the node and the skin's average normal
            // Note that for a positive dot product the node is saved as being on the positive side of the boundary, negative dot product = negative side
            auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            Vector sides_vector(n_nodes);
            int majority_side = 0;
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                const auto& r_node = r_geom[i_node];
                array_1d<double, 3> skin_pt_to_node = r_node.Coordinates() - avg_position;
                const double dot_product = inner_prod(skin_pt_to_node, avg_normal);
                if (dot_product > 0.0) {
                    sides_vector[i_node] =  1.0;
                    majority_side++;
                } else if (dot_product < 0.0) {
                    sides_vector[i_node] = -1.0;
                    majority_side--;
                } else {
                    sides_vector[i_node] = 0.0;
                }
            }

            // If one of the values is zero, then the (average) will go through that node, so we put the node on the other side of where the majority of the element's nodes lie
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                if (sides_vector[i_node] == 0.0) {
                    if (majority_side > 0) {
                        sides_vector[i_node] = -1.0;
                    } else{
                        sides_vector[i_node] = 1.0;
                    }
                }
            }

            // Store the vector deciding on the positive and negative side of the element's nodes
            rSidesVectorMap.insert(std::make_pair(p_element, sides_vector));
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetExtensionOperatorsForSplitElementNodes(
        const SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap,
        NodesCloudMapType& rExtensionOperatorMap)
    {
        // Get the MLS shape functions function
        auto p_meshless_sh_func = GetMLSShapeFunctionsFunction();

        // Get support nodes for both sides of all split elements and calculate extension operators for all nodes of split elements
        // NOTE that only extension operators are calculated and added to the map if a sufficient number of support nodes was found
        for (const auto& [p_element, sides_vector]: rSidesVectorMap) {
            const auto r_geom = p_element->GetGeometry();

            // Get averaged position and normal of the skin points located inside the element
            auto avg_position_and_normal = rAvgSkinMap[p_element];
            const array_1d<double, 3> avg_position = std::get<0>(avg_position_and_normal);
            const array_1d<double, 3> avg_normal = std::get<1>(avg_position_and_normal);

            // Check for which nodes and sides the extension operator has not been calculated (successfully) yet
            const std::size_t n_nodes = r_geom.PointsNumber();
            Vector ex_op_was_calculated(n_nodes);
            bool calculatePositiveSupport = false;
            bool calculateNegativeSupport = false;
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                const auto p_node = r_geom(i_node);
                const std::size_t found_in_map = rExtensionOperatorMap.count(p_node);
                if (found_in_map) {
                    ex_op_was_calculated[i_node] =  1.0;
                // Extension operator for node still needs to be calculated
                } else {
                    ex_op_was_calculated[i_node] =  0.0;
                    // If node is on the positive side, then its extension operator needs support nodes on the negative side
                    if (sides_vector[i_node] > 0.0) {
                        calculateNegativeSupport = true;
                    // If node is on the negative side, then its extension operator needs support nodes on the positive side
                    } else {
                        calculatePositiveSupport = true;
                    }
                }
            }

            // Calculate support for the positive side using support nodes from the positive side
            if (calculatePositiveSupport) {
                // Initialize the storage for the support/ cloud nodes and their coordinates
                Matrix cloud_nodes_coordinates;
                PointerVector<NodeType> cloud_nodes;

                // Get the element's positive nodes as search basis and support nodes on the positive side
                std::vector<NodeType::Pointer> elem_pos_nodes;
                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                    if (sides_vector[i_node] > 0.0) {
                        elem_pos_nodes.push_back(r_geom(i_node));
                    }
                }
                SetLateralSupportCloud(elem_pos_nodes, avg_position, avg_normal, cloud_nodes, cloud_nodes_coordinates, /*ConsiderPositiveSide=*/true);

                // Continue if the number of support nodes is sufficient for the calculation of the extension operator
                // NOTE that an extension operator might be calculated successfully for the respective node(s) with another element as basis
                const std::size_t n_cloud_nodes = cloud_nodes.size();
                if (n_cloud_nodes >= GetRequiredNumberOfPoints()) {

                    // Add the extension operator for all nodes on the negative side of the element for which is has not been calculated yet
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        if (ex_op_was_calculated[i_node] < 1.0 && sides_vector[i_node] < 0.0) {
                            const auto p_node = r_geom(i_node);

                            // Calculate the extension basis in the current node (MLS shape functions)
                            Vector N_container;
                            const array_1d<double,3> r_coords = p_node->Coordinates();
                            const double kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
                            p_meshless_sh_func(cloud_nodes_coordinates, r_coords, kernel_rad, N_container);

                            // Save the extension operator nodal data to the extension operator map
                            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
                            for (std::size_t i_cl_nod = 0; i_cl_nod < n_cloud_nodes; ++i_cl_nod) {
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

            // Calculate support for the negative side using support nodes from the negative side
            if (calculateNegativeSupport) {
                // Initialize the storage for the support/ cloud nodes and their coordinates
                Matrix cloud_nodes_coordinates;
                PointerVector<NodeType> cloud_nodes;

                // Get the element's negative nodes as search basis and support nodes on the negative side
                std::vector<NodeType::Pointer> elem_neg_nodes;
                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                    if (sides_vector[i_node] < 0.0) {
                        elem_neg_nodes.push_back(r_geom(i_node));
                    }
                }
                SetLateralSupportCloud(elem_neg_nodes, avg_position, -avg_normal, cloud_nodes, cloud_nodes_coordinates, /*ConsiderPositiveSide=*/false);

                // Continue if the number of support nodes is sufficient for the calculation of the extension operator
                // NOTE that an extension operator might be calculated successfully for the respective node(s) with another element as basis
                const std::size_t n_cloud_nodes = cloud_nodes.size();
                if (n_cloud_nodes >= GetRequiredNumberOfPoints()) {

                    // Add the extension operator for all nodes on the positive side of the element for which is has not been calculated yet
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        if (ex_op_was_calculated[i_node] < 1.0 && sides_vector[i_node] > 0.0) {
                            const auto p_node = r_geom(i_node);

                            // Calculate the extension basis in the current node (MLS shape functions)
                            Vector N_container;
                            const array_1d<double,3> r_coords = p_node->Coordinates();
                            const double kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
                            p_meshless_sh_func(cloud_nodes_coordinates, r_coords, kernel_rad, N_container);

                            // Save the extension operator nodal data to the extension operator map
                            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
                            for (std::size_t i_cl_nod = 0; i_cl_nod < n_cloud_nodes; ++i_cl_nod) {
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
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetLateralSupportCloud(
        const std::vector<NodeType::Pointer>& rSameSideNodes,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates,
        const bool ConsiderPositiveSide)
    {
        // Find the same side support cloud of nodes
        // NOTE that we use an unordered_set to ensure that these are unique
        // NOTE that we might add deactivated nodes to the set here to facilitate the search
        NodesCloudSetType aux_set;
        std::vector<NodeType::Pointer> cur_layer_nodes;
        std::vector<NodeType::Pointer> prev_layer_nodes;
        const std::size_t n_layers = mMLSExtensionOperatorOrder + 1;

        // Check starting points
        const std::size_t n_same_side_nodes = rSameSideNodes.size();
        KRATOS_ERROR_IF(n_same_side_nodes == 0) << "No starting nodes given for computing the lateral support cloud." << std::endl;

        // Add first given nodes to cloud nodes set (1st layer of sampling points)
        for (std::size_t i_node = 0; i_node < n_same_side_nodes; ++i_node) {
            const auto p_node = rSameSideNodes[i_node];
            // Calculate dot product of average skin normal of the element and the normalized vector between averaged skin point and the element's node
            array_1d<double, 3> avg_skin_pt_to_node = p_node->Coordinates() - rAvgSkinPosition;
            avg_skin_pt_to_node /= norm_2(avg_skin_pt_to_node);
            const double dot_product = inner_prod(avg_skin_pt_to_node, rAvgSkinNormal);
            if (dot_product > 0.1) {
                aux_set.insert(p_node);
            }
            prev_layer_nodes.push_back(p_node);
            // Check if same side node is active, only then the condition will be created
            if (!p_node->Is(ACTIVE)) {
                return;
            }
        }

        // Add second layer of nodal neighbors. These are the nodal neighbors of the given nodes in direction of the average skin normal
        // NOTE that if mCrossBoundaryNeighbors is true, then nodes can be taken from BOUNDARY elements (direct neighbors) if the dot product is sufficiently high.
        if (mCrossBoundaryNeighbors) {
            AddLateralSupportLayerCrossingBoundary(rAvgSkinPosition, rAvgSkinNormal, prev_layer_nodes, cur_layer_nodes, aux_set);
        } else {
            AddLateralSupportLayer(rAvgSkinPosition, rAvgSkinNormal, prev_layer_nodes, cur_layer_nodes, aux_set);
        }
        prev_layer_nodes = cur_layer_nodes;
        cur_layer_nodes.clear();

        // Add more layers of nodal neighbors of the current nodes to the cloud of nodes
        // NOTE that we check the order of the MLS interpolation to add nodes from enough interior
        // NOTE that we start from 2 here as the first and second layer already have been added (so for 2D linear operators there is no iteration)
        for (std::size_t i_layer = 2; i_layer < n_layers; ++i_layer) {
            AddLateralSupportLayer(rAvgSkinPosition, rAvgSkinNormal, prev_layer_nodes, cur_layer_nodes, aux_set);
            prev_layer_nodes = cur_layer_nodes;
            cur_layer_nodes.clear();
        }

        // Check how many of the current support nodes are active
        std::size_t n_cloud_nodes = 0;
        for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
            if ((*it_set)->Is(ACTIVE)) {
                n_cloud_nodes++;
            }
        }

        // If there are not enough active support nodes to perform the MLS calculation add another layer of neighboring nodes
        // Add maximal three extra layers
        std::size_t n_extra_layers = 0;
        while (n_cloud_nodes < GetRequiredNumberOfPoints()+1 && n_extra_layers < 3) {
            AddLateralSupportLayer(rAvgSkinPosition, rAvgSkinNormal, prev_layer_nodes, cur_layer_nodes, aux_set);
            n_extra_layers++;
            n_cloud_nodes = 0;
            for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
                if ((*it_set)->Is(ACTIVE)) {
                    n_cloud_nodes++;
                }
            }
            prev_layer_nodes = cur_layer_nodes;
            cur_layer_nodes.clear();
        }
        if (n_extra_layers > 1) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << n_extra_layers << " extra layers of points needed for MLS calculation." << std::endl;
        }

        // Add obtained cloud nodes to the cloud node vector if they are active and sort them by id
        //TODO sorting really necessary or helpful??
        rCloudNodes.resize(n_cloud_nodes);
        std::size_t aux_i = 0;
        for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
            if ((*it_set)->Is(ACTIVE)) {
                rCloudNodes(aux_i++) = *it_set;
                // For visualization: Mark active support nodes as BOUNDARY if they are on the positive side (support for negative side nodes) and INTERFACE for the negative side (support for positive side nodes)
                if (ConsiderPositiveSide) {
                    (*it_set)->Set(BOUNDARY, true);
                } else {
                    (*it_set)->Set(INTERFACE, true);
                }
            }
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

    void ShiftedBoundaryPointBasedInterfaceUtility::AddLateralSupportLayer(
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesCloudSetType& SupportNodesSet)
    {
        // Find elemental neighbors of the nodes of the previous layer
        for (auto& p_node : PreviousLayerNodes) {
            auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);

            // Add all nodes of neighboring elements to cloud nodes set if element is not (!) BOUNDARY and in normal direction
            // This way the boundary cannot be crossed (not 'ACTIVE' might be used instead here?)
            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                if (p_elem_neigh != nullptr && !p_elem_neigh->Is(BOUNDARY)) {
                    const auto& r_geom = p_elem_neigh->GetGeometry();
                    for (std::size_t i_neigh_node = 0; i_neigh_node < r_geom.PointsNumber(); ++i_neigh_node) {
                        // Add node of neighboring element only if it is in the inwards normal direction of the element's averaged skin points
                        // NOTE this is done for a more robust separation of both sides of the boundary
                        NodeType::Pointer p_neigh = r_geom(i_neigh_node);
                        // Calculate dot product of average skin normal of the element and the normalized vector between averaged skin point and the element's node
                        array_1d<double, 3> avg_skin_pt_to_node = p_neigh->Coordinates() - rAvgSkinPosition;
                        avg_skin_pt_to_node /= norm_2(avg_skin_pt_to_node);
                        const double dot_product = inner_prod(avg_skin_pt_to_node, rAvgSkinNormal);
                        if (dot_product > 0.1) {
                            auto set_return = SupportNodesSet.insert(p_neigh);
                            // If the node was inserted into the set as a new element, then add it to the current layer
                            if (set_return.second) {
                                CurrentLayerNodes.push_back(p_neigh);
                            }
                        }
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::AddLateralSupportLayerCrossingBoundary(
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesCloudSetType& SupportNodesSet)
    {
        // Find elemental neighbors of the nodes of the previous layer
        for (auto& p_node : PreviousLayerNodes) {
            auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);

            // Add all nodes of neighboring elements to cloud nodes set if node is in normal direction
            // NOTE that nodes of BOUNDARY elements can be taken as well here if the dot product and distance are sufficiently high.
            // This way the boundary can be crossed.
            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                if (p_elem_neigh != nullptr) {
                    const auto& r_geom = p_elem_neigh->GetGeometry();
                    for (std::size_t i_neigh_node = 0; i_neigh_node < r_geom.PointsNumber(); ++i_neigh_node) {
                        // Add node of neighboring element only if it is in the inwards normal direction of the element's averaged skin points
                        // NOTE this is done for a separation of both sides of the boundary here
                        NodeType::Pointer p_neigh = r_geom(i_neigh_node);
                        // Calculate the distance between average skin point and node and normalize the vector
                        array_1d<double, 3> avg_skin_pt_to_node = p_neigh->Coordinates() - rAvgSkinPosition;
                        const double d = norm_2(avg_skin_pt_to_node);
                        avg_skin_pt_to_node /= d;
                        const double dot_product = inner_prod(avg_skin_pt_to_node, rAvgSkinNormal);
                        const double h = GetElementSizeFunction(r_geom)(r_geom);
                        if ( (!p_elem_neigh->Is(BOUNDARY) && dot_product > 0.1) || (dot_product > 0.01 && d > 0.01*h) ) {
                            auto set_return = SupportNodesSet.insert(p_neigh);
                            // Only add nodes as basis for the next layer if they have not been added to the set already (because then they would be considered again and again)
                            if (set_return.second) {
                                CurrentLayerNodes.push_back(p_neigh);
                            }
                        }
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::CreateCloudNodeVectorsForSplitElement(
        const ElementType& rElement,
        const Vector& rSidesVector,
        NodesCloudMapType& rExtensionOperatorMap,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide)
    {
        // Create an auxiliary set with all the cloud nodes that affect the current element for each side separately
        // NOTE that a node can only be found if sufficient cloud nodes were found for the creation of the extension basis
        // NOTE that only active nodes are part of the node cloud
        NodesCloudSetType cloud_nodes_set_pos;
        NodesCloudSetType cloud_nodes_set_neg;
        const auto& r_geom = rElement.GetGeometry();
        for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            NodeType::Pointer p_node = r_geom(i_node);
            if (rSidesVector(i_node) > 0.0) {
                // Add positive side node to cloud nodes set of positive side of the boundary
                // NOTE they might not be part of the negative node's support because they are too close to the other side or not active
                cloud_nodes_set_pos.insert(p_node);
                // Add positive side's node's cloud nodes to cloud nodes set of negative side of the boundary
                const std::size_t found = rExtensionOperatorMap.count(p_node);
                if (found) {
                    auto& r_ext_op_data = rExtensionOperatorMap[p_node];
                    for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                        auto& p_node = std::get<0>(*it_data);
                        cloud_nodes_set_neg.insert(p_node);
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
                        auto& p_node = std::get<0>(*it_data);
                        cloud_nodes_set_pos.insert(p_node);
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

    void ShiftedBoundaryPointBasedInterfaceUtility::GetDataForSplitElementIntegrationPoint(
        const ElementType& rElement,
        const array_1d<double,3>& rIntPtCoordinates,
        Vector& rIntPtShapeFunctionValues,
        Matrix& rIntPtShapeFunctionDerivatives)
    {
        const auto& r_geom = rElement.GetGeometry();

        // Compute the local coordinates of the integration point in the element's geometry
        //Geometry<Node>::CoordinatesArrayType aux_global_coords = ZeroVector(3);
        //r_bd_geom.GlobalCoordinates(aux_global_coords, r_integration_points[i_int_pt].Coordinates());
        array_1d<double, 3> int_pt_local_coords = ZeroVector(3);
        r_geom.PointLocalCoordinates(int_pt_local_coords, rIntPtCoordinates);

        // Get N of the element at the integration point
        r_geom.ShapeFunctionsValues(rIntPtShapeFunctionValues, int_pt_local_coords);

        // Get DN_DX of the element at the integration point
        Matrix aux_DN_DXi_parent, aux_J_parent, aux_J_inv_parent;
        double aux_detJ_parent;
        r_geom.ShapeFunctionsLocalGradients(aux_DN_DXi_parent, int_pt_local_coords);
        r_geom.Jacobian(aux_J_parent, int_pt_local_coords);
        MathUtils<double>::InvertMatrix(aux_J_parent, aux_J_inv_parent, aux_detJ_parent);
        rIntPtShapeFunctionDerivatives = prod(aux_DN_DXi_parent, aux_J_inv_parent);
    }

    bool ShiftedBoundaryPointBasedInterfaceUtility::AddIntegrationPointCondition(
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
            const auto p_node = r_geom(i_node);
            // If node is on the side that is being considered, then add the standard shape function contribution of the node at the position of the skin point
            if (ConsiderPositiveSide != (rSidesVector[i_node] <= 0.0)) {
                // If a node on the side that is being considered is not active, then no wall condition is created
                //TODO NEXT redistribute shape functions? Or calculate extension basis for skin point directly??
                if (!p_node->Is(ACTIVE)) {
                    KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "No wall condition will be created for one side of the skin point because Node No." << p_node->Id() << " is not active." << std::endl;
                    return false;
                }
                // Note that we need to check for the ids to match as we do not know the node's position in the node vector
                for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                    auto& p_cl_node = rCloudNodeVector(i_cl);
                    if (p_node->Id() == p_cl_node->Id()) {
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

                // If node on the other side does not have an extension basis, then no wall condition is created
                const std::size_t found = rExtensionOperatorMap.count(p_node);
                if (!found) {
                    KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "No wall condition will be created for one side of the skin point because no extension operator was available for Node No." << p_node->Id() << std::endl;
                    return false;
                }

                // Get the node's extension operator data
                auto& ext_op_data = rExtensionOperatorMap[p_node];

                // Iterate over the node's extension operator data and apply the cloud node weight (i_cl_node_N) to make the basis conformant
                // Note that we need to check for the ids to match as we do not know the node's position in the node vector
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
        const double skin_pt_weight = norm_2(rIntPtAreaNormal);
        p_cond->SetValue(NORMAL, rIntPtAreaNormal/ skin_pt_weight);
        p_cond->SetValue(INTEGRATION_WEIGHT, skin_pt_weight);
        p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
        p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);

        return true;
    }

    bool ShiftedBoundaryPointBasedInterfaceUtility::SetEnclosedNodesPressure(
        ElementType& rElement,
        const Vector& rSidesVector,
        const array_1d<double, 3>& rAvgSkinPosition,
        const array_1d<double, 3>& rAvgSkinNormal)
    {
        auto& r_geom = rElement.GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            if ((mPositiveSideIsEnclosed && rSidesVector[i_node] > 0) or (mNegativeSideIsEnclosed && rSidesVector[i_node] < 0)) {

                auto& r_node = r_geom[i_node];
                // Calculate dot product of average skin normal of the element and the normalized vector between averaged skin point and the element's node
                array_1d<double, 3> avg_skin_pt_to_node = r_node.Coordinates() - rAvgSkinPosition;
                avg_skin_pt_to_node /= norm_2(avg_skin_pt_to_node);
                double dot_product = inner_prod(avg_skin_pt_to_node, rAvgSkinNormal);
                if (mNegativeSideIsEnclosed) {
                    dot_product *= -1.0;
                }
                if (dot_product > 0.1 && r_node.Is(ACTIVE)) {
                    r_node.Fix(PRESSURE);
                    r_node.FastGetSolutionStepValue(PRESSURE) = 0.0;
                    return true;
                }
            }
        }
        return false;
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
