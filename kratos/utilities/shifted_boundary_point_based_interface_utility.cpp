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
#include "includes/exception.h"
#include "includes/kratos_flags.h"
#include "includes/lock_object.h"
#include "includes/smart_pointers.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
#include "input_output/logger.h"
#include "shifted_boundary_point_based_interface_utility.h"
#include "utilities/element_size_calculator.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/shifted_boundary_point_based_interface_utility.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "geometries/plane_3d.h"
#include "processes/find_intersected_geometrical_objects_process.h"


namespace Kratos
{
    namespace ShiftedBoundaryUtilityInternals {

    template <>
    double CalculatePointDistance<2>(
        const Geometry<Node>& rObjectGeometry,
        const Point& rPoint)
    {
        return GeometryUtils::PointDistanceToLineSegment3D(
            rObjectGeometry[0],
            rObjectGeometry[1],
            rPoint);
    }

    template <>
    double CalculatePointDistance<3>(
        const Geometry<Node>& rObjectGeometry,
        const Point& rPoint)
    {
        return GeometryUtils::PointDistanceToTriangle3D(
            rObjectGeometry[0],
            rObjectGeometry[1],
            rObjectGeometry[2],
            rPoint);
    }

    template <>
    Plane3D CreateIntersectionPlane<2>(const std::vector< array_1d<double,3> >& rIntPtsVector)
    {
        // Since the Plane3D object only works in 3D, in 2D we define the intersection plane
        // by extruding the intersection point 0 in the z-direction.
        array_1d<double,3> z_coord_pt = rIntPtsVector[0];
        z_coord_pt[2] = 1.0;
        return Plane3D(Point{rIntPtsVector[0]}, Point{rIntPtsVector[1]}, Point{z_coord_pt});
    }

    template <>
    Plane3D CreateIntersectionPlane<3>(const std::vector< array_1d<double,3> >& rIntPtsVector)
    {
        return Plane3D(Point{rIntPtsVector[0]}, Point{rIntPtsVector[1]}, Point{rIntPtsVector[2]});
    }

    // see FluidElementUtilities
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

    // see FluidElementUtilities
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

    // see FluidElementUtilities
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

    }  // namespace ShiftedBoundaryUtilityInternals

    const Parameters ShiftedBoundaryPointBasedInterfaceUtility::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "skin_model_part_name" : "",
            "boundary_sub_model_part_name" : "",
            "boundary_wall_condition_name" : "",
            "extension_operator_type" : "MLS",
            "mls_extension_operator_order" : 1,
            "active_side_of_skin" : "both",
            "enclosed_area" : "none",
            "use_tessellated_boundary" : true
        })" );

        return default_parameters;
    }

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
        // Get model part for skin points
        const std::string skin_points_model_part_name = mSkinModelPartName + "Points";
        mpSkinPointsModelPart = &rModel.GetModelPart(skin_points_model_part_name);
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", mpSkinPointsModelPart->NumberOfNodes() != 0) << "Provided SBM model part has nodes." << std::endl;
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", mpSkinPointsModelPart->NumberOfElements() != 0) << "Provided SBM model part has elements." << std::endl;
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", mpSkinPointsModelPart->NumberOfConditions() != 0) << "Provided SBM model part has conditions." << std::endl;

        // Set the order of the MLS extension operator used in the MLS shape functions utility
        mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();
        // If true, the basis is created such that the surrogate boundary gradient is kept
        const std::string ext_op_type = ThisParameters["extension_operator_type"].GetString();
        if (ext_op_type == "MLS") {
            mExtensionOperator = ExtensionOperator::MLS;
        } else if (ext_op_type == "RBF") {
            mExtensionOperator = ExtensionOperator::RBF;
        } else {
            KRATOS_ERROR << "Wrong 'extension_operator_type' provided. Only 'MLS' 'extension_operator_type' is supported by point based shifted boundary interface utility." << std::endl;
        }

        // Set the shifted-boundary condition prototype to be used in the condition creation
        const std::string wall_condition_name = ThisParameters["boundary_wall_condition_name"].GetString();
        KRATOS_ERROR_IF(wall_condition_name == "") << "SBM boundary wall condition has not been provided." << std::endl;
        mpConditionPrototype = &KratosComponents<Condition>::Get(wall_condition_name);

        // Set which side of the skin model part is considered
        const std::string active_side_of_skin = ThisParameters["active_side_of_skin"].GetString();
        if (active_side_of_skin != "positive" && active_side_of_skin != "negative" && active_side_of_skin != "both") {
            KRATOS_ERROR << "Unknown 'active_side_of_skin' given. 'positive', 'negative' or 'both' sides are supported by point based shifted boundary interface utility." << std::endl;
        }
        //TODO actually deactivate side that is not considered --> flood fill and majority vote like Manuel
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", active_side_of_skin != "both") << "So far always both sides of the skin are being considered." << std::endl;
        //--- actually being used?? - No extension operators will be calculated and no wall conditions will be created for inactive side
        //NOTE that this does not have any effect right now
        if (active_side_of_skin == "positive") {
            mNegativeSideIsActive = false;
        } else if (active_side_of_skin == "negative") {
            mPositiveSideIsActive = false;
        }

        // Set which side of the skin model part is an enclosed area
        //TODO use flood fill process (--> Manuel) and Octree to check and search for enclosed areas for setting the pressure! (robustness)
        const std::string enclosed_area = ThisParameters["enclosed_area"].GetString();
        if (enclosed_area != "positive" && enclosed_area != "negative" && enclosed_area != "none") {
            KRATOS_ERROR << "Unknown 'enclosed_area' keyword given. 'positive' or 'negative' side or 'none' are supported by point based shifted boundary interface utility." << std::endl;
        }
        if (active_side_of_skin == "both") {
            if (enclosed_area == "positive") {
                mPositiveSideIsEnclosed = true;
            } else if (enclosed_area == "negative") {
                mNegativeSideIsEnclosed = true;
            }
        } else if (enclosed_area != "none") {
            KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "The keyword 'enclosed_area' will be ignored, because one side of the skin is not active as given in the 'active_side_of_skin' keyword." << std::endl;
        }

        // If true, then elements will be declared SBM_BOUNDARY which are intersected by the tessellated skin geometry
        mUseTessellatedBoundary = ThisParameters["use_tessellated_boundary"].GetBool();
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateAndAddPointBasedInterface(
        const bool DeactivateUnstableClusters
    )
    {
        // NOTE that its necessary to call these methods from outside for multiple interface utilities

        // Reset all SBM flags
        ResetFlags();
        // Use tessellated boundary description of skin model part to find SBM_BOUNDARY element in which no skin points are located
        // NOTE that nodes which are very close to the skin might be relocated here
        if (mUseTessellatedBoundary) SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes();
        // Locate skin points AFTER possible node relocation of SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes
        LocateSkinPoints();
        // Set the SBM_INTERFACE flag AFTER setting SBM_BOUNDARY flags is completed
        SetInterfaceFlags();
        // Deactivate SBM_BOUNDARY elements and nodes which are surrounded by deactivated elements
        DeactivateElementsAndNodes(DeactivateUnstableClusters);
        // Calculate the extension operators for nodes of elements in which skin points are located
        // and add conditions for one or both sides of all skin points
        CalculateAndAddSkinIntegrationPointConditions();

        //if (mFindEnclosedVolumes) FixEnclosedVolumesPressure();
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::ResetFlags()
    {
        // Activate all elements and initialize flags to false
        // NOTE Resetting the SBM flags will eliminate previously embedded model parts except for their wall conditions!
        // TODO container for relocated nodes? so that they can be set back to their original position if needed
        block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
            rNode.Set(ACTIVE, true);             // Nodes that belong to the elements to be assembled
            rNode.Set(SBM_BOUNDARY, false);      // Nodes that belong to the support clouds of the positive side  //TODO use variables with correct names?
            rNode.Set(SBM_INTERFACE, false);     // Nodes that belong to the support clouds of the negative side  //TODO use variables with correct names?
            rNode.Set(MODIFIED, false);          // Nodes that have been relocated because of a small distance to the skin geometry
        });
        block_for_each(mpModelPart->Elements(), [](ElementType& rElement){
            rElement.Set(ACTIVE, true);      // Elements in the positive distance region (the ones to be assembled)
            rElement.Set(SBM_BOUNDARY, false);   // Elements in which the skin geometry is located (true boundary gamma)
            rElement.Set(SBM_INTERFACE, false);  // Positive distance elements owning the surrogate boundary nodes
        });

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Boundary and interface flags were reset and all elements and nodes re-activated." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes()
    {
        // Reset SELECTED flag which is used to mark elements which are intersected by the skin geometry
        VariableUtils().SetFlag(SELECTED, false, mpModelPart->Elements());

        // Create and initialize find_intersected_objects_process
        //TODO SPEED UP only search in sub model part where it's likely that skin points are located? --> several layers of elements surrounding the previous SBM_BOUNDARY
        Flags find_intersected_options;
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_CONDITIONS, false);
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_ELEMENTS, true);
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS, false);
        find_intersected_options.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS, true);
        FindIntersectedGeometricalObjectsProcess find_intersected_objects_process = FindIntersectedGeometricalObjectsProcess(*mpModelPart, *mpSkinModelPart, find_intersected_options);
        find_intersected_objects_process.ExecuteInitialize();

        // Find intersections
        find_intersected_objects_process.FindIntersections();  // FindIntersectedSkinObjects marks intersected elements as SELECTED
        const std::vector<PointerVector<GeometricalObject>>&  r_intersected_objects = find_intersected_objects_process.GetIntersections();

        // Get elements in the same order as find_intersected_objects_process
        const std::size_t n_elements = (find_intersected_objects_process.GetModelPart1()).NumberOfElements();
        const auto& r_elements = (find_intersected_objects_process.GetModelPart1()).ElementsArray();

        // Define the minimal distance to the skin geometry
        // NOTE that a minimal distance of <5e-7 can already cause an ill-defined deactivated layer because of the tolerances of the geometrical operations (intersected elements and localization of skin points)
        const double distance_threshold = 5e-6;
        // Initialize a multiplier for the the distance by which nodes should be relocated --> length/1000 for a smooth deactivated element layer
        const double relocation_multiplier = 1e-3;

        // Get function pointer for calculating the distance between a node and an intersecting object
        PointDistanceFunctionType p_point_distance_function;
        IntersectionPlaneConstructorType p_intersection_plane_constructor;
        const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
        switch (n_dim) {
            case 2:
                p_point_distance_function = ShiftedBoundaryUtilityInternals::CalculatePointDistance<2>;
                p_intersection_plane_constructor = ShiftedBoundaryUtilityInternals::CreateIntersectionPlane<2>;
                break;
            case 3:
                p_point_distance_function = ShiftedBoundaryUtilityInternals::CalculatePointDistance<3>;
                p_intersection_plane_constructor = ShiftedBoundaryUtilityInternals::CreateIntersectionPlane<3>;
                break;
            default:
                KRATOS_ERROR << "Wrong domain size.";
        }

        // Calculate nodal distances of all nodes of elements that are intersected and relocate nodes that are too close to the skin geometry
        std::size_t n_relocated_nodes = 0;
        LockObject mutex;
        IndexPartition<std::size_t>(n_elements).for_each([&](std::size_t i_ele) {
            const auto& r_element_intersecting_objects = r_intersected_objects[i_ele];
            if (!r_element_intersecting_objects.empty()) {
                auto& p_element = r_elements[i_ele];
                auto& r_geom = p_element->GetGeometry();
                const std::size_t n_nodes = r_geom.PointsNumber();

                // Calculate distance below which nodes get relocated and by which nodes get relocated
                const double length = r_geom.Length();
                const double relocation_distance = std::max(distance_threshold, relocation_multiplier * length);

                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                    auto& r_node = r_geom[i_node];
                    // Check that node was not already modified from the iteration of another element
                    if (!r_node.Is(MODIFIED)) {
                        // Calculate minimal distance of node to skin geometry
                        // NOTE that this distance value can be negative to move in negative normal direction
                        const double distance_node = CalculateSkinDistanceToNode(r_node, r_element_intersecting_objects,
                                                                                 p_point_distance_function, p_intersection_plane_constructor,
                                                                                 relocation_distance, 0.01*relocation_distance);

                        // Relocate node in direction of normal if the skin geometry is too close
                        //TODO do not move outer boundary elements?
                        if (std::abs(distance_node) < relocation_distance) {
                            // Get normal of first intersecting object
                            const auto& r_int_obj = r_intersected_objects[i_ele][0];
                            array_1d<double,3> local_coords;
                            array_1d<double,3> unit_normal = r_int_obj.GetGeometry().Normal(local_coords);
                            unit_normal /= norm_2(unit_normal);

                            // Use negative normal direction for a negative distance value (so it moves away from skin geometry)
                            if (distance_node < 0.0) {
                                unit_normal *= (-1);
                            }

                            {
                                std::scoped_lock<LockObject> lock(mutex);

                                // Relocate node into the direction of the intersected object's (negative) normal
                                r_node.X()  += relocation_distance * unit_normal[0];
                                r_node.Y()  += relocation_distance * unit_normal[1];
                                r_node.Z()  += relocation_distance * unit_normal[2];

                                // NOTE that initial position (X0,Y0,Z0) is mesh output for visualization (ParaView)
                                //TODO store also in deformation etc instead of initial position for visualization for moving geometry?
                                r_node.X0() += relocation_distance * unit_normal[0];
                                r_node.Y0() += relocation_distance * unit_normal[1];
                                r_node.Z0() += relocation_distance * unit_normal[2];
                                r_node.Set(MODIFIED, true);

                                //TODO correction needed for acceleration of the node?? correction needed for FM-ALE??
                                // mesh velocity += dx/dt (1st order approximation)

                                n_relocated_nodes++;
                            }
                            //TODO: Check detJ of surrounding elements - if negative, then revert?
                        }
                    }
                }
            }
        });

        // Recalculate intersections if nodes were moved
        // NOTE that entirely new process is necessary for taking update node positions into consideration (TODO change that??)
        //TODO SPEED UP by calculating intersections only for elements surrounding the relocated nodes?
        if (n_relocated_nodes > 0) {
            FindIntersectedGeometricalObjectsProcess new_find_intersected_objects_process = FindIntersectedGeometricalObjectsProcess(*mpModelPart, *mpSkinModelPart, find_intersected_options);
            new_find_intersected_objects_process.ExecuteInitialize();
            new_find_intersected_objects_process.FindIntersections();
            const std::vector<PointerVector<GeometricalObject>>&  r_new_intersected_objects = new_find_intersected_objects_process.GetIntersections();

            // // Mark elements as SBM_BOUNDARY that are intersected by the tessellated skin geometry
            // IndexPartition<std::size_t>(n_elements).for_each([&](std::size_t i_ele) {
            //     if (!r_new_intersected_objects[i_ele].empty()) {
            //         r_elements[i_ele]->Set(SBM_BOUNDARY, true);
            //     }
            // });
            for (std::size_t i_ele = 0; i_ele < n_elements; ++i_ele) {
                if (!r_new_intersected_objects[i_ele].empty()) {
                    auto p_elem = r_elements[i_ele];
                    p_elem->Set(SBM_BOUNDARY, true);
                    for (auto& r_node : p_elem->GetGeometry()) {
                        r_node.Set(SBM_BOUNDARY, false);
                        r_node.Set(SBM_INTERFACE, false);
                    }
                }
            }
        } else {
            // // Mark elements as SBM_BOUNDARY that are intersected by the tessellated skin geometry
            // IndexPartition<std::size_t>(n_elements).for_each([&](std::size_t i_ele) {
            //     if (!r_intersected_objects[i_ele].empty()) {
            //         r_elements[i_ele]->Set(SBM_BOUNDARY, true);
            //     }
            // });
            for (std::size_t i_ele = 0; i_ele < n_elements; ++i_ele) {
                if (!r_intersected_objects[i_ele].empty()) {
                    auto p_elem = r_elements[i_ele];
                    p_elem->Set(SBM_BOUNDARY, true);
                    for (auto& r_node : p_elem->GetGeometry()) {
                        r_node.Set(SBM_BOUNDARY, false);
                        r_node.Set(SBM_INTERFACE, false);
                    }
                }
            }
        }

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "'" << mSkinModelPartName << "' tessellated boundary flags were set. " << n_relocated_nodes << " nodes were relocated." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::LocateSkinPoints()
    {
        mSkinPointsMap.clear();
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

        // Mark elements as SBM_BOUNDARY (gamma) in which skin points were located
        // NOTE that these elements should already be SBM_BOUNDARY because of the tessellated interface
        // NOTE that SBM_BOUNDARY elements will get deactivated and that the split elements BC is applied by means of the extension operators
        // block_for_each(mSkinPointsMap.begin(), mSkinPointsMap.end(), [](const std::pair<ElementType::Pointer, SkinPointsDataVectorType>& rKeyData){
        //    rKeyData.first->Set(SBM_BOUNDARY, true);
        // });
        for (auto& [p_elem, skin_pt_data] : mSkinPointsMap) {
            p_elem->Set(SBM_BOUNDARY, true);
            for (auto& r_node : p_elem->GetGeometry()) {
                r_node.Set(SBM_BOUNDARY, false);
                r_node.Set(SBM_INTERFACE, false);
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetInterfaceFlags()
    {
        // Find the surrogate boundary elements and mark them as SBM_INTERFACE (gamma_tilde)
        // Note that we rely on the fact that the neighbors are sorted according to the faces
        // TODO faster to store SBM_BOUNDARY element pointers somewhere?
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
        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Interface flags were set."  << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::DeactivateElementsAndNodes(
        const bool DeactivateUnstableClusters
    )
    {
        // Deactivate elements in which the (true) boundary is located
        block_for_each(mpModelPart->Elements(), [](ElementType& rElement){
            if (rElement.Is(SBM_BOUNDARY)) {
                rElement.Set(ACTIVE, false);
            }
        });
        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Boundary elements were deactivated." << std::endl;

        if (DeactivateUnstableClusters) {
            FindAndDeactivateUnstableClusters();
            KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Unstable Clusters were deactivated." << std::endl;
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

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Single nodes were deactivated." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateAndAddSkinIntegrationPointConditions()
    {
        // Iterate over the split elements to create a vector for each element defining both sides using the element's skin points normals and positions
        // The resulting vector is as long as the number of nodes of the element and a positive value stands for the positive side of the boundary, a negative one for the negative side.
        // Also store the average position and average normal of the skin points located in the element
        // Nodes on the positive side will be declared SBM_BOUNDARY and nodes on the negative side SBM_INTERFACE which is used in the creation of the support clouds
        mSidesVectorMap.clear();
        AverageSkinToElementsMapType avg_skin_map;
        SetSidesVectorsAndSkinNormalsForSplitElements(mSkinPointsMap, mSidesVectorMap, avg_skin_map);
        //KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Sides vectors and skin normals were set." << std::endl;

        // Iterate over the split elements to create an extension basis for each node of the element (MLS shape functions values for support cloud of node)
        // NOTE that no extension bases will be calculated and added for a node for which not a sufficient number of support nodes were found
        mExtensionOperatorMap.clear();
        SetExtensionOperatorsForSplitElementNodes(mSidesVectorMap, avg_skin_map, mExtensionOperatorMap);
        //KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "Extension operators were set." << std::endl;

        // Set the pressure of the first node of an enclosed volume to zero if one side is enclosed.
        auto skin_pt_element_iter = mSkinPointsMap.begin();
        bool enclosed_pressure_is_set = false;
        //while (skin_pt_element_iter != mSkinPointsMap.end()) {
        if (mPositiveSideIsEnclosed || mNegativeSideIsEnclosed) {
            while (!enclosed_pressure_is_set && skin_pt_element_iter != mSkinPointsMap.end()) {
                auto p_elem = skin_pt_element_iter->first;
                const array_1d<double,3> avg_skin_position = avg_skin_map[p_elem].first;
                const array_1d<double,3> avg_skin_normal = avg_skin_map[p_elem].second;
                enclosed_pressure_is_set = SetEnclosedNodesPressure(*p_elem, mSidesVectorMap[p_elem], avg_skin_position, avg_skin_normal);
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
            CreateCloudNodeVectorsForSplitElement(*p_element, mSidesVectorMap[p_element], mExtensionOperatorMap, cloud_nodes_vector_pos, cloud_nodes_vector_neg);

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

                // Get the split element's shape function values and derivatives at the skin/ integration point
                Vector skin_pt_N(n_nodes);
                Matrix skin_pt_DN_DX = ZeroMatrix(n_nodes, n_dim);
                GetDataForSplitElementSkinPoint(*p_element, skin_pt_position, skin_pt_N, skin_pt_DN_DX);

                // Add skin pt. condition for positive side of boundary - using support cloud data for negative nodes
                // NOTE that the boundary normal is negative in order to point outward (from positive to negative side),
                // because positive side is where dot product of vector to node with average normal is positive
                n_skin_pt_conditions_added_pos += AddIntegrationPointCondition(*p_element, sides_vector, h, skin_pt_position, -skin_pt_area_normal,
                mExtensionOperatorMap, cloud_nodes_vector_pos, skin_pt_N, skin_pt_DN_DX, max_cond_id, /*ConsiderPositiveSide=*/true);

                // Add skin pt. condition for negative side of boundary - using support cloud data for positive nodes
                // NOTE that boundary normal is pointing outward (from negative to positive side)
                n_skin_pt_conditions_added_neg += AddIntegrationPointCondition(*p_element, sides_vector, h, skin_pt_position, skin_pt_area_normal,
                mExtensionOperatorMap, cloud_nodes_vector_neg, skin_pt_N, skin_pt_DN_DX, max_cond_id, /*ConsiderPositiveSide=*/false);
            }
        }
        if (n_skin_pt_conditions_added_pos != n_skin_points) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "Integration point conditions were NOT successfully added for the positive side of "
                << n_skin_points-n_skin_pt_conditions_added_pos << " skin points." << std::endl;
        }
        if (n_skin_pt_conditions_added_neg != n_skin_points) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "Integration point conditions were NOT successfully added for the negative side of "
                << n_skin_points-n_skin_pt_conditions_added_neg << " skin points." << std::endl;
        }
        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "'" << mSkinModelPartName << "' skin point conditions were added." << std::endl;
    }

    // void ShiftedBoundaryPointBasedInterfaceUtility::CalculatePressureAtSkinNodes()
    // {
    //     const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
    //     switch (n_dim) {
    //         case 2:
    //             CalculatePressureAtSkinNodesTemplated<2>();
    //             break;
    //         case 3:
    //             CalculatePressureAtSkinNodesTemplated<3>();
    //             break;
    //         default:
    //             KRATOS_ERROR << "Wrong domain size.";
    //     }
    // }

    // template <std::size_t TDim>
    // void ShiftedBoundaryPointBasedInterfaceUtility::CalculatePressureAtSkinNodesTemplated()
    // {
    //     // Create the bin-based point locator (TDim here makes template necessary)
    //     const std::size_t point_locator_max_results = 10000;
    //     const double point_locator_tolerance = 1.0e-5;
    //     BinBasedFastPointLocator<TDim> point_locator(*mpModelPart);
    //     point_locator.UpdateSearchDatabase();

    //     // Initialize counter for number of skin nodes that can not be calculated properly
    //     std::size_t n_skin_nodes_without_correct_extension = 0;

    //     // For each skin node interpolate the positive and negative side pressure from the fluid mesh
    //     LockObject mutex_valid_ex;
    //     block_for_each(mpSkinModelPart->Nodes(), [&](NodeType& rSkinNode){
    //         // Locate node in fluid mesh (should be inside a deactivated element)
    //         Vector skin_node_N(TDim+1); //TODO this restricts the use to tri and tetra
    //         Element::Pointer p_element = nullptr;
    //         typename BinBasedFastPointLocator<TDim>::ResultContainerType search_results(point_locator_max_results);
    //         const bool is_found = point_locator.FindPointOnMesh(
    //             rSkinNode.Coordinates(), skin_node_N, p_element,
    //             search_results.begin(), point_locator_max_results, point_locator_tolerance);

    //         // If the skin node is found, interpolate the positive and negative face pressure
    //         if (is_found) {
    //             // Get positive and negative side pressure variables of skin node
    //             double& r_skin_p_pos = rSkinNode.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
    //             double& r_skin_p_neg = rSkinNode.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);

    //             // Calculate pressure at skin node for positive and negative side of Gamma
    //             const bool p_calculated_successfully = CalculatePressureAtSplitElementSkinPoint(p_element, skin_node_N, r_skin_p_pos, r_skin_p_neg);
    //             if (!p_calculated_successfully) {
    //                 std::scoped_lock<LockObject> lock(mutex_valid_ex);
    //                 n_skin_nodes_without_correct_extension++;
    //             }
    //         }
    //     });
    //     KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", n_skin_nodes_without_correct_extension > 0)
    //         << "The pressure at " << n_skin_nodes_without_correct_extension << " skin nodes was calculated without valid extension." << std::endl;
    // }

    // void ShiftedBoundaryPointBasedInterfaceUtility::CalculateVelocityAtSkinNodes()
    // {
    //     const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
    //     switch (n_dim) {
    //         case 2:
    //             CalculateVelocityAtSkinNodesTemplated<2>();
    //             break;
    //         case 3:
    //             CalculateVelocityAtSkinNodesTemplated<3>();
    //             break;
    //         default:
    //             KRATOS_ERROR << "Wrong domain size.";
    //     }
    // }

    // template <std::size_t TDim>
    // void ShiftedBoundaryPointBasedInterfaceUtility::CalculateVelocityAtSkinNodesTemplated()
    // {
    //     // Create the bin-based point locator (TDim here makes template necessary)
    //     const std::size_t point_locator_max_results = 10000;
    //     const double point_locator_tolerance = 1.0e-5;
    //     BinBasedFastPointLocator<TDim> point_locator(*mpModelPart);
    //     point_locator.UpdateSearchDatabase();

    //     // Initialize counter for number of skin nodes that can not be calculated properly
    //     std::size_t n_skin_nodes_without_correct_extension = 0;

    //     // For each skin node interpolate the positive and negative side velocity from the fluid mesh
    //     LockObject mutex_valid_ex;
    //     block_for_each(mpSkinModelPart->Nodes(), [&](NodeType& rSkinNode){
    //         // Locate node in fluid mesh (should be inside a deactivated element)
    //         Vector skin_node_N(TDim+1); //TODO this restricts the use to tri and tetra
    //         Element::Pointer p_element = nullptr;
    //         typename BinBasedFastPointLocator<TDim>::ResultContainerType search_results(point_locator_max_results);
    //         const bool is_found = point_locator.FindPointOnMesh(
    //             rSkinNode.Coordinates(), skin_node_N, p_element,
    //             search_results.begin(), point_locator_max_results, point_locator_tolerance);

    //         // If the skin node is found, interpolate the positive and negative face fluid velocity
    //         if (is_found) {
    //             // Get positive and negative side velocity variables of skin node
    //             array_1d<double,3>& r_skin_u_pos = rSkinNode.FastGetSolutionStepValue(POSITIVE_FACE_FLUID_VELOCITY);
    //             array_1d<double,3>& r_skin_u_neg = rSkinNode.FastGetSolutionStepValue(NEGATIVE_FACE_FLUID_VELOCITY);

    //             // Calculate velocity at skin node for positive and negative side of Gamma
    //             const bool u_calculated_successfully = CalculateVelocityAtSplitElementSkinPoint(p_element, skin_node_N, r_skin_u_pos, r_skin_u_neg);
    //             if (!u_calculated_successfully) {
    //                 std::scoped_lock<LockObject> lock(mutex_valid_ex);
    //                 n_skin_nodes_without_correct_extension++;
    //             }
    //         }
    //     });
    //     KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", n_skin_nodes_without_correct_extension > 0)
    //         << "The velocity at " << n_skin_nodes_without_correct_extension << " skin nodes was calculated without valid extension." << std::endl;
    // }

    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPoints()
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

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPointsTemplated()
    {
        const std::size_t voigt_size = 3 * (TDim-1);
        const std::size_t block_size = TDim +1;

        array_1d<double, 3> force_on_skin = ZeroVector(3);

        std::size_t n_split_elements_without_correct_extension = 0;

        // Loop over all elements containing skin points to integrate traction=sigma*n over the interface
        // NOTE that both interface sides need to be integrated
        LockObject mutex, mutex_valid_ex;
        std::for_each(mSkinPointsMap.begin(), mSkinPointsMap.end(), [&](const std::pair<ElementType::Pointer, SkinPointsDataVectorType>& rKeyData){
            const auto p_element = rKeyData.first;
            const auto& skin_points_data_vector = rKeyData.second;
            const auto& r_geom = p_element->GetGeometry();
            const std::size_t n_nodes = r_geom.PointsNumber();
            const std::size_t local_size = n_nodes * block_size;

            // Calculate unknowns at the nodes of the element for the positive and the negative side of Gamma
            Vector unknowns_pos = ZeroVector(local_size);
            Vector unknowns_neg = ZeroVector(local_size);
            const bool unknowns_calculated_successfully = CalculateUnknownsForBothSidesOfSplitElement<TDim>(p_element, unknowns_pos, unknowns_neg);
            if (!unknowns_calculated_successfully) {
                std::scoped_lock<LockObject> lock(mutex_valid_ex);
                n_split_elements_without_correct_extension++;
            }

            // Get constitutive law, see FluidElementData and FluidElement::CalculateMaterialResponse
            const auto p_constitutive_law =  p_element->GetProperties()[CONSTITUTIVE_LAW];
            ConstitutiveLaw::Parameters constitutive_law_values(r_geom, p_element->GetProperties(), mpModelPart->GetProcessInfo());
            Flags& r_cl_options = constitutive_law_values.GetOptions();
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);  // not being used (?)
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            // Iterate over the element's skin points adding a positive side and a negative side drag contribution
            for (std::size_t i_skin_pt = 0; i_skin_pt < skin_points_data_vector.size(); ++i_skin_pt) {
                // Get the skin point's position and area normal
                const auto skin_pt_data = skin_points_data_vector[i_skin_pt];
                const array_1d<double,3> skin_pt_position = std::get<0>(skin_pt_data);
                const array_1d<double,3> skin_pt_area_normal = std::get<1>(skin_pt_data);
                const std::size_t skin_pt_node_id = std::get<2>(skin_pt_data);

                // Get normal and represented area of skin point
                // NOTE that the normal points in fluid inward direction for the positive side of the skin and fluid outward for the negative side of the skin
                const double skin_pt_area = norm_2(skin_pt_area_normal);
                const array_1d<double,3> aux_unit_normal = skin_pt_area_normal / skin_pt_area;

                // Get shape function values and derivatives of the element at the skin point
                Vector skin_pt_N(n_nodes);
                Matrix skin_pt_DN_DX = ZeroMatrix(n_nodes, TDim);
                GetDataForSplitElementSkinPoint(*p_element, skin_pt_position, skin_pt_N, skin_pt_DN_DX);

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

                // Add the shear stress and pressure drag contribution of the skin point to the skin drag of the skin geometry
                {
                    std::scoped_lock<LockObject> lock(mutex);
                    force_on_skin += skin_pt_force;
                }

                // Store velocity, pressure and traction in a skin point model part
                auto& skin_pt_in_model_part = mpSkinPointsModelPart->GetNode(skin_pt_node_id);
                skin_pt_in_model_part.FastGetSolutionStepValue(POSITIVE_FACE_FLUID_VELOCITY) = u_pos;
                skin_pt_in_model_part.FastGetSolutionStepValue(NEGATIVE_FACE_FLUID_VELOCITY) = u_neg;
                skin_pt_in_model_part.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = p_pos;
                skin_pt_in_model_part.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) = p_neg;
                skin_pt_in_model_part.FastGetSolutionStepValue(TRACTION_FROM_FLUID_PRESSURE) = traction_p;
                skin_pt_in_model_part.FastGetSolutionStepValue(TRACTION_FROM_FLUID_STRESS) = traction_tau;
                skin_pt_in_model_part.FastGetSolutionStepValue(DRAG_FORCE) = skin_pt_force;
            }
        });
        KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", n_split_elements_without_correct_extension > 0)
        << "The unknowns inside " << n_split_elements_without_correct_extension << " split elements were calculated without valid extension." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPointsAndNodes()
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

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPointsAndNodesTemplated()
    {
        // Calculate variables at all skin points (integration points) and store them in the skin points model part
        CalculateVariablesAtSkinPointsTemplated<TDim>();

        // Create dictionaries to accumulate the weighted values and weights for each skin node for interpolation
        std::unordered_map<std::size_t, array_1d<double,2>> dict_sum_weighted_p; //positive and negative side pressure
        std::unordered_map<std::size_t, array_1d<array_1d<double,3>,2>> dict_sum_weighted_u; //positive and negative side velocity
        // traction from pressure, traction from stress, drag force
        std::unordered_map<std::size_t, double> dict_sum_weights;
        for (auto& r_skin_node : mpSkinModelPart->Nodes()) {
            dict_sum_weighted_p[r_skin_node.Id()] = ZeroVector(2);
            dict_sum_weighted_u[r_skin_node.Id()][0] = ZeroVector(3);
            dict_sum_weighted_u[r_skin_node.Id()][1] = ZeroVector(3);
            dict_sum_weights[r_skin_node.Id()] = 0.0;
        }

        // Set the bin-based fast point locator
        //TODO faster to keep pointer_locator for all iterations?
        const std::size_t point_locator_max_results = 10000;
        const double point_locator_tolerance = 1.0e-5;
        BinBasedFastPointLocator<TDim> point_locator_skin(*mpSkinModelPart);
        point_locator_skin.UpdateSearchDatabase();


        // Loop over all skin points and locate them in the skin model part to add their contribution to the skin nodes
        LockObject mutex;
        block_for_each(mpSkinPointsModelPart->Nodes(), [&](NodeType& rSkinPoint){
            // Search for the skin point in the skin mesh to get the element containing the point and the shape function values
            array_1d<double,3> skin_pt_position = rSkinPoint.Coordinates();
            Vector aux_N(TDim+1); //Note that this restricts the use to tri and tetra
            Element::Pointer p_element = nullptr;
            typename BinBasedFastPointLocator<TDim>::ResultContainerType search_results(point_locator_max_results);
            const bool is_found = point_locator_skin.FindPointOnMesh(
                skin_pt_position, aux_N, p_element,
                search_results.begin(), point_locator_max_results, point_locator_tolerance);

            if (is_found) {
                // Get skin point variables
                const double& r_p_pos = rSkinPoint.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                const double& r_p_neg = rSkinPoint.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
                const array_1d<double,3>& r_u_pos = rSkinPoint.FastGetSolutionStepValue(POSITIVE_FACE_FLUID_VELOCITY);
                const array_1d<double,3>& r_u_neg = rSkinPoint.FastGetSolutionStepValue(NEGATIVE_FACE_FLUID_VELOCITY);
                // const array_1d<double,3>& r_traction_p = rSkinPoint.FastGetSolutionStepValue(TRACTION_FROM_FLUID_PRESSURE);
                // const array_1d<double,3>& r_traction_tau = rSkinPoint.FastGetSolutionStepValue(TRACTION_FROM_FLUID_STRESS);
                // const array_1d<double,3>& r_drag_force = rSkinPoint.FastGetSolutionStepValue(DRAG_FORCE);

                // Loop over the nodes of the element to add the skin point's contribution to each node
                const auto& r_geom = p_element->GetGeometry();
                for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                    const IndexType node_id = r_geom[i_node].Id();
                    const double N_i = aux_N(i_node);

                    std::scoped_lock<LockObject> lock(mutex);
                    dict_sum_weighted_p[node_id][0] += N_i * r_p_pos;
                    dict_sum_weighted_p[node_id][1] += N_i * r_p_neg;
                    dict_sum_weighted_u[node_id][0] += N_i * r_u_pos;
                    dict_sum_weighted_u[node_id][1] += N_i * r_u_neg;
                    // dict_sum_weighted_vectors[node_id][2] += N_i * r_traction_p;
                    // dict_sum_weighted_vectors[node_id][3] += N_i * r_traction_tau;
                    // dict_sum_weighted_vectors[node_id][4] += N_i * r_drag_force;
                    dict_sum_weights[node_id] += N_i;
                }
            } else {
                KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "Skin point " << rSkinPoint.Id() << " could not be located in the skin mesh." << std::endl;
            }
        });

        // Set the skin nodes variables based on the skin points contributions
        block_for_each(mpSkinModelPart->Nodes(), [&](NodeType& rSkinNode){
            const auto& sum_weighted_p = dict_sum_weighted_p[rSkinNode.Id()];
            const auto& sum_weighted_u = dict_sum_weighted_u[rSkinNode.Id()];
            const double& sum_weights = dict_sum_weights[rSkinNode.Id()];

            rSkinNode.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE)       = sum_weighted_p[0] / sum_weights;
            rSkinNode.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE)       = sum_weighted_p[1] / sum_weights;
            rSkinNode.FastGetSolutionStepValue(POSITIVE_FACE_FLUID_VELOCITY) = sum_weighted_u[0] / sum_weights;
            rSkinNode.FastGetSolutionStepValue(NEGATIVE_FACE_FLUID_VELOCITY) = sum_weighted_u[1] / sum_weights;
            // rSkinNode.FastGetSolutionStepValue(TRACTION_FROM_FLUID_PRESSURE) = sum_weighted_vectors[2]   / sum_weights;
            // rSkinNode.FastGetSolutionStepValue(TRACTION_FROM_FLUID_STRESS)   = sum_weighted_vectors[3]   / sum_weights;
            // rSkinNode.FastGetSolutionStepValue(DRAG_FORCE)                   = sum_weighted_vectors[4]  / sum_weights;
        });
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::CalculateExtensionError()
    {
        /*
        std::vector<double> x_pos, y_pos, error_u_pos, error_u_neg, error_p_pos, error_p_neg, error_dpX_pos, error_dpY_pos, error_dpX_neg, error_dpY_neg;

        std::size_t n_skin_points_without_correct_extension = 0;

        // Loop over all elements containing skin points
        for (auto& rKeyData : mSkinPointsMap) {
            const auto p_element = rKeyData.first;
            const auto& skin_points_data_vector = rKeyData.second;
            const auto& r_geom = p_element->GetGeometry();
            const std::size_t n_dim = r_geom.WorkingSpaceDimension();
            const std::size_t n_nodes = r_geom.PointsNumber();

            // Iterate over the element's skin points adding a positive side and a negative side drag contribution
            for (std::size_t i_skin_pt = 0; i_skin_pt < skin_points_data_vector.size(); ++i_skin_pt) {
                // Get the skin point's position and area normal
                const auto skin_pt_data = skin_points_data_vector[i_skin_pt];
                const array_1d<double,3> skin_pt_position = std::get<0>(skin_pt_data);
                const array_1d<double,3> skin_pt_area_normal = std::get<1>(skin_pt_data);

                // Get shape function values of the element at the skin point
                Vector skin_pt_N(n_nodes);
                Matrix skin_pt_DN_DX = ZeroMatrix(n_nodes, n_dim);
                GetDataForSplitElementSkinPoint(*p_element, skin_pt_position, skin_pt_N, skin_pt_DN_DX);

                // Calculate the velocity at the skin point for positive and negative side of Gamma
                array_1d<double,3> u_pos, u_neg;
                const bool u_calculated_successfully = CalculateVelocityAtSplitElementSkinPoint(p_element, skin_pt_N, u_pos, u_neg);
                if (!u_calculated_successfully) {
                    n_skin_points_without_correct_extension++;
                }

                // Calculate the pressure at the skin point for positive and negative side of Gamma
                double p_pos, p_neg;
                const bool p_calculated_successfully = CalculatePressureAtSplitElementSkinPoint(p_element, skin_pt_N, p_pos, p_neg);
                if (!p_calculated_successfully) {
                    n_skin_points_without_correct_extension++;
                }

                // Calculate the pressure gradient at the skin point for positive and negative side of Gamma
                Vector dp_pos, dp_neg;
                const bool dp_calculated_successfully = CalculatePressureGradientAtSplitElementSkinPoint(p_element, skin_pt_DN_DX, dp_pos, dp_neg);

                // Analytical solution of the hydrostatic reservoir
                const array_1d<double,3> u_analyt = ZeroVector(3);
                const double p_analyt = 100.0 * (1.0 - skin_pt_position[1]);
                const double dpX_analyt =    0.0;
                const double dpY_analyt = -100.0;

                // Calculate error of velocity and pressure for positive and negative side
                const double skin_pt_error_u_pos   = norm_2(u_pos - u_analyt);
                const double skin_pt_error_u_neg   = norm_2(u_neg - u_analyt);
                const double skin_pt_error_p_pos   = std::abs(p_pos - p_analyt);
                const double skin_pt_error_p_neg   = std::abs(p_neg - p_analyt);
                const double skin_pt_error_dpX_pos = std::abs(dp_pos(0) - dpX_analyt);
                const double skin_pt_error_dpY_pos = std::abs(dp_pos(1) - dpY_analyt);
                const double skin_pt_error_dpX_neg = std::abs(dp_neg(0) - dpX_analyt);
                const double skin_pt_error_dpY_neg = std::abs(dp_neg(1) - dpY_analyt);

                // Add data to vectors
                x_pos.push_back(skin_pt_position[0]);
                y_pos.push_back(skin_pt_position[1]);
                error_u_pos.push_back(skin_pt_error_u_pos);
                error_u_neg.push_back(skin_pt_error_u_neg);
                error_p_pos.push_back(skin_pt_error_p_pos);
                error_p_neg.push_back(skin_pt_error_p_neg);
                error_dpX_pos.push_back(skin_pt_error_dpX_pos);
                error_dpY_pos.push_back(skin_pt_error_dpY_pos);
                error_dpX_neg.push_back(skin_pt_error_dpX_neg);
                error_dpY_neg.push_back(skin_pt_error_dpY_neg);
            }
        }
            KRATOS_WARNING_IF("ShiftedBoundaryPointBasedInterfaceUtility", n_skin_points_without_correct_extension > 0)
            << "The pressure at " << n_skin_points_without_correct_extension << " skin points was calculated without valid extension." << std::endl;

        KRATOS_WATCH(x_pos);
        KRATOS_WATCH(y_pos);
        KRATOS_WATCH(error_u_pos);
        KRATOS_WATCH(error_u_neg);
        KRATOS_WATCH(error_p_pos);
        KRATOS_WATCH(error_p_neg);
        KRATOS_WATCH(error_dpX_pos);
        KRATOS_WATCH(error_dpY_pos);
        KRATOS_WATCH(error_dpX_neg);
        KRATOS_WATCH(error_dpY_neg);
        */
    }

    template <std::size_t TDim>
    bool ShiftedBoundaryPointBasedInterfaceUtility::CalculateUnknownsForBothSidesOfSplitElement(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns)
    {
        //TODO ? Check whether the element is actually a part of the boundary
        // if (!pElement->Is(SBM_BOUNDARY)) {
        //     return 0;
        // }
        const std::size_t block_size = TDim +1;
        const auto& r_geom = pElement->GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();

        Vector sides_vector(n_nodes);
        const std::size_t sides_found = mSidesVectorMap.count(pElement);
        // Some deactivated elements containing Gamma might not have skin points and therefore no side vector (originally)
        if (!sides_found) {
            // Create sides vector based on the element's nodes belonging to SBM_BOUNDARY or SBM_INTERFACE, if possible.
            // Helpful for sparse distributions of skin points.
            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                auto& r_node = r_geom[i_node];
                if (r_node.Is(SBM_BOUNDARY) && !r_node.Is(SBM_INTERFACE)) {
                    sides_vector[i_node] =  1.0;
                } else if (r_node.Is(SBM_INTERFACE) && !r_node.Is(SBM_BOUNDARY)) {
                    sides_vector[i_node] = -1.0;
                } else { return false; }
            }
            // Store the created sides vector for future calculations
            mSidesVectorMap.insert(std::make_pair(pElement, sides_vector));
        } else {
            // Get sides vector for the element the skin node is located in
            sides_vector = mSidesVectorMap[pElement];
        }

        // Calculate positive and negative side unknowns at all nodes of the element
        bool element_is_without_positive_extension = false;
        bool element_is_without_negative_extension = false;
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
                } else { element_is_without_negative_extension = true; }
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
                } else { element_is_without_negative_extension = true; }
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

        if (element_is_without_positive_extension || element_is_without_negative_extension) { return false; }
        else { return true; }
    }

    // bool ShiftedBoundaryPointBasedInterfaceUtility::CalculatePressureAtSplitElementSkinPoint(
    //     const ElementType::Pointer pElement,
    //     const Vector& rPointShapeFunctionValues,
    //     double& rPositiveSidePressure,
    //     double& rNegativeSidePressure)
    // {
    //     //TODO ? Check whether the element is actually a part of the boundary
    //     // if (!pElement->Is(SBM_BOUNDARY)) {
    //     //     return 0;
    //     // }
    //     const auto& r_geom = pElement->GetGeometry();
    //     const std::size_t n_nodes = r_geom.PointsNumber();

    //     // Initialize pressure variables
    //     rPositiveSidePressure = 0.0;
    //     rNegativeSidePressure = 0.0;

    //     Vector sides_vector(n_nodes);
    //     const std::size_t sides_found = mSidesVectorMap.count(pElement);
    //     // Some deactivated elements containing Gamma might not have skin points and therefore no side vector (originally)
    //     if (!sides_found) {
    //         // Create sides vector based on the element's nodes belonging to SBM_BOUNDARY or SBM_INTERFACE, if possible.
    //         // Helpful for sparse distributions of skin points.
    //         for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
    //             auto& r_node = r_geom[i_node];
    //             if (r_node.Is(SBM_BOUNDARY) && !r_node.Is(SBM_INTERFACE)) {
    //                 sides_vector[i_node] =  1.0;
    //             } else if (r_node.Is(SBM_INTERFACE) && !r_node.Is(SBM_BOUNDARY)) {
    //                 sides_vector[i_node] = -1.0;
    //             } else { return false; }
    //         }
    //         // Store the created sides vector for future calculations
    //         mSidesVectorMap.insert(std::make_pair(pElement, sides_vector));
    //     } else {
    //         // Get sides vector for the element the skin node is located in
    //         sides_vector = mSidesVectorMap[pElement];
    //     }

    //     // Calculate positive and negative side pressure at skin node
    //     bool point_is_without_positive_extension = false;
    //     bool point_is_without_negative_extension = false;
    //     for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto p_node = r_geom(i_node);
    //         const double i_node_N = rPointShapeFunctionValues[i_node];

    //         // Check whether extension exists for node
    //         const std::size_t extension_found = mExtensionOperatorMap.count(p_node);

    //         // The extension needs to be used for nodes of the element which are located on the other side of the pressure that is being calculated
    //         if (sides_vector[i_node] > 0.0) {
    //             rPositiveSidePressure += i_node_N * p_node->FastGetSolutionStepValue(PRESSURE);
    //             if (extension_found) {
    //                 const auto& ext_op_data = mExtensionOperatorMap[p_node];
    //                 // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the pressure at the support node
    //                 for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
    //                     //const auto& r_support_node_data = *it_data;
    //                     const auto p_support_node = std::get<0>(*it_data);
    //                     const double weight_support_node = std::get<1>(*it_data);
    //                     rNegativeSidePressure += i_node_N * weight_support_node * p_support_node->FastGetSolutionStepValue(PRESSURE);
    //                 }
    //             } else { point_is_without_negative_extension = true; }
    //         } else {
    //             if (extension_found) {
    //                 const auto& ext_op_data = mExtensionOperatorMap[p_node];
    //                 // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the pressure at the support node
    //                 for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
    //                     //const auto& r_support_node_data = *it_data;
    //                     const auto p_support_node = std::get<0>(*it_data);
    //                     const double weight_support_node = std::get<1>(*it_data);
    //                     rPositiveSidePressure += i_node_N * weight_support_node * p_support_node->FastGetSolutionStepValue(PRESSURE);
    //                 }
    //             } else { point_is_without_positive_extension = true; }
    //             rNegativeSidePressure += i_node_N * p_node->FastGetSolutionStepValue(PRESSURE);
    //         }
    //     }
    //     if (point_is_without_positive_extension || point_is_without_negative_extension) { return false; }
    //     else { return true; }
    // }

    // bool ShiftedBoundaryPointBasedInterfaceUtility::CalculatePressureGradientAtSplitElementSkinPoint(
    //     const ElementType::Pointer pElement,
    //     const Matrix& rPointShapeFunctionDerivatives,
    //     Vector& rPositiveSidePressureGradient,
    //     Vector& rNegativeSidePressureGradient)
    // {
    //     //TODO ? Check whether the element is actually a part of the boundary
    //     // if (!pElement->Is(SBM_BOUNDARY)) {
    //     //     return 0;
    //     // }
    //     const auto& r_geom = pElement->GetGeometry();
    //     const std::size_t n_dim = r_geom.WorkingSpaceDimension();
    //     const std::size_t n_nodes = r_geom.PointsNumber();

    //     // Initialize pressure variables
    //     rPositiveSidePressureGradient = ZeroVector(n_dim);
    //     rNegativeSidePressureGradient = ZeroVector(n_dim);

    //     Vector sides_vector(n_nodes);
    //     const std::size_t sides_found = mSidesVectorMap.count(pElement);
    //     // Some deactivated elements containing Gamma might not have skin points and therefore no side vector (originally)
    //     // NOTE that this might be wrong if another skin geometry is embedded close by!
    //     if (!sides_found) {
    //         // Create sides vector based on the element's nodes belonging to SBM_BOUNDARY or SBM_INTERFACE, if possible.
    //         // Helpful for sparse distributions of skin points.
    //         for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
    //             auto& r_node = r_geom[i_node];
    //             if (r_node.Is(SBM_BOUNDARY) && !r_node.Is(SBM_INTERFACE)) {
    //                 sides_vector[i_node] =  1.0;
    //             } else if (r_node.Is(SBM_INTERFACE) && !r_node.Is(SBM_BOUNDARY)) {
    //                 sides_vector[i_node] = -1.0;
    //             } else { return false; }
    //         }
    //         // Store the created sides vector for future calculations
    //         mSidesVectorMap.insert(std::make_pair(pElement, sides_vector));
    //     } else {
    //         // Get sides vector for the element the skin node is located in
    //         sides_vector = mSidesVectorMap[pElement];
    //     }

    //     // Calculate positive and negative side pressure gradient at skin node
    //     bool point_is_without_positive_extension = false;
    //     bool point_is_without_negative_extension = false;
    //     for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto p_node = r_geom(i_node);

    //         // Check whether extension exists for node
    //         const std::size_t extension_found = mExtensionOperatorMap.count(p_node);

    //         // The extension needs to be used for nodes of the element which are located on the other side of the pressure that is being calculated
    //         if (sides_vector[i_node] > 0.0) {
    //             for (std::size_t d = 0; d < n_dim; ++d) {
    //                 rPositiveSidePressureGradient(d) += rPointShapeFunctionDerivatives(i_node, d) * p_node->FastGetSolutionStepValue(PRESSURE);
    //             }
    //             if (extension_found) {
    //                 const auto& ext_op_data = mExtensionOperatorMap[p_node];
    //                 // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the pressure at the support node
    //                 for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
    //                     //const auto& r_support_node_data = *it_data;
    //                     const auto p_support_node = std::get<0>(*it_data);
    //                     const double weight_support_node = std::get<1>(*it_data);
    //                     for (std::size_t d = 0; d < n_dim; ++d) {
    //                         rNegativeSidePressureGradient(d) += rPointShapeFunctionDerivatives(i_node, d) * weight_support_node * p_node->FastGetSolutionStepValue(PRESSURE);
    //                     }
    //                 }
    //             } else { point_is_without_negative_extension = true; }
    //         } else {
    //             if (extension_found) {
    //                 const auto& ext_op_data = mExtensionOperatorMap[p_node];
    //                 // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the pressure at the support node
    //                 for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
    //                     //const auto& r_support_node_data = *it_data;
    //                     const auto p_support_node = std::get<0>(*it_data);
    //                     const double weight_support_node = std::get<1>(*it_data);
    //                     for (std::size_t d = 0; d < n_dim; ++d) {
    //                         rPositiveSidePressureGradient(d) += rPointShapeFunctionDerivatives(i_node, d) * weight_support_node * p_node->FastGetSolutionStepValue(PRESSURE);
    //                     }
    //                 }
    //             } else { point_is_without_positive_extension = true; }
    //             for (std::size_t d = 0; d < n_dim; ++d) {
    //                 rNegativeSidePressureGradient(d) += rPointShapeFunctionDerivatives(i_node, d) * p_node->FastGetSolutionStepValue(PRESSURE);
    //             }
    //         }
    //     }
    //     if (point_is_without_positive_extension || point_is_without_negative_extension) { return false; }
    //     else { return true; }
    // }

    // bool ShiftedBoundaryPointBasedInterfaceUtility::CalculateVelocityAtSplitElementSkinPoint(
    //     const ElementType::Pointer pElement,
    //     const Vector& rPointShapeFunctionValues,
    //     array_1d<double,3>& rPositiveSideVelocity,
    //     array_1d<double,3>& rNegativeSideVelocity)
    // {
    //     //TODO ? Check whether the element is actually a part of the boundary
    //     // if (!pElement->Is(SBM_BOUNDARY)) {
    //     //     return 0;
    //     // }
    //     const auto& r_geom = pElement->GetGeometry();
    //     const std::size_t n_nodes = r_geom.PointsNumber();

    //     // Initialize velocity variables
    //     rPositiveSideVelocity = ZeroVector(3);
    //     rNegativeSideVelocity = ZeroVector(3);

    //     Vector sides_vector(n_nodes);
    //     const std::size_t sides_found = mSidesVectorMap.count(pElement);
    //     // Some deactivated elements containing Gamma might not have skin points and therefore no side vector (originally)
    //     // NOTE that this might be wrong if another skin geometry is embedded close by!
    //     if (!sides_found) {
    //         // Create sides vector based on the element's nodes belonging to SBM_BOUNDARY or SBM_INTERFACE, if possible.
    //         // Helpful for sparse distributions of skin points.
    //         for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
    //             auto& r_node = r_geom[i_node];
    //             if (r_node.Is(SBM_BOUNDARY) && !r_node.Is(SBM_INTERFACE)) {
    //                 sides_vector[i_node] =  1.0;
    //             } else if (r_node.Is(SBM_INTERFACE) && !r_node.Is(SBM_BOUNDARY)) {
    //                 sides_vector[i_node] = -1.0;
    //             } else { return false; }
    //         }
    //         // Store the created sides vector for future calculations
    //         mSidesVectorMap.insert(std::make_pair(pElement, sides_vector));
    //     } else {
    //         // Get sides vector for the element the skin node is located in
    //         sides_vector = mSidesVectorMap[pElement];
    //     }

    //     // Calculate positive and negative side velocity at skin node
    //     bool point_is_without_positive_extension = false;
    //     bool point_is_without_negative_extension = false;
    //     for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto p_node = r_geom(i_node);
    //         const double i_node_N = rPointShapeFunctionValues[i_node];

    //         // Check whether extension exists for node
    //         const std::size_t extension_found = mExtensionOperatorMap.count(p_node);

    //         // The extension needs to be used for nodes of the element which are located on the other side of the velocity that is being calculated
    //         if (sides_vector[i_node] > 0.0) {
    //             rPositiveSideVelocity += i_node_N * p_node->FastGetSolutionStepValue(VELOCITY);
    //             if (extension_found) {
    //                 const auto& ext_op_data = mExtensionOperatorMap[p_node];
    //                 // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the velocity at the support node
    //                 for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
    //                     //const auto& r_support_node_data = *it_data;
    //                     const auto p_support_node = std::get<0>(*it_data);
    //                     const double weight_support_node = std::get<1>(*it_data);
    //                     rNegativeSideVelocity += i_node_N * weight_support_node * p_support_node->FastGetSolutionStepValue(VELOCITY);
    //                 }
    //             } else { point_is_without_negative_extension = true; }
    //         } else {
    //             if (extension_found) {
    //                 const auto& ext_op_data = mExtensionOperatorMap[p_node];
    //                 // Iterate over the node's extension operator data and add the support node weight (i_cl_node_N) times the velocity at the support node
    //                 for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
    //                     //const auto& r_support_node_data = *it_data;
    //                     const auto p_support_node = std::get<0>(*it_data);
    //                     const double weight_support_node = std::get<1>(*it_data);
    //                     rPositiveSideVelocity += i_node_N * weight_support_node * p_support_node->FastGetSolutionStepValue(VELOCITY);
    //                 }
    //             } else { point_is_without_positive_extension = true; }
    //             rNegativeSideVelocity += i_node_N * p_node->FastGetSolutionStepValue(VELOCITY);
    //         }
    //     }
    //     if (point_is_without_positive_extension || point_is_without_negative_extension) { return false; }
    //     else { return true; }
    // }

    double ShiftedBoundaryPointBasedInterfaceUtility::CalculateSkinDistanceToNode(
        const Node& rNode,
        const PointerVector<GeometricalObject>& rIntersectingObjects,
        PointDistanceFunctionType pPointDistanceFunction,
        IntersectionPlaneConstructorType pIntersectionPlaneConstructor,
        const double& DistanceThreshold,
        const double& ThresholdForSignedness)
    {
        double minimal_distance = 1.0;

        // Compute distance of node to each given intersecting object
        for (const auto& it_int_obj : rIntersectingObjects.GetContainer()) {
            const auto& r_int_obj_geom = it_int_obj->GetGeometry();
            double distance = pPointDistanceFunction(r_int_obj_geom, rNode);

            // Check that the computed distance is the minimal one so far
            if (distance < minimal_distance) {
                minimal_distance = distance;
                // Create plane to calculate signedness of distance if the distance value is
                // below the minimal accepted distance and above the threshold at which the signedness makes a difference for resulting geometry
                // In between these values the signedness ensures that the node moves away and not towards the skin geometry
                if (minimal_distance < DistanceThreshold && minimal_distance > ThresholdForSignedness) {
                    std::vector<array_1d<double,3>> plane_pts;
                    for (std::size_t i_obj_node = 0; i_obj_node < r_int_obj_geom.PointsNumber(); ++i_obj_node){
                        plane_pts.push_back(r_int_obj_geom[i_obj_node]);
                    }
                    Plane3D plane = pIntersectionPlaneConstructor(plane_pts);

                    if (plane.CalculateSignedDistance(rNode) < 0.0){
                        minimal_distance *= (-1);
                    }
                }
            }
        }

        return minimal_distance;
    }

    template <std::size_t TDim>
    void ShiftedBoundaryPointBasedInterfaceUtility::MapSkinPointsToElements(
        SkinPointsToElementsMapType& rSkinPointsMap)
    {
        //TODO SPEED UP
        //TODO faster to only search in sub model part where it's likely that skin points are located?
        // --> several layers of elements surrounding the previous SBM_BOUNDARY

        // Check that the skin model part has elements
        KRATOS_ERROR_IF_NOT(mpSkinModelPart->NumberOfElements())
            << "There are no elements in skin model part (boundary) '" << mpSkinModelPart->FullName() << "'." << std::endl;

        // Set the bin-based fast point locator
        //TODO faster to keep pointer_locator for all skin model parts and iterations?
        //TODO UpdateSearchDatabase necessary whenever there is a node relocation?
        const std::size_t point_locator_max_results = 10000;
        const double point_locator_tolerance = 1.0e-5;
        BinBasedFastPointLocator<TDim> point_locator(*mpModelPart);
        point_locator.UpdateSearchDatabase();

        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;

        //const std::size_t n_gp_per_element = mpSkinModelPart->ElementsBegin()->GetGeometry().IntegrationPointsNumber(integration_method);
        //const std::size_t n_skin_points = mpSkinModelPart->NumberOfElements() * n_gp_per_element;
        std::vector<Element::Pointer> skin_point_located_elements;
        std::vector<array_1d<double,3>> skin_point_positions;
        std::vector<array_1d<double,3>> skin_point_normals;

        // Search the skin points (skin model part integration points) in the volume mesh elements
        std::size_t n_skin_points_not_found = 0;
        std::size_t n_skin_points_found = 0;
        LockObject mutex;
        block_for_each(mpSkinModelPart->Elements(), [&](ElementType& rSkinElement){
            const auto& r_skin_geom = rSkinElement.GetGeometry();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_skin_geom.IntegrationPoints(integration_method);
            // Get detJ for all integration points of the skin element
            Vector int_pt_detJs;
            r_skin_geom.DeterminantOfJacobian(int_pt_detJs, integration_method);

            for (std::size_t i_int_pt = 0; i_int_pt < integration_points.size(); ++i_int_pt) {
                // Get position of skin point
                const array_1d<double,3> skin_pt_local_coords = integration_points[i_int_pt].Coordinates();
                array_1d<double,3> skin_pt_position = ZeroVector(3);
                r_skin_geom.GlobalCoordinates(skin_pt_position, skin_pt_local_coords);

                // Get normal at the skin point and make its length a measure of the area/ integration weight
                array_1d<double,3> skin_pt_area_normal = r_skin_geom.Normal(skin_pt_local_coords);
                // Normalize normal
                skin_pt_area_normal /= std::max(norm_2(skin_pt_area_normal), 1e-10);  // tolerance = std::pow(1e-3 * h, Dim-1)
                // Scale normal with integration weight
                const double int_pt_weight = int_pt_detJs[i_int_pt] * integration_points[i_int_pt].Weight();
                skin_pt_area_normal *= int_pt_weight;

                // Search for the skin point in the volume mesh to get the element containing the point
                Vector aux_N(TDim+1); //TODO this restricts the use to tri and tetra
                Element::Pointer p_element = nullptr;
                typename BinBasedFastPointLocator<TDim>::ResultContainerType search_results(point_locator_max_results);
                const bool is_found = point_locator.FindPointOnMesh(
                    skin_pt_position, aux_N, p_element,
                    search_results.begin(), point_locator_max_results, point_locator_tolerance);

                // Add data to vectors (only one thread is allowed here at a time)
                {
                    std::scoped_lock<LockObject> lock(mutex);
                    if (is_found){
                        skin_point_located_elements.push_back(p_element);
                        skin_point_positions.push_back(skin_pt_position);
                        skin_point_normals.push_back(skin_pt_area_normal);
                        n_skin_points_found++;
                    } else {
                        n_skin_points_not_found++;
                    }
                }
            }
        });
        if (n_skin_points_not_found > 0) {
            KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility")
                << n_skin_points_not_found << " skin points have not been found in any volume model part element." << std::endl;
        }

        const std::size_t n_nodes_skin_point_model_part = mpSkinPointsModelPart->NumberOfNodes();
        for (std::size_t i_skin_pt = 0; i_skin_pt < skin_point_located_elements.size(); ++i_skin_pt) {
            auto p_element =  skin_point_located_elements[i_skin_pt];
            auto& skin_pt_position = skin_point_positions[i_skin_pt];
            auto& skin_pt_normal = skin_point_normals[i_skin_pt];
            const std::size_t skin_pt_node_id = n_nodes_skin_point_model_part + i_skin_pt;

            // Check if skin point data already exists for the element in which the skin point is located in the map
            if (rSkinPointsMap.find(p_element) == rSkinPointsMap.end()) {
                // If element is not found, add the first skin point to the new key
                SkinPointsDataVectorType skin_points_data_vector(1);
                skin_points_data_vector[0] = std::make_tuple(skin_pt_position, skin_pt_normal, skin_pt_node_id);
                auto skin_points_key_data = std::make_pair(p_element, skin_points_data_vector);
                rSkinPointsMap.insert(skin_points_key_data);
            } else {
                // If element is found, resize the skin point data vector of the element and add the new skin point
                SkinPointsDataVectorType& skin_points_data_vector = rSkinPointsMap[p_element];
                const std::size_t n_skin_points_in_element = skin_points_data_vector.size();
                skin_points_data_vector.resize(n_skin_points_in_element+1);
                skin_points_data_vector[n_skin_points_in_element] = std::make_tuple(skin_pt_position, skin_pt_normal, skin_pt_node_id);
            }

            // Add skin point to skin point model part
            mpSkinPointsModelPart->CreateNewNode(skin_pt_node_id, skin_pt_position[0], skin_pt_position[1], skin_pt_position[2]);
        }

        KRATOS_INFO("ShiftedBoundaryPointBasedInterfaceUtility") << "'" << mSkinModelPartName << "' skin points ("
            << n_skin_points_found << ") were mapped to volume mesh elements (" << rSkinPointsMap.size() << ")." << std::endl;
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::FindAndDeactivateUnstableClusters()
    {
        // Initialize the cluster ID of all nodes and its counter
        // NOTE ACTIVATION_LEVEL is (misused?) here as cluster ID
        block_for_each(mpModelPart->Nodes(), [&](NodeType& rNode){
            rNode.FastGetSolutionStepValue(ACTIVATION_LEVEL) = 0;
        });

        // Initialize a map for all clusters and a set of starting nodes
        ClustersMapType clusters_map;
        std::size_t max_cluster_id = 1;
        NodesSetType starting_nodes;

        // Find all nodes of boundary elements and increase the ACTIVATION_LEVEL and take some nodes at the boundary as starting nodes
        // NOTE flag ACTIVE could be used here instead of SBM_BOUNDARY.
        std::size_t seed_count = 0;
        LockObject mutex;
        block_for_each(mpModelPart->Elements(), [&](ElementType& rElement){
            if (rElement.Is(SBM_BOUNDARY)) {
                for (NodeType& rNode : rElement.GetGeometry()) {
                    {
                        std::scoped_lock<LockObject> lock(mutex); //TODO is the locking necessary?
                        rNode.FastGetSolutionStepValue(ACTIVATION_LEVEL) = 1;
                        ++seed_count;
                        if (seed_count%100 == 0) {
                            starting_nodes.insert(&rNode);
                        }
                    }
                }
            }
        });

        // Create a new cluster for each starting node
        while (!starting_nodes.empty()) {
            // Get a starting node and erase it from the set
            auto it_node = starting_nodes.begin();
            NodeType::Pointer p_node = *it_node;
            starting_nodes.erase(it_node);

            // Create and add a new cluster for the starting node
            AddNewCluster(clusters_map, max_cluster_id, p_node);
        }

        // Frontier expansion from starting nodes for clustering the boundary
        //TODO parallel? is this performant?
        //std::for_each(clusters_map.begin(), clusters_map.end(), [&](std::pair< std::size_t, std::tuple<NodesSetType, ElementsSetType, std::unordered_set<std::size_t>> >& rKeyData){
            //const std::size_t cluster_id = rKeyData.first;
            // NodesSetType& boundary_nodes = std::get<0>(rKeyData.second);
            // std::unordered_set<std::size_t>& found_ids = std::get<2>(rKeyData.second);
            // AdvanceClusterAlongBoundary(cluster_id, boundary_nodes, found_ids, mutex);
        //});
        for (auto& [cluster_id, cluster_data]: clusters_map) {
            NodesSetType& boundary_nodes = std::get<0>(cluster_data);
            std::unordered_set<std::size_t>& found_ids = std::get<2>(cluster_data);
            AdvanceClusterAlongBoundary(cluster_id, boundary_nodes, found_ids);
        }


        // Find boundary nodes that are still not part of a cluster and collect them as starting nodes
        block_for_each(mpModelPart->Nodes(), [&](NodeType& rNode){
            if (rNode.FastGetSolutionStepValue(ACTIVATION_LEVEL) == 1) {
                std::scoped_lock<LockObject> lock(mutex);
                starting_nodes.insert(&rNode);
            }
        });

        // Create a new cluster for each starting node and expand to get the entire cluster
        while (!starting_nodes.empty()) {
            // Get a starting node
            auto it_node = starting_nodes.begin();
            NodeType::Pointer p_node = *it_node;

            // Create and add a new cluster for the starting node
            AddNewCluster(clusters_map, max_cluster_id, p_node);

            // Use frontier expansion to find all connected nodes
            auto& cluster_data = clusters_map[max_cluster_id];
            AdvanceClusterAlongBoundary(max_cluster_id, std::get<0>(cluster_data), std::get<2>(cluster_data));

            // Remove all seeds that are part of a cluster by now
            for (auto it = starting_nodes.begin(); it != starting_nodes.end(); ) {
                if ((*it)->FastGetSolutionStepValue(ACTIVATION_LEVEL) > 1.0) {
                    it = starting_nodes.erase(it);  // returns next valid iterator
                } else {
                    ++it;
                }
            }
        }


        // Merge connected clusters via found IDs to get unique clusters and remove clusters that consist of one single node
        // NOTE that only boundary nodes are merged here, because cluster elements are expected to be still empty.
        MergeConnectedClusters(clusters_map);

        // Loop over all clusters to flood fill the respective volume and check if any node of the cluster has a fixed dof and otherwise deactivate the cluster elements
        std::size_t biggest_cluster_id = 0;
        std::size_t biggest_cluster_n = 0;
        for (auto& [cluster_id, cluster_data]: clusters_map) {
            const NodesSetType& boundary_nodes = std::get<0>(cluster_data);
            ElementsSetType& cluster_elements = std::get<1>(cluster_data);

            // Flood fill the volume of the cluster until a fixed dof or all cluster elements were found
            const bool fixed_dof_found = FindClusterElementsUntilFixedDof(cluster_id, boundary_nodes, cluster_elements);
            if (fixed_dof_found) {
                KRATOS_WATCH(fixed_dof_found);
            }

            //TODO Right now this is called in the "Initialize" of the solver before any boundary conditions have been applied. Therefore, no fixed dofs are found for any cluster.
            //HACK Instead of deactivating clusters without any fixed dof, now all clusters will get deactivated except for the biggest one.
            if (cluster_elements.size() > biggest_cluster_n) {
                biggest_cluster_id = cluster_id;
                biggest_cluster_n = cluster_elements.size();
            }

            // Deactivate all elements of the cluster if no fixed degree of freedom was found
            // if (!fixed_dof_found) {
            //     for (auto& p_elem : cluster_elements) {
            //         p_elem->Set(ACTIVE, false);
            //     }
            // }
        }

        // KRATOS_WATCH(biggest_cluster_id);
        // KRATOS_WATCH(biggest_cluster_n);

        // Deactivate the elements of all clusters except for the elements of the biggest cluster
        for (auto& [cluster_id, cluster_data]: clusters_map) {
            if (cluster_id != biggest_cluster_id) {
                ElementsSetType& cluster_elements = std::get<1>(cluster_data);
                for (auto& p_elem : cluster_elements) {
                    p_elem->Set(ACTIVE, false);
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::AddNewCluster(
        ClustersMapType& ClustersMap,
        std::size_t& MaxClusterId,
        NodeType::Pointer pNode)
    {
        // NodesSetType boundary_nodes;
        // boundary_nodes.insert(pNode);
        // ElementsSetType cluster_elements;
        // std::unordered_set<std::size_t> found_ids;
        // pNode->FastGetSolutionStepValue(ACTIVATION_LEVEL) = ++MaxClusterId;
        // ClustersMap[MaxClusterId] = std::make_tuple(boundary_nodes, cluster_elements, found_ids);
        pNode->FastGetSolutionStepValue(ACTIVATION_LEVEL) = ++MaxClusterId;
        auto& cluster_data = ClustersMap[MaxClusterId];
        NodesSetType& boundary_nodes = std::get<0>(cluster_data);
        boundary_nodes.insert(pNode);
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::AdvanceClusterAlongBoundary(
        const std::size_t ClusterId,
        NodesSetType& rBoundaryNodes,
        std::unordered_set<std::size_t>& rFoundIds,
        LockObject& Mutex)
    {
        // Start frontier with starting node
        NodesSetType current_frontier;
        current_frontier.insert(*(rBoundaryNodes.begin()));

        // Search for new nodes of the cluster while there are nodes in the frontier
        while (!current_frontier.empty()) {
            // Get a node from the frontier and erase it from the set
            auto it_node = current_frontier.begin();
            NodeType::Pointer p_node = *it_node;
            current_frontier.erase(*it_node);  //TODO can I erase the entry without disturbing the pointer?

            // Get neighboring elements of the frontier node and check the cluster ID of the element's nodes if it is active
            // NOTE that boundary elements are required to be already deactivated
            auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);
            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                if (p_elem_neigh != nullptr) {
                    if (p_elem_neigh->Is(ACTIVE)) {
                        for (auto& r_neigh_node : p_elem_neigh->GetGeometry()) {
                            std::size_t node_cluster_id = 0;
                            // Check cluster ID of the node and add it to the cluster if it is not part of a cluster yet (locked)
                            {
                                std::scoped_lock<LockObject> lock(Mutex);
                                node_cluster_id = r_neigh_node.FastGetSolutionStepValue(ACTIVATION_LEVEL);
                                // If node is at the boundary and not part of a cluster yet, mark it as part of the current cluster
                                if (node_cluster_id == 1) {
                                    r_neigh_node.FastGetSolutionStepValue(ACTIVATION_LEVEL) = ClusterId;
                                }
                            }
                            // Add new cluster node to the set and the frontier
                            if (node_cluster_id == 1) {
                                rBoundaryNodes.insert(&r_neigh_node);
                                current_frontier.insert(&r_neigh_node);
                            // If node is already part of a cluster, store the information that both clusters are connected in found IDs
                            // NOTE that found ID is only added here if it is smaller, which is not necessary, but speeds up merging
                            } else if (node_cluster_id < ClusterId && node_cluster_id > 1) {
                                rFoundIds.insert(node_cluster_id);
                            }
                        }
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::AdvanceClusterAlongBoundary(
        const std::size_t ClusterId,
        NodesSetType& rBoundaryNodes,
        std::unordered_set<std::size_t>& rFoundIds)
    {
        // Start frontier with starting node
        NodesSetType current_frontier;
        current_frontier.insert(*(rBoundaryNodes.begin()));

        // Search for new nodes of the cluster while there are nodes in the frontier
        while (!current_frontier.empty()) {
            // Get a node from the frontier and erase it from the set
            auto it_node = current_frontier.begin();
            NodeType::Pointer p_node = *it_node;
            current_frontier.erase(*it_node);

            // Get neighboring elements of the frontier node and check the cluster ID of the element's nodes if it is active
            // NOTE that boundary elements are required to be already deactivated
            auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);
            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                if (p_elem_neigh != nullptr) {
                    if (p_elem_neigh->Is(ACTIVE)) {
                        for (auto& r_neigh_node : p_elem_neigh->GetGeometry()) {
                            std::size_t node_cluster_id = 0;
                            node_cluster_id = r_neigh_node.FastGetSolutionStepValue(ACTIVATION_LEVEL);

                            // If node is at the boundary and not part of a cluster yet, add it to the current cluster and the frontier
                            if (node_cluster_id == 1) {
                                r_neigh_node.FastGetSolutionStepValue(ACTIVATION_LEVEL) = ClusterId;
                                rBoundaryNodes.insert(&r_neigh_node);
                                current_frontier.insert(&r_neigh_node);
                            // If node is already part of a cluster, store the information that both clusters are connected in found IDs
                            } else if (node_cluster_id != ClusterId && node_cluster_id > 1) {
                                rFoundIds.insert(node_cluster_id);
                            }
                        }
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::MergeConnectedClusters(ClustersMapType& rClustersMap)
    {
        struct UnionFind {
            std::unordered_map<std::size_t, std::size_t> parent;
            std::unordered_map<std::size_t, std::size_t> rank;

            std::size_t find_parent(std::size_t id) {
                // If the ID is not yet in the map, initialize the parent as itself.
                if (parent.count(id) == 0) {
                    parent[id] = id;
                }
                // If the parent points to another Id assign that IDs parent recursively until the ID is its own parent (path compression).
                if (parent[id] != id) {
                    parent[id] = find_parent(parent[id]);
                }
                return parent[id];
            }

            void unite_parents(std::size_t id_A, std::size_t id_B) {
                std::size_t p_A = find_parent(id_A), p_B = find_parent(id_B);
                // Do nothing if parents are already the same.
                if (p_A == p_B) { return; }
                // Assign higher ranked parent as lower ranked one as well.
                if (rank[p_A] < rank[p_B]) { parent[p_A] = p_B; }
                else if (rank[p_A] > rank[p_B]) { parent[p_B] = p_A; }
                else {
                    parent[p_A] = p_B;
                    ++rank[p_B];
                }
            }
        };

        // Find and unite parent IDs using an union find.
        UnionFind union_find;
        for (auto& [cluster_id, cluster_data]: rClustersMap) {
            const auto& found_ids = std::get<2>(cluster_data);
            for (std::size_t connected_id : found_ids) {
                if (rClustersMap.find(connected_id) != rClustersMap.end()) {
                    union_find.unite_parents(cluster_id, connected_id);
                }
            }
        }

        // Merge clusters with the same parents.
        ClustersMapType merged_map;
        for (auto& [cluster_id, cluster_data]: rClustersMap) {
            // Get cluster parent and boundary nodes.
            const std::size_t parent_id = union_find.find_parent(cluster_id);
            NodesSetType& boundary_nodes = std::get<0>(cluster_data);

            // Merge boundary nodes into parent (created on first call).
            auto& merged_cluster_data =  merged_map[parent_id];
            NodesSetType& merged_boundary_nodes = std::get<0>(merged_cluster_data);
            merged_boundary_nodes.insert(boundary_nodes.begin(), boundary_nodes.end());
        }

        // Remove clusters that consist of one single node and set their nodes cluster ID to 1.
        for (auto it = merged_map.begin(); it != merged_map.end(); ) {
            NodesSetType& r_boundary_nodes = std::get<0>((*it).second);
            const std::size_t n_boundary_nodes = r_boundary_nodes.size();
            if (n_boundary_nodes < 2) {
                for (auto p_node : r_boundary_nodes) {
                    p_node->FastGetSolutionStepValue(ACTIVATION_LEVEL) = 1;
                }
                it = merged_map.erase(it);
            } else {
                ++it;
            }
        }

        // Assign new cluster IDs in ascending order starting from 1 for easier identification
        // and adapt the ACTIVATION_LEVEL of all cluster nodes accordingly.
        ClustersMapType new_ids_map;
        std::size_t new_id = 1;
        for (auto& [merged_id, merged_data]: merged_map) {
            NodesSetType& merged_boundary_nodes = std::get<0>(merged_data);

            // Create cluster with new ID and move boundary nodes to new cluster
            auto& new_id_data =  new_ids_map[++new_id];
            NodesSetType& new_boundary_nodes = std::get<0>(new_id_data);
            new_boundary_nodes = std::move(merged_boundary_nodes);

            // Change the nodal variable for the cluster id.
            std::for_each(new_boundary_nodes.begin(), new_boundary_nodes.end(), [&new_id](NodeType::Pointer pNode){
                pNode->FastGetSolutionStepValue(ACTIVATION_LEVEL) = new_id;
            });
        }

        // Replace old map with merged and cleansed one.
        rClustersMap.swap(new_ids_map);
    }

    const bool ShiftedBoundaryPointBasedInterfaceUtility::FindClusterElementsUntilFixedDof(
        const std::size_t ClusterId,
        const NodesSetType& rBoundaryNodes,
        ElementsSetType& rClusterElements)
    {
        NodesSetType current_frontier;

        // Add all boundary nodes and check them for a fixed dof.
        for (auto p_node : rBoundaryNodes) {
            current_frontier.insert(p_node);
            const Node::DofsContainerType& nodal_dofs = p_node->GetDofs();
            for(auto it_dof = nodal_dofs.begin() ; it_dof != nodal_dofs.end() ; it_dof++){
                if((*it_dof)->IsFixed()){
                    //fixed_dof_found = true;
                    return true;
                }
            }
        }

        //TODO use parallelization?
        // #pragma omp parallel for schedule(dynamic) shared(fixed_dof_found)
        // for (int k = 1; k <= max_cluster_id; ++k) {
        //     if(fixed_dof_found) continue;

        // Search for new elements of the cluster while there are nodes in the frontier.
        while (!current_frontier.empty()) {
            // Get a node from the frontier and erase it from the set
            auto it_node = current_frontier.begin();
            NodeType::Pointer p_node = *it_node;
            current_frontier.erase(*it_node);  //TODO can I erase the entry without disturbing the pointer?

            // Check the frontier node's neighboring elements.
            auto& r_elem_neigh_vect = p_node->GetValue(NEIGHBOUR_ELEMENTS);
            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                if (p_elem_neigh != nullptr) {
                    if (p_elem_neigh->Is(ACTIVE)) {
                        // Add neighboring elements if it is not part of the boundary.
                        rClusterElements.insert(p_elem_neigh);
                        for (auto& r_neigh_node : p_elem_neigh->GetGeometry()) {

                            // Check cluster ID of the node and add it to the cluster and the frontier if it is not part of a cluster yet.
                            auto& r_node_cluster_id = r_neigh_node.FastGetSolutionStepValue(ACTIVATION_LEVEL);
                            if (r_node_cluster_id == 0) {
                                r_node_cluster_id = ClusterId;
                                current_frontier.insert(&r_neigh_node);

                                // Check for new node if a degree of freedom is fixed
                                const Node::DofsContainerType& nodal_dofs = (&r_neigh_node)->GetDofs();
                                for(auto it_dof = nodal_dofs.begin() ; it_dof != nodal_dofs.end() ; it_dof++){
                                    if((*it_dof)->IsFixed()){
                                        //fixed_dof_found = true;
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Return that no fixed dof was found, when all cluster elements were found successfully without early return.
        return false;
    }


    void ShiftedBoundaryPointBasedInterfaceUtility::SetSidesVectorsAndSkinNormalsForSplitElements(
        const SkinPointsToElementsMapType& rSkinPointsMap,
        SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap)
    {
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

    void ShiftedBoundaryPointBasedInterfaceUtility::SetExtensionOperatorsForSplitElementNodes(
        const SidesVectorToElementsMapType& rSidesVectorMap,
        AverageSkinToElementsMapType& rAvgSkinMap,
        NodesCloudMapType& rExtensionOperatorMap)
    {
        // Get the extension operator shape functions function
        auto p_meshless_sh_func = mExtensionOperator == ExtensionOperator::MLS ? GetMLSShapeFunctionsFunction() : GetRBFShapeFunctionsFunction();

        // Get support node clouds for all nodes of all split elements and calculate their extension operators
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
                            cloud_data_vector(i_cl_nod) = i_data;
                        }
                        const auto ext_op_key_data = std::make_pair(p_node, cloud_data_vector);
                        //TODO make this threadsafe for parallelization
                        rExtensionOperatorMap.insert(ext_op_key_data);
                    // } else {
                    //     KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility")
                    //     << "No enough support nodes were found for node " << p_node->Id() << ". Extension basis can not be calculated." << std::endl;
                    }
                }
            }
        }
    }

    void ShiftedBoundaryPointBasedInterfaceUtility::SetLateralSupportCloud(
        const NodeType::Pointer pOtherSideNode,
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates,
        const Kratos::Flags& rSearchSideFlag)
    {
        // Find the support cloud of nodes on the search side (other side as the given node)
        // NOTE that we use an unordered_set to ensure that these are unique
        // NOTE that we check the order of the MLS interpolation to add nodes from enough layers
        NodesSetType aux_set;
        std::vector<NodeType::Pointer> cur_layer_nodes;
        std::vector<NodeType::Pointer> prev_layer_nodes;
        const std::size_t n_layers = mMLSExtensionOperatorOrder + 1;

        // Find elemental neighbors of the given node and add their nodes to the cloud nodes set if they are located on the search side
        // This is the first layer of sampling/ support points
        // NOTE that the sides of the first layer of nodes at gamma_tilde already need to be defined in their flags (done by SetSidesVectorsAndSkinNormalsForSplitElements)
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
            // KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility")
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
        //    KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << n_extra_layers << " extra layers of points needed for MLS calculation." << std::endl;
        // }

        // Add obtained cloud nodes to the cloud node vector and sort them by id
        //TODO sorting really necessary or helpful??
        rCloudNodes.resize(n_cloud_nodes);
        std::size_t aux_i = 0;
        for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
            rCloudNodes(aux_i++) = *it_set;
            // Mark support node for visualization - TODO make threadsafe/ put somewhere else for parallelization
            (*it_set)->Set(rSearchSideFlag, true);
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
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet)
    {
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

    void ShiftedBoundaryPointBasedInterfaceUtility::AddLateralSupportLayer(
        const array_1d<double,3>& rAvgSkinPosition,
        const array_1d<double,3>& rAvgSkinNormal,
        const std::vector<NodeType::Pointer>& PreviousLayerNodes,
        std::vector<NodeType::Pointer>& CurrentLayerNodes,
        NodesSetType& SupportNodesSet)
    {
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

    void ShiftedBoundaryPointBasedInterfaceUtility::CreateCloudNodeVectorsForSplitElement(
        const ElementType& rElement,
        const Vector& rSidesVector,
        NodesCloudMapType& rExtensionOperatorMap,
        PointerVector<NodeType>& rCloudNodeVectorPositiveSide,
        PointerVector<NodeType>& rCloudNodeVectorNegativeSide)
    {
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

    void ShiftedBoundaryPointBasedInterfaceUtility::GetDataForSplitElementSkinPoint(
        const ElementType& rElement,
        const array_1d<double,3>& rIntPtCoordinates,
        Vector& rIntPtShapeFunctionValues,
        Matrix& rIntPtShapeFunctionDerivatives)
    {
        const auto& r_geom = rElement.GetGeometry();

        // Compute the local coordinates of the integration point in the element's geometry
        array_1d<double,3> int_pt_local_coords = ZeroVector(3);
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
        std::size_t& r_ConditionId,
        const bool ConsiderPositiveSide)
    {
        const auto& r_geom = rElement.GetGeometry();

        // Initialize the extension operator containers
        const std::size_t n_cl_nodes = rCloudNodeVector.size();
        const std::size_t n_dim = r_geom.WorkingSpaceDimension();
        Vector N_container = ZeroVector(n_cl_nodes);
        Matrix DN_DX_container = ZeroMatrix(n_cl_nodes, n_dim);

        array_1d<double,3> area_normal = rIntPtAreaNormal;
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
                    //     KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "No wall condition will be created for positive side of the skin point because Node No." << p_node->Id() << " is not active." << std::endl;
                    // } else {
                    //     KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "No wall condition will be created for negative side of the skin point because Node No." << p_node->Id() << " is not active." << std::endl;
                    // }
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
                    //KRATOS_WARNING("ShiftedBoundaryPointBasedInterfaceUtility") << "No wall condition will be created for one side of the skin point because no extension operator was available for Node No." << p_node->Id() << std::endl;
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
            N_container_2(i_node) = rIntPtShapeFunctionValues(i_node);
            for (std::size_t d = 0; d < n_dim; ++d) {
                DN_DX_container_2(i_node, d) += rIntPtShapeFunctionDerivatives(i_node, d);
            }
        }*/

        // Create a new condition with a geometry made up with the basis nodes
        auto p_prop = rElement.pGetProperties();
        auto p_cond = mpConditionPrototype->Create(++r_ConditionId, rCloudNodeVector, p_prop);
        p_cond->Set(ACTIVE, true);
        mpBoundarySubModelPart->AddCondition(p_cond);

        //TODO for laplacian static heat testing:
        // Set Dirichlet boundary condition
        //const double dirichlet_value = std::pow(rIntPtCoordinates(0),2) + std::pow(rIntPtCoordinates(1),2);
        //p_cond->SetValue(TEMPERATURE, dirichlet_value);

        // Store the SBM BC data in the condition database
        p_cond->SetValue(ELEMENT_H, ElementSize);
        p_cond->SetValue(INTEGRATION_COORDINATES, rIntPtCoordinates);
        p_cond->SetValue(NORMAL, area_normal/ skin_pt_weight);
        p_cond->SetValue(INTEGRATION_WEIGHT, skin_pt_weight);
        p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
        p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);

        return true;
    }

    bool ShiftedBoundaryPointBasedInterfaceUtility::SetEnclosedNodesPressure(
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
                // Calculate dot product of average skin normal of the element and the normalized vector between averaged skin point and the element's node
                // array_1d<double,3> avg_skin_pt_to_node = r_node.Coordinates() - rAvgSkinPosition;
                // avg_skin_pt_to_node /= norm_2(avg_skin_pt_to_node);
                // double dot_product = inner_prod(avg_skin_pt_to_node, rAvgSkinNormal);
                // if (mNegativeSideIsEnclosed) {
                //     dot_product *= -1.0;
                // }
                // if (dot_product > 0.1 &&
                if (r_node.Is(ACTIVE)) {
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

    ShiftedBoundaryPointBasedInterfaceUtility::MeshlessShapeFunctionsFunctionType ShiftedBoundaryPointBasedInterfaceUtility::GetRBFShapeFunctionsFunction() const
    {
        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
            RBFShapeFunctionsUtility::CalculateShapeFunctions(rPoints, rX, h, rN);
        };
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

    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPointsTemplated<2>();
    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPointsTemplated<3>();

    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPointsAndNodesTemplated<2>();
    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::CalculateVariablesAtSkinPointsAndNodesTemplated<3>();

    template bool KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::CalculateUnknownsForBothSidesOfSplitElement<2>(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns);
    template bool KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::CalculateUnknownsForBothSidesOfSplitElement<3>(
        const ElementType::Pointer pElement,
        Vector& rPositiveSideUnknowns,
        Vector& rNegativeSideUnknowns);

    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::MapSkinPointsToElements<2>(SkinPointsToElementsMapType& rSkinPointsMap);
    template void KRATOS_API(KRATOS_CORE) ShiftedBoundaryPointBasedInterfaceUtility::MapSkinPointsToElements<3>(SkinPointsToElementsMapType& rSkinPointsMap);

}  // namespace Kratos.
