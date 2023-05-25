//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/find_intersected_geometrical_objects_with_obb_process.h"

namespace Kratos
{
/// Local Flags
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedGeometricalObjectsWithOBBProcess, DEBUG_OBB, 4);
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedGeometricalObjectsWithOBBProcess, SEPARATING_AXIS_THEOREM, 5);
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedGeometricalObjectsWithOBBProcess, BUILD_OBB_FROM_BB, 6);

/***********************************************************************************/
/***********************************************************************************/

FindIntersectedGeometricalObjectsWithOBBProcess::FindIntersectedGeometricalObjectsWithOBBProcess(
    ModelPart& rModelPartIntersected,
    ModelPart& rModelPartIntersecting,
    const double BoundingBoxFactor,
    const Flags Options
    ) : BaseType(rModelPartIntersected, rModelPartIntersecting, Options),
        mBoundingBoxFactor(BoundingBoxFactor)
{
    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);

#ifdef  KRATOS_DEBUG
    // We create new properties for debugging
    if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
        rModelPartIntersected.CreateNewProperties(10001);
        rModelPartIntersected.CreateSubModelPart(rModelPartIntersected.Name() + "_AUXILIAR_DEBUG_OBB");
        rModelPartIntersecting.CreateNewProperties(10002);
        rModelPartIntersecting.CreateSubModelPart(rModelPartIntersecting.Name() + "_AUXILIAR_DEBUG_OBB");
    }
#endif
}

/***********************************************************************************/
/***********************************************************************************/

FindIntersectedGeometricalObjectsWithOBBProcess::FindIntersectedGeometricalObjectsWithOBBProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : BaseType(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString()),
        rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    mThisParameters = ThisParameters;

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_intersected_model_part_name = mThisParameters["intersected_model_part_name"].GetString();
    const std::string& r_intersecting_model_part_name = mThisParameters["intersecting_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_intersected_model_part_name == "") << "intersected_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_intersecting_model_part_name == "") << "intersecting_model_part_name must be defined on parameters" << std::endl;

    // Setting flags
    const bool intersecting_conditions = ThisParameters["intersecting_conditions"].GetBool();
    const bool intersecting_elements = ThisParameters["intersecting_elements"].GetBool();
    const bool intersected_conditions = ThisParameters["intersected_conditions"].GetBool();
    const bool intersected_elements = ThisParameters["intersected_elements"].GetBool();

    BaseType::mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_CONDITIONS, intersecting_conditions);
    BaseType::mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_ELEMENTS, intersecting_elements);
    BaseType::mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS, intersected_conditions);
    BaseType::mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS, intersected_elements);

    // Setting the bounding box factor
    mBoundingBoxFactor = mThisParameters["bounding_box_factor"].GetDouble();

    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);

    // If we debug OBB
    BaseType::mOptions.Set(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB, mThisParameters["debug_obb"].GetBool());

    // If we build the OBB from the geometry BB
    BaseType::mOptions.Set(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB, mThisParameters["build_from_bounding_box"].GetBool());

    // The intersection type
    ConvertIntersection(mThisParameters["OBB_intersection_type"].GetString());

    // We create new properties for debugging
#ifdef  KRATOS_DEBUG
    if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
        this->GetModelPart1().CreateNewProperties(1001);
        this->GetModelPart1().CreateSubModelPart(this->GetModelPart1().Name() + "_AUXILIAR_DEBUG_OBB");
        this->GetModelPart2().CreateNewProperties(1002);
        this->GetModelPart2().CreateSubModelPart(this->GetModelPart2().Name() + "_AUXILIAR_DEBUG_OBB");
    }
#endif
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBProcess::SetOctreeBoundingBox()
{
    // Getting first iterators
    const auto it_node_begin_1 = BaseType::mrModelPartIntersected.NodesBegin();
    const auto it_node_begin_2 = BaseType::mrModelPartIntersecting.NodesBegin();

    // Setting initial guess
    PointType low(it_node_begin_1->Coordinates());
    PointType high(it_node_begin_1->Coordinates());

    // Loop over all nodes in first modelpart
    for (IndexType i_node = 0 ; i_node < BaseType::mrModelPartIntersected.NumberOfNodes(); ++i_node) {
        auto it_node = it_node_begin_1 + i_node;
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Loop over all skin nodes
    for (IndexType i_node = 0 ; i_node < BaseType::mrModelPartIntersecting.NumberOfNodes(); ++i_node) {
        auto it_node = it_node_begin_2 + i_node;
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Slightly increase the bounding box size to avoid problems with geometries in the borders
    // Note that std::numeric_limits<double>::double() is added for the 2D cases. Otherwise, the
    // third component will be 0, breaking the octree behaviour.
    for(IndexType i = 0 ; i < 3; i++) {
        low[i] -= std::abs(high[i] - low[i])*1e-3 + std::numeric_limits<double>::epsilon();
        high[i] += std::abs(high[i] - low[i])*1e-3 + std::numeric_limits<double>::epsilon();
    }

    // In case we consider the bounding box we extend box
    if (mBoundingBoxFactor > 0.0) {
        const std::size_t dimension = BaseType::mrModelPartIntersected.GetProcessInfo()[DOMAIN_SIZE];
        for(IndexType i = 0 ; i < dimension; i++) {
            if (std::abs(high[i] - low[i]) < mBoundingBoxFactor) {
                low[i] -= 1.0/3.0 * mBoundingBoxFactor;
                high[i] += 2.0/3.0 * mBoundingBoxFactor;
            }
        }
    }

    // TODO: Octree needs refactoring to work with BoundingBox. Pooyan.
    GetOctreePointer()->SetBoundingBox(low.data().data(), high.data().data());
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBProcess::HasIntersection(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    const IndexType working_space_dimension = rFirstGeometry.WorkingSpaceDimension(); // TODO: DOMAIN_SIZE should be considered for consistency with other implementations
    const IndexType local_space_dimension = rFirstGeometry.LocalSpaceDimension();
    if (this->Is(BOUNDARY)) {
        if (working_space_dimension == 2) {
            if (local_space_dimension == 2) {
                return this->HasIntersection2D(rFirstGeometry, rSecondGeometry);
            } else {
                return this->HasDirectIntersection2D(rFirstGeometry, rSecondGeometry);
            }
        } else {
            if (local_space_dimension == 3) {
                return this->HasIntersection3D(rFirstGeometry, rSecondGeometry);
            } else {
                return this->HasDirectIntersection3D(rFirstGeometry, rSecondGeometry);
            }
        }
    } else {
        return rFirstGeometry.HasIntersection(rSecondGeometry);
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBProcess::HasIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // The edges
    const PointerVector<GeometryType> r_edges_1 = rFirstGeometry.GenerateEdges();
    const std::size_t number_of_edges_1 = r_edges_1.size();
    const PointerVector<GeometryType> r_edges_2 = rSecondGeometry.GenerateEdges();
    const std::size_t number_of_edges_2 = r_edges_2.size();

    // First geometry
    for (std::size_t i_1 = 0; i_1 < number_of_edges_1; ++i_1) {
        auto& r_edge_1 = *(r_edges_1.begin() + i_1);

        OrientedBoundingBox<2> first_obb(r_edge_1, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

    #ifdef  KRATOS_DEBUG
        // We create new elements for debugging
        if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
            auto p_prop = this->GetModelPart1().pGetProperties(1001);
            CreateDebugOBB2D(this->GetModelPart1(), p_prop, first_obb);
        }
    #endif

        // Second geometry
        for (std::size_t i_2 = 0; i_2 < number_of_edges_2; ++i_2) {
            auto& r_edge_2 = *(r_edges_2.begin() + i_2);

            OrientedBoundingBox<2> second_obb(r_edge_2, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

        #ifdef  KRATOS_DEBUG
            // We create new elements for debugging
            if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
                auto p_prop = this->GetModelPart2().pGetProperties(1002);
                CreateDebugOBB2D(this->GetModelPart2(), p_prop, second_obb);
            }
        #endif

            // Computing intersection OBB
            if (first_obb.HasIntersection(second_obb, GetOBBHasIntersectionType())){
                return true;
            }
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBProcess::HasDirectIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // First geometry
    OrientedBoundingBox<2> first_obb(rFirstGeometry, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

#ifdef  KRATOS_DEBUG
    // We create new elements for debugging
    if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
        auto p_prop = this->GetModelPart1().pGetProperties(1001);
        CreateDebugOBB2D(this->GetModelPart1(), p_prop, first_obb);
    }
#endif

    // Second geometry
    OrientedBoundingBox<2> second_obb(rSecondGeometry, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

#ifdef  KRATOS_DEBUG
    // We create new elements for debugging
    if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
        auto p_prop = this->GetModelPart2().pGetProperties(1002);
        CreateDebugOBB2D(this->GetModelPart2(), p_prop, second_obb);
    }
#endif

    // Computing intersection OBB
    if (first_obb.HasIntersection(second_obb, GetOBBHasIntersectionType())){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBProcess::HasIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // The faces
    const PointerVector<GeometryType> r_faces_1 = rFirstGeometry.GenerateFaces();
    const std::size_t number_of_faces_1 = r_faces_1.size();

    const PointerVector<GeometryType> r_faces_2 = rSecondGeometry.GenerateFaces();
    const std::size_t number_of_faces_2 = r_faces_2.size();

    // First geometry
    for (std::size_t i_1 = 0; i_1 < number_of_faces_1; ++i_1) {
        auto& r_face_1 = *(r_faces_1.begin() + i_1);

        // Creating OBB
        OrientedBoundingBox<3> first_obb(r_face_1, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

    #ifdef  KRATOS_DEBUG
        // We create new elements for debugging
        if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
            auto p_prop = this->GetModelPart1().pGetProperties(1001);
            CreateDebugOBB3D(this->GetModelPart1(), p_prop, first_obb);
        }
    #endif

        // Second geometry
        for (std::size_t i_2 = 0; i_2 < number_of_faces_2; ++i_2) {
            auto& r_face_2 = *(r_faces_2.begin() + i_2);

            OrientedBoundingBox<3> second_obb(r_face_2, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

        #ifdef  KRATOS_DEBUG
            // We create new elements for debugging
            if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
                auto p_prop = this->GetModelPart2().pGetProperties(1002);
                CreateDebugOBB3D(this->GetModelPart2(), p_prop, second_obb);
            }
        #endif

            // Computing intersection OBB
            if (first_obb.HasIntersection(second_obb, GetOBBHasIntersectionType())){
                return true;
            }
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBProcess::HasDirectIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // First geometry
    OrientedBoundingBox<3> first_obb(rFirstGeometry, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

#ifdef  KRATOS_DEBUG
    // We create new elements for debugging
    if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
        auto p_prop = this->GetModelPart1().pGetProperties(1001);
        CreateDebugOBB3D(this->GetModelPart1(), p_prop, first_obb);
    }
#endif

    // Second geometry
    OrientedBoundingBox<3> second_obb(rSecondGeometry, mBoundingBoxFactor, BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::BUILD_OBB_FROM_BB));

#ifdef  KRATOS_DEBUG
    // We create new elements for debugging
    if (BaseType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
        auto p_prop = this->GetModelPart2().pGetProperties(1002);
        CreateDebugOBB3D(this->GetModelPart2(), p_prop, second_obb);
    }
#endif

    // Computing intersection OBB
    if (first_obb.HasIntersection(second_obb, GetOBBHasIntersectionType())){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBProcess::CreateDebugOBB2D(
    ModelPart& rModelPart,
    Properties::Pointer pProperties,
    OrientedBoundingBox<2>& rOrientedBoundingBox
    )
{
    ModelPart& r_sub_model_part = rModelPart.GetSubModelPart(rModelPart.Name() + "_AUXILIAR_DEBUG_OBB");

    const std::size_t initial_node_id = rModelPart.GetRootModelPart().NumberOfNodes();// NOTE: We assume ordered nodes
    const auto quad = rOrientedBoundingBox.GetEquivalentGeometry();
    const array_1d<double, 3>& r_center = rOrientedBoundingBox.GetCenter();
    const array_1d<array_1d<double, 3>, 2>& r_orientation_vectors = rOrientedBoundingBox.GetOrientationVectors();
    std::vector<NodeType::Pointer> element_nodes (4);
    for (int i = 0; i < 4; ++i) {
        element_nodes[i] = r_sub_model_part.CreateNewNode(initial_node_id + i + 1, quad[i].X(), quad[i].Y(), quad[i].Z());
    }
    auto p_node_center = r_sub_model_part.CreateNewNode(initial_node_id + 4 + 1, r_center[0], r_center[1], r_center[2]);
    p_node_center->SetValue(NORMAL, r_orientation_vectors[0]);
    p_node_center->SetValue(TANGENT_XI, r_orientation_vectors[1]);

    const std::size_t initial_element_id = rModelPart.GetRootModelPart().NumberOfElements();// NOTE: We assume ordered elements
    r_sub_model_part.CreateNewElement("Element2D4N", initial_element_id + 1, PointerVector<NodeType>{element_nodes}, pProperties);

//     // More debug
//     const array_1d<double, 2>& r_half_lenghts = rOrientedBoundingBox.GetHalfLength();
//     std::cout << "array_1d<double, 3> center;\ncenter[0] = " << r_center[0] << ";\ncenter[1] = " << r_center[1] << ";\ncenter[2] = " << r_center[2] << ";\n\narray_1d<array_1d<double, 3>, 2> directions;\ndirections[0][0] = " << r_orientation_vectors[0][0] << ";\ndirections[0][1] = " << r_orientation_vectors[0][1] << ";\ndirections[0][2] = " << r_orientation_vectors[0][2] << ";\ndirections[1][0] = " << r_orientation_vectors[1][0] << ";\ndirections[1][1] = " << r_orientation_vectors[1][1] << ";\ndirections[1][2] = " << r_orientation_vectors[1][2] << ";\n\narray_1d<double, 2> half_lenghts;\nhalf_lenghts[0] = " << r_half_lenghts[0] << ";\nhalf_lenghts[1] = " << r_half_lenghts[1] << ";\nOrientedBoundingBox<2> obb(center, directions, half_lenghts);" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBProcess::CreateDebugOBB3D(
    ModelPart& rModelPart,
    Properties::Pointer pProperties,
    OrientedBoundingBox<3>& rOrientedBoundingBox
    )
{
    ModelPart& r_sub_model_part = rModelPart.GetSubModelPart(rModelPart.Name() + "_AUXILIAR_DEBUG_OBB");

    const std::size_t initial_node_id = rModelPart.GetRootModelPart().NumberOfNodes();// NOTE: We assume ordered nodes
    const auto hexa = rOrientedBoundingBox.GetEquivalentGeometry();
    const array_1d<double, 3>& r_center = rOrientedBoundingBox.GetCenter();
    const array_1d<array_1d<double, 3>, 3>& r_orientation_vectors = rOrientedBoundingBox.GetOrientationVectors();
    std::vector<NodeType::Pointer> element_nodes (8);
    for (int i = 0; i < 8; ++i) {
        element_nodes[i] = r_sub_model_part.CreateNewNode(initial_node_id + i + 1, hexa[i].X(), hexa[i].Y(), hexa[i].Z());
    }
    auto p_node_center = r_sub_model_part.CreateNewNode(initial_node_id + 8 + 1, r_center[0], r_center[1], r_center[2]);
    p_node_center->SetValue(NORMAL, r_orientation_vectors[0]);
    p_node_center->SetValue(TANGENT_XI, r_orientation_vectors[1]);
    p_node_center->SetValue(TANGENT_ETA, r_orientation_vectors[2]);

    const std::size_t initial_element_id = rModelPart.GetRootModelPart().NumberOfElements();// NOTE: We assume ordered elements
    r_sub_model_part.CreateNewElement("Element3D8N", initial_element_id + 1, PointerVector<NodeType>{element_nodes}, pProperties);

//     // More debug
//     const array_1d<double, 3>& r_half_lenghts = rOrientedBoundingBox.GetHalfLength();
//     std::cout << "array_1d<double, 3> center;\ncenter[0] = " << r_center[0] << ";\ncenter[1] = " << r_center[1] << ";\ncenter[2] = " << r_center[2] << ";\n\narray_1d<array_1d<double, 3>, 3> directions;\ndirections[0][0] = " << r_orientation_vectors[0][0] << ";\ndirections[0][1] = " << r_orientation_vectors[0][1] << ";\ndirections[0][2] = " << r_orientation_vectors[0][2] << ";\ndirections[1][0] = " << r_orientation_vectors[1][0] << ";\ndirections[1][1] = " << r_orientation_vectors[1][1] << ";\ndirections[1][2] = " << r_orientation_vectors[1][2] << ";\ndirections[2][0] = " << r_orientation_vectors[2][0] << ";\ndirections[2][1] = " << r_orientation_vectors[2][1] << ";\ndirections[2][2] = " << r_orientation_vectors[2][2] << ";\n\narray_1d<double, 3> half_lenghts;\nhalf_lenghts[0] = " << r_half_lenghts[0] << ";\nhalf_lenghts[1] = " << r_half_lenghts[1] << ";\nhalf_lenghts[2] = " << r_half_lenghts[2] <<  ";\nOrientedBoundingBox<3> obb(center, directions, half_lenghts);" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters FindIntersectedGeometricalObjectsWithOBBProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "intersected_model_part_name"  : "",
        "intersecting_model_part_name" : "",
        "bounding_box_factor"          : -1.0,
        "debug_obb"                    : false,
        "OBB_intersection_type"        : "SeparatingAxisTheorem",
        "build_from_bounding_box"      : true,
        "intersecting_conditions"      : true,
        "intersecting_elements"        : true,
        "intersected_conditions"       : true,
        "intersected_elements"         : true
    })" );

    return default_parameters;
}

}  // namespace Kratos.
