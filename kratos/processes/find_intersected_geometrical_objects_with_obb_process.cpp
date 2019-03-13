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
template<class TEntity>
FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::FindIntersectedGeometricalObjectsWithOBBProcess(
    ModelPart& rPart1,
    ModelPart& rPart2,
    const double BoundingBoxFactor,
    const bool DebugOBB
    ) : BaseType(rPart1, rPart2),
        mBoundingBoxFactor(BoundingBoxFactor),
        mDebugOBB(DebugOBB)
{
    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);

    // We create new properties for debugging
    if (mDebugOBB) {
        rPart1.CreateNewProperties(10001);
        rPart1.CreateSubModelPart(rPart1.Name() + "_AUXILIAR_DEBUG_OBB");
        rPart2.CreateNewProperties(10002);
        rPart2.CreateSubModelPart(rPart2.Name() + "_AUXILIAR_DEBUG_OBB");
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::FindIntersectedGeometricalObjectsWithOBBProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : BaseType(rModel.GetModelPart(ThisParameters["first_model_part_name"].GetString()),
        rModel.GetModelPart(ThisParameters["second_model_part_name"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    BaseType::mThisParameters = ThisParameters;

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_first_model_part_name = BaseType::mThisParameters["first_model_part_name"].GetString();
    const std::string& r_second_model_part_name = BaseType::mThisParameters["second_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_first_model_part_name == "") << "first_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_second_model_part_name == "") << "second_model_part_name must be defined on parameters" << std::endl;

    // Setting the bounding box factor
    mBoundingBoxFactor = BaseType::mThisParameters["bounding_box_factor"].GetDouble();

    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);

    // If we debug OBB
    mDebugOBB = BaseType::mThisParameters["debug_obb"].GetBool();

    // We create new properties for debugging
    if (mDebugOBB) {
        this->GetModelPart1().CreateNewProperties(1001);
        this->GetModelPart1().CreateSubModelPart(this->GetModelPart1().Name() + "_AUXILIAR_DEBUG_OBB");
        this->GetModelPart2().CreateNewProperties(1002);
        this->GetModelPart2().CreateSubModelPart(this->GetModelPart2().Name() + "_AUXILIAR_DEBUG_OBB");
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::HasIntersection(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    const std::size_t work_dim = rFirstGeometry.WorkingSpaceDimension(); // TODO: DOMAIN_SIZE should be considered for consistency with other implementations
    if (this->Is(BOUNDARY)) {
        if (work_dim == 2) {
            return this->HasIntersection2D(rFirstGeometry, rSecondGeometry);
        } else {
            return this->HasIntersection3D(rFirstGeometry, rSecondGeometry);
        }
    } else {
        if (work_dim == 2) {
            return BaseType::HasIntersection2D(rFirstGeometry, rSecondGeometry);
        } else {
            return BaseType::HasIntersection3D(rFirstGeometry, rSecondGeometry);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::HasIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // The local dimensions
    const std::size_t local_dimension_1 = rFirstGeometry.LocalSpaceDimension();
    const std::size_t local_dimension_2 = rSecondGeometry.LocalSpaceDimension();

    // The edges
    PointerVector<GeometryType> r_edges_1 = rFirstGeometry.Edges();
    const std::size_t number_of_edges_1 = (local_dimension_1 < 2) ? 1 : r_edges_1.size();
    PointerVector<GeometryType> r_edges_2 = rSecondGeometry.Edges();
    const std::size_t number_of_edges_2 = (local_dimension_2 < 2) ? 1 : r_edges_2.size();

    // First geometry
    for (std::size_t i_1 = 0; i_1 < number_of_edges_1; ++i_1) {
        auto& r_edge_1 = *(r_edges_1.begin() + i_1);

        // Check the intersection of each edge of the object bounding box against the intersecting object bounding box
        NodeType first_geometry_low_node, first_geometry_high_node;
        r_edge_1.BoundingBox(first_geometry_low_node, first_geometry_high_node);

        // Creating OBB
        array_1d<array_1d<double, 3>, 2> first_direction_vector;
        noalias(first_direction_vector[0]) = first_geometry_high_node - first_geometry_low_node;
        const double norm_first_direction_vector = norm_2(first_direction_vector[0]);
        if (norm_first_direction_vector > std::numeric_limits<double>::epsilon())
            first_direction_vector[0] /= norm_first_direction_vector;
        else
            KRATOS_ERROR << "Zero norm on OrientedBoundingBox direction" << std::endl;
        first_direction_vector[1][0] = first_direction_vector[0][1];
        first_direction_vector[1][1] = - first_direction_vector[0][0];
        first_direction_vector[1][2] = 0.0;
        const array_1d<double, 3> first_center_point = 0.5 * (first_geometry_low_node.Coordinates() + first_geometry_high_node.Coordinates());
        array_1d<double, 2> first_half_distances;
        first_half_distances[0] = 0.5 * norm_first_direction_vector + mBoundingBoxFactor;
        first_half_distances[1] = mBoundingBoxFactor;
        OrientedBoundingBox<2> first_obb(first_center_point, first_direction_vector, first_half_distances);

        // We create new elements for debugging
        if (mDebugOBB) {
            auto p_prop = this->GetModelPart1().pGetProperties(1001);
            CreateDebugOBB2D(this->GetModelPart1(), p_prop, first_obb);
        }

        // Second geometry
        for (std::size_t i_2 = 0; i_2 < number_of_edges_2; ++i_2) {
            auto& r_edge_2 = *(r_edges_2.begin() + i_2);

            // Check the intersection of each edge of the object bounding box against the intersecting object bounding box
            NodeType second_geometry_low_node, second_geometry_high_node;
            r_edge_2.BoundingBox(second_geometry_low_node, second_geometry_high_node);

            // Creating OBB
            array_1d<array_1d<double, 3>, 2> second_direction_vector;
            noalias(second_direction_vector[0]) = second_geometry_high_node - second_geometry_low_node;
            const double norm_second_direction_vector = norm_2(second_direction_vector[0]);
            second_direction_vector[0] /= norm_second_direction_vector;
            second_direction_vector[1][0] = second_direction_vector[0][1];
            second_direction_vector[1][1] = - second_direction_vector[0][0];
            second_direction_vector[1][2] = 0.0;
            const array_1d<double, 3> second_center_point = 0.5 * (second_geometry_low_node.Coordinates() + second_geometry_high_node.Coordinates());
            array_1d<double, 2> second_half_distances;
            second_half_distances[0] = 0.5 * norm_second_direction_vector + mBoundingBoxFactor;
            second_half_distances[1] = mBoundingBoxFactor;
            OrientedBoundingBox<2> second_obb(second_center_point, second_direction_vector, second_half_distances);

            // We create new elements for debugging
            if (mDebugOBB) {
                auto p_prop = this->GetModelPart2().pGetProperties(1002);
                CreateDebugOBB2D(this->GetModelPart2(), p_prop, second_obb);
            }

            // Computing intersection OBB
            if (first_obb.HasIntersection(second_obb)){
                return true;
            }
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::HasIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // The local dimensions
    const std::size_t local_dimension_1 = rFirstGeometry.LocalSpaceDimension();
    const std::size_t local_dimension_2 = rSecondGeometry.LocalSpaceDimension();

    // The faces
    PointerVector<GeometryType> r_faces_1 = rFirstGeometry.Faces();
    const std::size_t number_of_faces_1 = (local_dimension_1 < 3) ? 1 : r_faces_1.size();

    PointerVector<GeometryType> r_faces_2 = rSecondGeometry.Faces();
    const std::size_t number_of_faces_2 = (local_dimension_2 < 3) ? 1 : r_faces_2.size();

    // First geometry
    for (std::size_t i_1 = 0; i_1 < number_of_faces_1; ++i_1) {
        auto& r_face_1 = *(r_faces_1.begin() + i_1);

        // Check the intersection of each face of the object bounding box against the intersecting object bounding box
        NodeType first_geometry_low_node, first_geometry_high_node;
        r_face_1.BoundingBox(first_geometry_low_node, first_geometry_high_node);

        // Creating OBB
        array_1d<array_1d<double, 3>, 3> first_direction_vector;
        noalias(first_direction_vector[0]) = first_geometry_high_node - first_geometry_low_node;
        const double norm_first_direction_vector = norm_2(first_direction_vector[0]);
        if (norm_first_direction_vector > std::numeric_limits<double>::epsilon())
            first_direction_vector[0] /= norm_first_direction_vector;
        else
            KRATOS_ERROR << "Zero norm on OrientedBoundingBox direction" << std::endl;
        MathUtils<double>::OrthonormalBasis(first_direction_vector[0], first_direction_vector[1], first_direction_vector[2]);
        const array_1d<double, 3> first_center_point = r_face_1.Center().Coordinates();// 0.5 * (first_geometry_low_node.Coordinates() + first_geometry_high_node.Coordinates());
        array_1d<double, 3> first_half_distances;
        double distance_0 = 0.0;
        double distance_1 = 0.0;
        double distance_2 = 0.0;
        for (auto& r_node : r_face_1) {
            const array_1d<double, 3> vector_points = r_node.Coordinates() - first_center_point;
            distance_0 = std::max(distance_0, std::abs(inner_prod(vector_points, first_direction_vector[0])));
            distance_1 = std::max(distance_1, std::abs(inner_prod(vector_points, first_direction_vector[1])));
            distance_2 = std::max(distance_2, std::abs(inner_prod(vector_points, first_direction_vector[2])));
        }
        first_half_distances[0] = distance_0 + mBoundingBoxFactor;
        first_half_distances[1] = distance_1 + mBoundingBoxFactor;
        first_half_distances[2] = distance_2 + mBoundingBoxFactor;
        OrientedBoundingBox<3> first_obb(first_center_point, first_direction_vector, first_half_distances);

        // We create new elements for debugging
        if (mDebugOBB) {
            auto p_prop = this->GetModelPart1().pGetProperties(1001);
            CreateDebugOBB3D(this->GetModelPart1(), p_prop, first_obb);
        }

        // Second geometry
        for (std::size_t i_2 = 0; i_2 < number_of_faces_2; ++i_2) {
            auto& r_face_2 = *(r_faces_2.begin() + i_2);

            // Check the intersection of each face of the object bounding box against the intersecting object bounding box
            NodeType second_geometry_low_node, second_geometry_high_node;
            r_face_2.BoundingBox(second_geometry_low_node, second_geometry_high_node);

            // Creating OBB
            array_1d<array_1d<double, 3>, 3> second_direction_vector;
            noalias(second_direction_vector[0]) = second_geometry_high_node - second_geometry_low_node;
            const double norm_second_direction_vector = norm_2(second_direction_vector[0]);
            second_direction_vector[0] /= norm_second_direction_vector;
            MathUtils<double>::OrthonormalBasis(second_direction_vector[0], second_direction_vector[1], second_direction_vector[2]);
            const array_1d<double, 3> second_center_point = r_face_2.Center().Coordinates();// 0.5 * (second_geometry_low_node.Coordinates() + second_geometry_high_node.Coordinates());
            array_1d<double, 3> second_half_distances;
            distance_0 = 0.0;
            distance_1 = 0.0;
            distance_2 = 0.0;
            for (auto& r_node : r_face_2) {
                const array_1d<double, 3> vector_points = r_node.Coordinates() - second_center_point;
                distance_0 = std::max(distance_0, std::abs(inner_prod(vector_points, second_direction_vector[0])));
                distance_1 = std::max(distance_1, std::abs(inner_prod(vector_points, second_direction_vector[1])));
                distance_2 = std::max(distance_2, std::abs(inner_prod(vector_points, second_direction_vector[2])));
            }
            second_half_distances[0] = distance_0 + mBoundingBoxFactor;
            second_half_distances[1] = distance_1 + mBoundingBoxFactor;
            second_half_distances[2] = distance_2 + mBoundingBoxFactor;
            OrientedBoundingBox<3> second_obb(second_center_point, second_direction_vector, second_half_distances);

            // We create new elements for debugging
            if (mDebugOBB) {
                auto p_prop = this->GetModelPart2().pGetProperties(1002);
                CreateDebugOBB3D(this->GetModelPart2(), p_prop, second_obb);
            }

            // Computing intersection OBB
            if (first_obb.HasIntersection(second_obb)){
                return true;
            }
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::CreateDebugOBB2D(
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
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::CreateDebugOBB3D(
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
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
Parameters FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "first_model_part_name"  : "",
        "second_model_part_name" : "",
        "bounding_box_factor"    : -1.0,
        "debug_obb"              : false
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class FindIntersectedGeometricalObjectsWithOBBProcess<Condition>;
template class FindIntersectedGeometricalObjectsWithOBBProcess<Element>;

}  // namespace Kratos.
