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
#include "includes/obb.h"
#include "processes/find_intersected_geometrical_objects_with_obb_process.h"

namespace Kratos
{

template<class TEntity>
FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::FindIntersectedGeometricalObjectsWithOBBProcess(
    ModelPart& rPart1,
    ModelPart& rPart2,
    const double BoundingBoxFactor
    ) : BaseType(rPart1, rPart2),
        mBoundingBoxFactor(BoundingBoxFactor)
{
    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);
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

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_first_model_part_name = ThisParameters["first_model_part_name"].GetString();
    const std::string& r_second_model_part_name = ThisParameters["second_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_first_model_part_name == "") << "first_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_second_model_part_name == "") << "second_model_part_name must be defined on parameters" << std::endl;

    // Setting the bounding box factor
    mBoundingBoxFactor = ThisParameters["bounding_box_factor"].GetDouble();

    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);
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
//     // Check the intersection of each edge of the object bounding box against the intersecting object bounding box
//     NodeType first_geometry_low_node, first_geometry_high_node;
//     rFirstGeometry.BoundingBox(first_geometry_low_node, first_geometry_high_node);
//     NodeType second_geometry_low_node, second_geometry_high_node;
//     rSecondGeometry.BoundingBox(second_geometry_low_node, second_geometry_high_node);
//
//     // In order to compute the extended bounding box we need to compute the tangents
//     array_1d<double, 3> tangent_1 = first_geometry_high_node.Coordinates() - first_geometry_low_node.Coordinates();
//     tangent_1 /= norm_2(tangent_1);
//
//     array_1d<double, 3> tangent_2 = second_geometry_high_node.Coordinates() - second_geometry_low_node.Coordinates();
//     tangent_2 /= norm_2(tangent_2);
//
//     // Computing the normal to this tangents
//     array_1d<double, 3> normal_1, normal_2;
//     normal_1[0] = first_geometry_high_node[1] - first_geometry_low_node[1];
//     normal_1[1] = first_geometry_low_node[0] - first_geometry_high_node[0];
//     normal_1[2] = 0.0;
//     normal_2[0] = second_geometry_high_node[1] - second_geometry_low_node[1];
//     normal_2[1] = second_geometry_low_node[0] - second_geometry_high_node[0];
//     normal_2[2] = 0.0;
//
//     // Creating the box
//     std::vector<array_1d<array_1d<double, 3>, 2>> edges_1(4);
//     std::vector<array_1d<array_1d<double, 3>, 2>> edges_2(4);
//
//     noalias(edges_1[0][0]) = first_geometry_high_node + mBoundingBoxFactor * (normal_1 + tangent_1);
//     noalias(edges_1[0][1]) = first_geometry_low_node + mBoundingBoxFactor * (normal_1 - tangent_1);
//     noalias(edges_1[1][0]) = edges_1[0][1];
//
//     // Now we do the check
//     array_1d<double,3> int_pt = ZeroVector(3);
//     for (auto& edge_1 : edges_1) {
//         for (auto& edge_2 : edges_2) {
//             const int int_id = IntersectionUtilities::ComputeLineLineIntersection(edge_1, edge_2[0], edge_2[1], int_pt);
//
//             if (int_id != 0){
//                 return true;
//             }
//         }
//     }
//
//     // Let check second geometry is inside the first one.
//     // Considering that there are no intersection, if one point is inside all of it is inside.
//     array_1d<double, 3> local_point;
//     if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
//         return true;
//     }

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
//     // Check the intersection of each face against the intersecting object
//     auto faces = rFirstGeometry.Faces();
//     for (auto& face : faces) {
//         if (face.HasIntersection(rSecondGeometry)){
//             return true;
//         }
//     }
//
//     // Let check second geometry is inside the first one.
//     // Considering that there are no intersection, if one point is inside all of it is inside.
//     array_1d<double, 3> local_point;
//     if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
//         return true;
//     }

    return false;
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
        "bounding_box_factor"    : -1.0
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class FindIntersectedGeometricalObjectsWithOBBProcess<Condition>;
template class FindIntersectedGeometricalObjectsWithOBBProcess<Element>;

}  // namespace Kratos.
