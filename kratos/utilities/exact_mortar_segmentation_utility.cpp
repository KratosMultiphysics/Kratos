//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "containers/model.h"
// DEBUG
#include "includes/gid_io.h"

namespace Kratos {
template <>
bool ExactMortarIntegrationUtility<2, 2, false>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // We take the geometry GP from the core
    const double tolerance = 1.0e3 * ZeroTolerance;

    double total_weight = 0.0;
    array_1d<double, 2> auxiliar_coordinates(2 , 0.0);

    // Declaring auxiliar values
    PointType projected_gp_global;
    GeometryNodeType::CoordinatesArrayType projected_gp_local;

    // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
    for (IndexType i_slave = 0; i_slave < 2; ++i_slave) {
        const array_1d<double, 3>& normal = rOriginalSlaveGeometry[i_slave].FastGetSolutionStepValue(NORMAL);

        const double distance = GeometricalProjectionUtilities::FastProjectDirection(rOriginalMasterGeometry, rOriginalSlaveGeometry[i_slave].Coordinates(), projected_gp_global, rMasterNormal, -normal ); // The opposite direction

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }

        const bool is_inside = rOriginalMasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local, tolerance );

        if (is_inside) {// The slave node belongs to the master
            if (i_slave == 0) {
                auxiliar_coordinates[0] = -1.0;  // First node
            } else {
                auxiliar_coordinates[1] = 1.0;  // Second node
            }
        }
    }

    // We check if the element is fully integrated
    if ((auxiliar_coordinates[0] == -1.0 && auxiliar_coordinates[1] == 1.0)) {
        total_weight = 2.0;
    } else { // If not then we proceed
        std::vector<double> auxiliar_xi;
        for (IndexType i_master = 0; i_master < 2; ++i_master) {
            projected_gp_local[0] = (i_master == 0) ? -1.0 : 1.0;
            double delta_xi = (i_master == 0) ? 0.5 : -0.5;
            const bool is_inside = GeometricalProjectionUtilities::ProjectIterativeLine2D(rOriginalSlaveGeometry, rOriginalMasterGeometry[i_master].Coordinates(), projected_gp_local, rSlaveNormal, tolerance, delta_xi);

            if (is_inside)
                auxiliar_xi.push_back(projected_gp_local[0]);
        }

        // In this case one edge of the slave belongs to the master and additionally one node of the master belongs to the slave
        if (auxiliar_xi.size() == 1 && ((auxiliar_coordinates[0] == - 1.0 || auxiliar_coordinates[1] == 1.0))) {
            if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) { // NOTE: Equivalent to == -1.0
                auxiliar_coordinates[1] = auxiliar_xi[0];
            } else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) { // NOTE: Equivalent to == 1.0
                auxiliar_coordinates[0] = auxiliar_xi[0];
            } else {
                KRATOS_ERROR << "THIS IS NOT SUPPOSED TO HAPPEN. Auxiliar local coordinate: " <<  auxiliar_xi[0] << ". Auxiliar coordinates: " << auxiliar_coordinates[0] << " , " << auxiliar_coordinates[1] << std::endl;
            }
        } else if ( auxiliar_xi.size() == 2) { // Both nodes of the master belong to the slave (and none of the nodes of the slave belong to the master, the nodes can coincide, there is no other possibility)
            if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) { // NOTE: Equivalent to == -1.0. In this case the node in the left edge is already assigned
                auxiliar_coordinates[1] = auxiliar_xi[0] < auxiliar_xi[1] ? auxiliar_xi[1] : auxiliar_xi[0]; // We set in the proper position
            } else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) { // NOTE: Equivalent to == 1.0. In this case the node in the right edge is already assigned
                auxiliar_coordinates[0] = auxiliar_xi[0] < auxiliar_xi[1] ? auxiliar_xi[0] : auxiliar_xi[1]; // We set in the proper position
            } else { // There isn't any coincidence with the edges
                if (auxiliar_xi[0] < auxiliar_xi[1]) { // We check that are in proper order
                    auxiliar_coordinates[0] = auxiliar_xi[0];
                    auxiliar_coordinates[1] = auxiliar_xi[1];
                } else {
                    auxiliar_coordinates[1] = auxiliar_xi[0];
                    auxiliar_coordinates[0] = auxiliar_xi[1];
                }
            }
        } else { // THIS IS NOT SUPPOSED TO HAPPEN
            KRATOS_DEBUG_ERROR << "THIS IS NOT SUPPOSED TO HAPPEN!!!!\n" << rOriginalSlaveGeometry << "\n" << rOriginalMasterGeometry << std::endl;
            return false;  // NOTE: Giving problems
        }

        total_weight = auxiliar_coordinates[1] - auxiliar_coordinates[0];
    }

    KRATOS_ERROR_IF(total_weight < 0.0) << "Wrong order of the coordinates: "<< auxiliar_coordinates << std::endl;
    KRATOS_ERROR_IF(total_weight > 2.0) << "Impossible, Weight higher than 2: "<< auxiliar_coordinates << std::endl;

    // We do the final assignmen
    if (total_weight > ZeroTolerance) {
        rConditionsPointsSlave.resize(1);
        array_1d<PointType, 2> list_points;
        list_points[0].Coordinates()[0] = auxiliar_coordinates[0];
        list_points[1].Coordinates()[0] = auxiliar_coordinates[1];
        rConditionsPointsSlave[0] = list_points;

        return true;
    } else {
        rConditionsPointsSlave.clear();
        return false;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 3, false>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[1].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(3);
    PointerVector<PointType> points_array_master(3);
    for (IndexType i_node = 0; i_node < 3; ++i_node) {
        PointType aux_point;
        double distance;

        aux_point.Coordinates() = rOriginalSlaveGeometry[i_node].Coordinates();  // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);

        aux_point = GeometricalProjectionUtilities::FastProject( slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node)= Kratos::make_shared<PointType>(aux_point);

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    Triangle3D3<PointType> slave_geometry(points_array_slave);
    Triangle3D3<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 3> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside)) {
        rConditionsPointsSlave.resize(1);

        for (IndexType i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[0][i_node] = point;
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry);

        return TriangleIntersections<GeometryNodeType>(rConditionsPointsSlave, point_list, rOriginalSlaveGeometry, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 4, false>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[2].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(4);
    PointerVector<PointType> points_array_slave_not_rotated(4);
    PointerVector<PointType> points_array_master(4);
    for (IndexType i_node = 0; i_node < 4; ++i_node) {
        PointType aux_point;
        double distance_slave, distance_master;

        aux_point = GeometricalProjectionUtilities::FastProject(slave_center, rOriginalSlaveGeometry[i_node], rSlaveNormal, distance_slave);
        points_array_slave_not_rotated(i_node) = Kratos::make_shared<PointType>(aux_point);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);

        aux_point = GeometricalProjectionUtilities::FastProject( slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance_master);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node) = Kratos::make_shared<PointType>(aux_point);

        if (distance_master > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    Quadrilateral3D4<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> slave_geometry_not_rotated(points_array_slave_not_rotated);
    Quadrilateral3D4<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 4> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry, we divide into two triangles
    if (CheckAllInside(all_inside)) {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry);

        return TriangleIntersections<GeometryPointType>(rConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center, true);
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry);

        return TriangleIntersections<GeometryPointType>(rConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 3, false, 4>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[1].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // Auxiliar values
    PointType aux_point;
    double distance;

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(3);
    PointerVector<PointType> points_array_master(4);

    for (IndexType i_node = 0; i_node < 3; ++i_node) {
        aux_point.Coordinates() = rOriginalSlaveGeometry[i_node].Coordinates();  // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);
    }

    for (IndexType i_node = 0; i_node < 4; ++i_node) {
        aux_point = GeometricalProjectionUtilities::FastProject( slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node)= Kratos::make_shared<PointType>(aux_point);

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    Triangle3D3<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 4> all_inside_master;

    // We check if the nodes are inside
    CheckInside<4>(all_inside_master, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside_master)) { // We decompose in two triangles
        rConditionsPointsSlave.resize(2);

        for (IndexType i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[0][i_node] = point;
        }
        for (IndexType i_node = 1; i_node < 4; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[1][i_node - 1] = point;
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints<4>(point_list, all_inside_master, master_geometry);

        // We check if the nodes are inside
        array_1d<bool, 3> all_inside_slave;
        CheckInside(all_inside_slave, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside_slave, slave_geometry);

        return TriangleIntersections<GeometryNodeType>(rConditionsPointsSlave, point_list, rOriginalSlaveGeometry, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 4, false, 3>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[2].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // Auxiliar values
    PointType aux_point;
    double distance;

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(4);
    PointerVector<PointType> points_array_slave_not_rotated(4);
    PointerVector<PointType> points_array_master(3);

    for (IndexType i_node = 0; i_node < 4; ++i_node) {
        aux_point = GeometricalProjectionUtilities::FastProject(slave_center, rOriginalSlaveGeometry[i_node], rSlaveNormal, distance);
        points_array_slave_not_rotated(i_node) = Kratos::make_shared<PointType>(aux_point);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    for (IndexType i_node = 0; i_node < 3; ++i_node) {
        aux_point = GeometricalProjectionUtilities::FastProject( slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node) = Kratos::make_shared<PointType>(aux_point);

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    Quadrilateral3D4<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> slave_geometry_not_rotated(points_array_slave_not_rotated);
    Triangle3D3<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 3> all_inside_master;

    // We check if the nodes are inside
    CheckInside<3>(all_inside_master, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside_master)) { // We generate only one triangle
        rConditionsPointsSlave.resize(1);

        for (IndexType i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[0][i_node] = point;
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints<3>(point_list, all_inside_master, master_geometry);

        // We check if the nodes are inside
        array_1d<bool, 4> all_inside_slave;
        CheckInside(all_inside_slave, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside_slave, slave_geometry);

        return TriangleIntersections<GeometryPointType>(rConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

// NOTE: The following are "hardcopies" of the previous ones, this is because C++ doesn't allow yet the partial specialization (CHECK THE ERROR TWICE!!!!!)

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<2, 2, true>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // We take the geometry GP from the core
    const double tolerance = 1.0e3 * ZeroTolerance;

    double total_weight = 0.0;
    array_1d<double, 2> auxiliar_coordinates(2, 0.0);
    array_1d<PointBelongsLine2D2N, 2> auxiliar_belong;

    // Declaring auxiliar values
    PointType projected_gp_global;
    GeometryNodeType::CoordinatesArrayType projected_gp_local;

    // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
    for (unsigned int i_slave = 0; i_slave < 2; ++i_slave) {
        const array_1d<double, 3>& normal = rOriginalSlaveGeometry[i_slave].FastGetSolutionStepValue(NORMAL);

        const double distance = GeometricalProjectionUtilities::FastProjectDirection(rOriginalMasterGeometry, rOriginalSlaveGeometry[i_slave].Coordinates(), projected_gp_global, rMasterNormal, -normal ); // The opposite direction

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }

        const bool is_inside = rOriginalMasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local, tolerance );

        if (is_inside) {// The slave node belongs to the master
            if (i_slave == 0) {  // First node
                auxiliar_coordinates[0] = -1.0;
                auxiliar_belong[0] = PointBelongsLine2D2N::SlaveLine2D2N0;
            } else { // Second node
                auxiliar_coordinates[1] = 1.0;
                auxiliar_belong[1] = PointBelongsLine2D2N::SlaveLine2D2N1;
            }
        }
    }

    // We check if the element is fully integrated
    if ((auxiliar_coordinates[0] == - 1.0 && auxiliar_coordinates[1] == 1.0)) {
        total_weight = 2.0;
    } else { // If not then we proceed
        std::vector<double> auxiliar_xi;
        std::vector<PointBelongsLine2D2N> auxiliar_master_belong;
        for (IndexType i_master = 0; i_master < 2; ++i_master) {
            projected_gp_local[0] = (i_master == 0) ? -1.0 : 1.0;
            double delta_xi = (i_master == 0) ? 0.5 : -0.5;
            const bool is_inside = GeometricalProjectionUtilities::ProjectIterativeLine2D(rOriginalSlaveGeometry,  rOriginalMasterGeometry[i_master].Coordinates(), projected_gp_local, rSlaveNormal, tolerance, delta_xi);

            if (is_inside) {
                auxiliar_xi.push_back(projected_gp_local[0]);
                auxiliar_master_belong.push_back( static_cast<PointBelongsLine2D2N>(2 + i_master));
            }
        }

        // In this case one edge of the slave belongs to the master and additionally one node of the master belongs to the slave
        if (auxiliar_xi.size() == 1 && ((auxiliar_coordinates[0] == -1.0 ||  auxiliar_coordinates[1] == 1.0))) {
            if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) { // NOTE: Equivalent to == -1.0
                auxiliar_coordinates[1] = auxiliar_xi[0];
                auxiliar_belong[1] = auxiliar_master_belong[0];
            } else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) { // NOTE: Equivalent to == 1.0
                auxiliar_coordinates[0] = auxiliar_xi[0];
                auxiliar_belong[0] = auxiliar_master_belong[0];
            } else {
                KRATOS_ERROR << "THIS IS NOT SUPPOSED TO HAPPEN. Auxiliar local coordinate: " <<  auxiliar_xi[0] << ". Auxiliar coordinates: " << auxiliar_coordinates[0] << " , " << auxiliar_coordinates[1] << std::endl;
            }
        } else if ( auxiliar_xi.size() == 2) { // Both nodes of the master belong to the slave (and none of the nodes of the slave belong to the master, the nodes can coincide, there is no other possibility)
            if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) { // NOTE: Equivalent to == -1.0. In this case the node in the left edge is already assigned
                auxiliar_coordinates[1] = auxiliar_xi[0] < auxiliar_xi[1]  ? auxiliar_xi[1] : auxiliar_xi[0];  // We set in the proper position
                auxiliar_belong[1] = auxiliar_xi[0] < auxiliar_xi[1]  ? auxiliar_master_belong[1] : auxiliar_master_belong[0];
            } else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) { // NOTE: Equivalent to == 1.0. In this case the node in the right edge is already assigned
                auxiliar_coordinates[0] = auxiliar_xi[0] < auxiliar_xi[1]  ? auxiliar_xi[0]  : auxiliar_xi[1];  // We set in the proper position
                auxiliar_belong[0] = auxiliar_xi[0] < auxiliar_xi[1]  ? auxiliar_master_belong[0] : auxiliar_master_belong[1];
            } else { // There isn't any coincidence with the edges
                if (auxiliar_xi[0] < auxiliar_xi[1]) { // We check that are in proper order
                    auxiliar_coordinates[0] = auxiliar_xi[0];
                    auxiliar_coordinates[1] = auxiliar_xi[1];
                    auxiliar_belong[0] = auxiliar_master_belong[0];
                    auxiliar_belong[1] = auxiliar_master_belong[1];
                } else {
                    auxiliar_coordinates[1] = auxiliar_xi[0];
                    auxiliar_coordinates[0] = auxiliar_xi[1];
                    auxiliar_belong[1] = auxiliar_master_belong[0];
                    auxiliar_belong[0] = auxiliar_master_belong[1];
                }
            }
        } else { // THIS IS NOT SUPPOSED TO HAPPEN
            KRATOS_DEBUG_ERROR << "THIS IS NOT SUPPOSED TO HAPPEN!!!!\n" << rOriginalSlaveGeometry << "\n" << rOriginalMasterGeometry << std::endl;
            return false;  // NOTE: Giving problems
        }

        total_weight = auxiliar_coordinates[1] - auxiliar_coordinates[0];
    }

    KRATOS_ERROR_IF(total_weight < 0.0) << "Wrong order of the coordinates: "<< auxiliar_coordinates << std::endl;
    KRATOS_ERROR_IF(total_weight > 2.0) << "Impossible, Weight higher than 2: "<< auxiliar_coordinates << std::endl;

    // We do the final assignmen
    if (total_weight > ZeroTolerance) {
        rConditionsPointsSlave.resize(1);
        array_1d<PointBelong<2>, 2> list_points;
        list_points[0].Coordinates()[0] = auxiliar_coordinates[0];
        list_points[0].SetBelong(auxiliar_belong[0]);
        list_points[1].Coordinates()[0] = auxiliar_coordinates[1];
        list_points[1].SetBelong(auxiliar_belong[1]);
        rConditionsPointsSlave[0] = list_points;

        return true;
    } else {
        rConditionsPointsSlave.clear();
        return false;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 3, true>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[1].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(3);
    PointerVector<PointType> points_array_master(3);
    for (IndexType i_node = 0; i_node < 3; ++i_node) {
        PointType aux_point;
        double distance;

        aux_point.Coordinates() = rOriginalSlaveGeometry[i_node].Coordinates();  // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);

        aux_point = GeometricalProjectionUtilities::FastProject(slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node) = Kratos::make_shared<PointType>(aux_point);
    }

    Triangle3D3<PointType> slave_geometry(points_array_slave);
    Triangle3D3<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 3> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside)) {
        rConditionsPointsSlave.resize(1);

        for (IndexType i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[0][i_node] = PointBelong<3>(point.Coordinates(), static_cast<PointBelongsTriangle3D3N>(i_node + 3));
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry, PointBelongs::Master);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry, PointBelongs::Slave);

        return TriangleIntersections<GeometryNodeType>(rConditionsPointsSlave, point_list, rOriginalSlaveGeometry, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 4, true>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[2].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(4);
    PointerVector<PointType> points_array_slave_not_rotated(4);
    PointerVector<PointType> points_array_master(4);
    for (IndexType i_node = 0; i_node < 4; ++i_node) {
        PointType aux_point;
        double distance_slave, distance_master;

        aux_point = GeometricalProjectionUtilities::FastProject(slave_center, rOriginalSlaveGeometry[i_node], rSlaveNormal, distance_slave);
        points_array_slave_not_rotated(i_node) = Kratos::make_shared<PointType>(aux_point);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);

        aux_point = GeometricalProjectionUtilities::FastProject(slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance_master);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node) = Kratos::make_shared<PointType>(aux_point);

        if (distance_master > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    Quadrilateral3D4<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> slave_geometry_not_rotated(points_array_slave_not_rotated);
    Quadrilateral3D4<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 4> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside)) {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry, PointBelongs::Master);

        return TriangleIntersections<GeometryPointType>(rConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry,  master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center,  true);
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry, PointBelongs::Master);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry, PointBelongs::Slave);

        return TriangleIntersections<GeometryPointType>(rConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 3, true, 4>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[1].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // Auxiliar values
    PointType aux_point;
    double distance;

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(3);
    PointerVector<PointType> points_array_master(4);

    for (IndexType i_node = 0; i_node < 3; ++i_node) {
        aux_point.Coordinates() = rOriginalSlaveGeometry[i_node].Coordinates();  // NOTE: We are in a linear triangle, ali_nodel the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);
    }

    for (IndexType i_node = 0; i_node < 4; ++i_node) {
        aux_point = GeometricalProjectionUtilities::FastProject( slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node)= Kratos::make_shared<PointType>(aux_point);

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    Triangle3D3<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 4> all_inside_master;

    // We check if the nodes are inside
    CheckInside<4>(all_inside_master, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside_master)) { // We decompose in two triangles
        rConditionsPointsSlave.resize(2);

        for (IndexType i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[0][i_node] = PointBelong<3, 4>(point.Coordinates(), static_cast<PointBelongsTriangle3D3NQuadrilateral3D4N>(i_node + 3));
        }
        for (IndexType i_node = 1; i_node < 4; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[1][i_node - 1] = PointBelong<3, 4>(point.Coordinates(), static_cast<PointBelongsTriangle3D3NQuadrilateral3D4N>(i_node + 3));
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints<4>(point_list, all_inside_master, master_geometry, PointBelongs::Master);

        // We check if the nodes are inside
        array_1d<bool, 3> all_inside_slave;
        CheckInside(all_inside_slave, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside_slave, slave_geometry, PointBelongs::Slave);

        return TriangleIntersections<GeometryNodeType>(rConditionsPointsSlave, point_list, rOriginalSlaveGeometry, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 4, true, 3>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    ConditionArrayListType& rConditionsPointsSlave
    )
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType slave_center = rOriginalSlaveGeometry.Center();

    // We define the condition tangents
    array_1d<double, 3> slave_tangent_xi = rOriginalSlaveGeometry[2].Coordinates() - rOriginalSlaveGeometry[0].Coordinates();
    slave_tangent_xi = slave_tangent_xi/norm_2(slave_tangent_xi);
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, rSlaveNormal, slave_tangent_xi);

    // Auxiliar values
    PointType aux_point;
    double distance;

    // We define the auxiliar geometry
    PointerVector<PointType> points_array_slave(4);
    PointerVector<PointType> points_array_slave_not_rotated(4);
    PointerVector<PointType> points_array_master(3);

    for (IndexType i_node = 0; i_node < 4; ++i_node) {
        aux_point = GeometricalProjectionUtilities::FastProject(slave_center, rOriginalSlaveGeometry[i_node], rSlaveNormal, distance);
        points_array_slave_not_rotated(i_node) = Kratos::make_shared<PointType>(aux_point);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave(i_node) = Kratos::make_shared<PointType>(aux_point);

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }

    for (IndexType i_node = 0; i_node < 3; ++i_node) {
        aux_point = GeometricalProjectionUtilities::FastProject( slave_center, rOriginalMasterGeometry[i_node], rSlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master(i_node) = Kratos::make_shared<PointType>(aux_point);

        if (distance > mDistanceThreshold) {
            rConditionsPointsSlave.clear();
            return false;
        }
    }


    Quadrilateral3D4<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> slave_geometry_not_rotated(points_array_slave_not_rotated);
    Triangle3D3<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 3> all_inside_master;

    // We check if the nodes are inside
    CheckInside<3>(all_inside_master, slave_geometry, master_geometry, ZeroTolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside_master)) { // We generate only one triangle
        rConditionsPointsSlave.resize(1);

        for (IndexType i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            rOriginalSlaveGeometry.PointLocalCoordinates( point, rOriginalMasterGeometry[i_node]);
            rConditionsPointsSlave[0][i_node] = PointBelong<4, 3>(point.Coordinates(), static_cast<PointBelongsQuadrilateral3D4NTriangle3D3N>(i_node + 4));;
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints<3>(point_list, all_inside_master, master_geometry, PointBelongs::Master);

        // We check if the nodes are inside
        array_1d<bool, 4> all_inside_slave;
        CheckInside(all_inside_slave, master_geometry, slave_geometry, ZeroTolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside_slave, slave_geometry, PointBelongs::Slave);

        return TriangleIntersections<GeometryPointType>(rConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetExactIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    IntegrationPointsType& rIntegrationPointsSlave
    )
{
    ConditionArrayListType conditions_points_slave;

    const bool is_inside = GetExactIntegration(rOriginalSlaveGeometry, rSlaveNormal, rOriginalMasterGeometry, rMasterNormal, conditions_points_slave);

    for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
        PointerVector<PointType> points_array(TDim);  // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            rOriginalSlaveGeometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
            points_array(i_node) = Kratos::make_shared<PointType>(global_point);
        }

        DecompositionType decomp_geom(points_array);

        const GeometryPointType::IntegrationPointsArrayType& local_integration_slave = decomp_geom.IntegrationPoints(mAuxIntegrationMethod);

        // Integrating the mortar operators
        for (unsigned int point_number = 0; point_number < local_integration_slave.size(); ++point_number) {
            const double weight =  local_integration_slave[point_number].Weight();
            const PointType local_point_decomp = local_integration_slave[point_number].Coordinates();
            PointType local_point_parent;
            PointType gp_global;
            decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
            rOriginalSlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);

            const double det_J = decomp_geom.DeterminantOfJacobian( local_point_decomp ) * (TDim == 2 ? 2.0 : 1.0);

            const array_1d<double, 3>& coordinates = local_point_parent.Coordinates();
            rIntegrationPointsSlave.push_back( IntegrationPointType( coordinates[0], coordinates[1], weight * det_J )); // TODO: Change push_back for a fix operation
        }
    }

    return is_inside;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetExactAreaIntegration(
    GeometryNodeType& rOriginalSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    GeometryNodeType& rOriginalMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    double& rArea
    )
{
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = GetExactIntegration(rOriginalSlaveGeometry, rSlaveNormal, rOriginalMasterGeometry, rMasterNormal, conditions_points_slave);
    if (is_inside) {
        GetTotalArea(rOriginalSlaveGeometry, conditions_points_slave, rArea);
        // Debugging
        if (mEchoLevel > 1) {
            MathematicaDebug(0, rOriginalSlaveGeometry, 0, rOriginalMasterGeometry, conditions_points_slave);
        }
    }

    return is_inside;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetTotalArea(
    GeometryNodeType& rOriginalSlaveGeometry,
    ConditionArrayListType& rConditionsPointsSlave,
    double& rArea
    )
{
    rArea = 0.0;

    for (IndexType i_geom = 0; i_geom < rConditionsPointsSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array( TDim);  // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            rOriginalSlaveGeometry.GlobalCoordinates(global_point, rConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>(global_point);
        }

        DecompositionType decomp_geom( PointerVector<PointType>{points_array} );

        rArea += decomp_geom.Area();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TestGetExactIntegration(
    Condition::Pointer pSlaveCond,
    Condition::Pointer pMasterCond,
    Matrix& rCustomSolution
    )
{
    IntegrationPointsType integration_points_slave;

    const bool solution = GetExactIntegration(pSlaveCond->GetGeometry(), pSlaveCond->GetValue(NORMAL), pMasterCond->GetGeometry(), pMasterCond->GetValue(NORMAL), integration_points_slave);

    rCustomSolution.resize(integration_points_slave.size(), TDim, false);

    for (IndexType GP = 0; GP < integration_points_slave.size(); ++GP) {
        // Solution save:
        const array_1d<double, 3>& coordinates_gp = integration_points_slave[GP].Coordinates();
        rCustomSolution(GP, 0) = coordinates_gp[0];
        if (TDim == 2)
            rCustomSolution(GP, 1) = integration_points_slave[GP].Weight();
        else {
            rCustomSolution(GP, 1) = coordinates_gp[1];
            rCustomSolution(GP, 2) = integration_points_slave[GP].Weight();
        }
    }

    return solution;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
double ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TestGetExactAreaIntegration(
    ModelPart& rMainModelPart,
    Condition::Pointer pSlaveCond
    )
{
    // Initalize values
    double area = 0.0;

    if ( pSlaveCond->Has( INDEX_MAP )) {
        IndexMap::Pointer indexes_map = pSlaveCond->GetValue( INDEX_MAP );

        for (auto it_pair = indexes_map->begin(); it_pair != indexes_map->end(); ++it_pair ) {
            double local_area = 0.0;
            Condition::Pointer p_master_cond = rMainModelPart.pGetCondition(it_pair->first);
            const bool is_inside = GetExactAreaIntegration(pSlaveCond->GetGeometry(), pSlaveCond->GetValue(NORMAL), p_master_cond->GetGeometry(), p_master_cond->GetValue(NORMAL), local_area);
            if (is_inside) area += local_area;
        }
    } else {
        IndexSet::Pointer indexes_set = pSlaveCond->GetValue( INDEX_SET );

        for (auto it_pair = indexes_set->begin(); it_pair != indexes_set->end(); ++it_pair ) {
            double local_area = 0.0;
            Condition::Pointer p_master_cond = rMainModelPart.pGetCondition(*it_pair);
            const bool is_inside = GetExactAreaIntegration(pSlaveCond->GetGeometry(), pSlaveCond->GetValue(NORMAL), p_master_cond->GetGeometry(), p_master_cond->GetValue(NORMAL), local_area);
            if (is_inside) area += local_area;
        }
    }

    return area;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TestGiDDebug(ModelPart& rMainModelPart)
{
    if (TDim == 3) {
        Model& current_model = rMainModelPart.GetOwnerModel();
        ModelPart& aux_model_part = current_model.CreateModelPart("exact_mortar_aux_model_part");

        IndexType node_counter = 1;
        IndexType cond_counter = 1;

        for (auto& cond : rMainModelPart.Conditions()) {
            if (cond.Is(SLAVE)) {
                auto& slave_geometry = cond.GetGeometry();

                if ( cond.Has( INDEX_MAP )) {
                    IndexMap::Pointer indexes_map = cond.GetValue( INDEX_MAP );

                    for (auto it_pair = indexes_map->begin(); it_pair != indexes_map->end(); ++it_pair ) {
                        Condition::Pointer p_master_cond = rMainModelPart.pGetCondition(it_pair->first);
                        ConditionArrayListType conditions_points_slave;
                        const bool is_inside = GetExactIntegration(slave_geometry, cond.GetValue(NORMAL), p_master_cond->GetGeometry(), p_master_cond->GetValue(NORMAL), conditions_points_slave);
                        if (is_inside) {
                            for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
                                std::vector<NodeType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                                for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                                    PointType global_point;
                                    slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                    points_array[i_node] = aux_model_part.CreateNewNode(node_counter, global_point.X(), global_point.Y(), global_point.Z());
                                    node_counter++;
                                }
                                aux_model_part.CreateNewCondition("Condition3D", cond_counter, Geometry< Node < 3 > >::PointsArrayType{points_array}, cond.pGetProperties());
                                cond_counter++;
                            }
                        }
                    }
                } else {
                    IndexSet::Pointer indexes_set = cond.GetValue( INDEX_SET );

                    for (auto it_pair = indexes_set->begin(); it_pair != indexes_set->end(); ++it_pair ) {
                        Condition::Pointer p_master_cond = rMainModelPart.pGetCondition(*it_pair);
                        ConditionArrayListType conditions_points_slave;
                        const bool is_inside = GetExactIntegration(slave_geometry, cond.GetValue(NORMAL), p_master_cond->GetGeometry(), p_master_cond->GetValue(NORMAL), conditions_points_slave);
                        if (is_inside) {
                            for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
                                std::vector<NodeType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                                for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                                    PointType global_point;
                                    slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                    points_array[i_node] = aux_model_part.CreateNewNode(node_counter, global_point.X(), global_point.Y(), global_point.Z());
                                    node_counter++;
                                }
                                aux_model_part.CreateNewCondition("Condition3D", cond_counter, Geometry< Node < 3 > >::PointsArrayType{points_array}, cond.pGetProperties());
                                cond_counter++;
                            }
                        }
                    }
                }
            }

            current_model.DeleteModelPart("exact_mortar_aux_model_part");
        }

        auto pgidio = Kratos::make_shared<GidIO<>>("ExactMortarIntegrationUtilityDEBUG", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
        pgidio->InitializeMesh(0);
        pgidio->WriteMesh(aux_model_part.GetMesh());
        pgidio->FinalizeMesh();
        pgidio->CloseResultFile();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetIntegrationMethod()
{
    // Setting the auxiliar integration points
    switch (mIntegrationOrder) {
        case 1:
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_1;
            break;
        case 2:
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_2;
            break;
        case 3:
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_3;
            break;
        case 4:
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_4;
            break;
        case 5:
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_5;
            break;
        default:
            mAuxIntegrationMethod = GeometryData::GI_GAUSS_2;
            break;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
GeometryNodeType::IntegrationPointsArrayType ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetIntegrationTriangle()
{
    // Setting the auxiliar integration points
    switch (mIntegrationOrder) {
        case 1:
            return Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 2:
            return Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 3:
            return Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 4:
            return Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 5:
            return Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        default:
            return Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
inline std::vector<IndexType> ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::ComputeAnglesIndexes(
    PointListType& rPointList,
    const array_1d<double, 3>& rNormal
    ) const
{
    const SizeType list_size = rPointList.size();

    // We reorder the nodes according with the angle they form with the first node
    std::vector<double> angles(list_size - 1);
    const array_1d<double, 3>& ref_point_coordinates = rPointList[0].Coordinates();
    array_1d<double, 3> v = rPointList[1].Coordinates() - ref_point_coordinates;

    v /= norm_2(v);
    array_1d<double, 3> n;
    MathUtils<double>::CrossProduct( n, rNormal, v);

    for (IndexType elem = 1; elem < list_size; ++elem) {
        angles[elem - 1] = AnglePoints(rPointList[0], rPointList[elem], v, n);
        if (angles[elem - 1] < 0.0) {
            v = rPointList[elem].Coordinates() - ref_point_coordinates;
            v /= norm_2(v);
            MathUtils<double>::CrossProduct( n, rNormal, v);
            for (IndexType aux_elem = 0; aux_elem <= (elem - 1); ++aux_elem)
                angles[aux_elem] -= angles[elem - 1];
        }
    }

    return MortarUtilities::SortIndexes<double>(angles);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
inline void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::ComputeClippingIntersections(
    PointListType& rPointList,
    GeometryPointType& rSlaveGeometry,
    GeometryPointType& rMasterGeometry,
    const PointType& rRefCenter
    )
{
    // We consider the Z coordinate constant
    const double z_ref = rRefCenter.Coordinates()[2];

    // We find the intersection in each side
    for (IndexType i_edge = 0; i_edge < TNumNodes; ++i_edge) {
        const IndexType ip_edge = (i_edge == (TNumNodes - 1)) ? 0 : i_edge + 1;
        for (IndexType j_edge = 0; j_edge < TNumNodesMaster; ++j_edge) {
            const IndexType jp_edge = (j_edge == (TNumNodesMaster - 1)) ? 0 : j_edge + 1;

            PointType intersected_point;
            const bool intersected = Clipping2D( intersected_point, rSlaveGeometry[i_edge], rSlaveGeometry[ip_edge], rMasterGeometry[j_edge], rMasterGeometry[jp_edge] );

            if (intersected) {
                // Set the coordinate
                intersected_point.Coordinates()[2] = z_ref;

                // Ititialize the check
                bool add_point = true;
                for (IndexType iter = 0; iter < rPointList.size(); ++iter) {
                    if (CheckPoints(intersected_point, rPointList[iter])) {
                        add_point = false;
                        break;
                    }
                }

                if (add_point) {
                    if (TBelong) { // NOTE: We do some kind of strange hash to know the intersected edges
                        const std::size_t hash = (TNumNodesMaster + TNumNodes) + 10 * i_edge + 100 * ip_edge + 1000 * j_edge + 10000 * jp_edge;
                        rPointList.push_back(PointBelong<TNumNodes, TNumNodesMaster>(intersected_point.Coordinates(), static_cast<BelongType>(hash)));
                    } else
                        rPointList.push_back(intersected_point);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
template <class TGeometryType>
inline bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TriangleIntersections(
    ConditionArrayListType& rConditionsPointsSlave,
    PointListType& rPointList,
    TGeometryType& rOriginalSlaveGeometry,
    GeometryPointType& rSlaveGeometry,
    GeometryPointType& rMasterGeometry,
    const array_1d<double, 3>& rSlaveTangentXi,
    const array_1d<double, 3>& rSlaveTangentEta,
    const PointType& rRefCenter,
    const bool IsAllInside
    )
{
    // We do the clipping
    if (IsAllInside == false)
        ComputeClippingIntersections(rPointList, rSlaveGeometry, rMasterGeometry, rRefCenter);

    // We compose the triangles
    const SizeType list_size = rPointList.size();
    if (list_size >  2) { // Technically the minimum is three, just in case I consider 2
        // Declaring auxiliar local point
        PointType local_point;
        PointListType aux_master_point_list(list_size);

        // We recover this point to the triangle plane and compute the local coordinates
        for (IndexType i_point_list = 0; i_point_list < rPointList.size(); ++i_point_list) {
            rMasterGeometry.PointLocalCoordinates(local_point, rPointList[i_point_list].Coordinates());
            aux_master_point_list[i_point_list].Coordinates() = local_point.Coordinates();
            MortarUtilities::RotatePoint(rPointList[i_point_list], rRefCenter, rSlaveTangentXi, rSlaveTangentEta, true);
            rOriginalSlaveGeometry.PointLocalCoordinates(local_point, rPointList[i_point_list].Coordinates());
            rPointList[i_point_list].Coordinates() = local_point.Coordinates();
        }

        // We will check if the triangle is inside the slave geometry, so we will compute an auxiliar shape function
        array_1d<double, 2> auxiliar_slave_center_local_coords, auxiliar_master_center_local_coords;

        // We compute the angles between the nodes
        rSlaveGeometry.PointLocalCoordinates(local_point, rSlaveGeometry.Center());
        const array_1d<double, 3>& normal = rSlaveGeometry.UnitNormal(local_point);
        const std::vector<IndexType> index_vector = ComputeAnglesIndexes(rPointList, normal);

        // We resize the array of points of the decomposed triangles
        rConditionsPointsSlave.resize((list_size - 2));

        IndexType aux_elem_index = 0;
        for (IndexType elem = 0; elem < list_size - 2; ++elem) { // NOTE: We always have two points less that the number of nodes
            ArrayTriangleType points_locals_slave, points_locals_master;

            points_locals_slave[0] = rPointList[0];
            points_locals_slave[1] = rPointList[index_vector[elem + 0] + 1];
            points_locals_slave[2] = rPointList[index_vector[elem + 1] + 1];
            points_locals_master[0] = aux_master_point_list[0];
            points_locals_master[1] = aux_master_point_list[index_vector[elem + 0] + 1];
            points_locals_master[2] = aux_master_point_list[index_vector[elem + 1] + 1];

            // We compute if the center is inside the slave geometry
            auxiliar_slave_center_local_coords[0] = 1.0/3.0 * (points_locals_slave[0].X() + points_locals_slave[1].X() + points_locals_slave[2].X());
            auxiliar_slave_center_local_coords[1] = 1.0/3.0 * (points_locals_slave[0].Y() + points_locals_slave[1].Y() + points_locals_slave[2].Y());
            auxiliar_master_center_local_coords[0] = 1.0/3.0 * (points_locals_master[0].X() + points_locals_master[1].X() + points_locals_master[2].X());
            auxiliar_master_center_local_coords[1] = 1.0/3.0 * (points_locals_master[0].Y() + points_locals_master[1].Y() + points_locals_master[2].Y());
            const bool center_is_inside = CheckCenterIsInside(auxiliar_slave_center_local_coords) && CheckCenterIsInside(auxiliar_master_center_local_coords, TNumNodesMaster);
            if (!center_is_inside) {
                rConditionsPointsSlave.erase(rConditionsPointsSlave.begin() + aux_elem_index);
                KRATOS_WARNING_IF("ExactMortarIntegrationUtility", mEchoLevel > 0) << "The generated intersection is probably a concave polygon. Check it out: \n" << rSlaveGeometry << "\n" << rMasterGeometry << std::endl;
                continue; // We skip this triangle
            }

            // We add the triangle to the vector
            rConditionsPointsSlave[aux_elem_index] = points_locals_slave;

            // We update the auxiliar index
            ++aux_elem_index;
        }

        if (rConditionsPointsSlave.size() > 0)
            return true;
        else
            return false;
    } else { // No intersection
        rConditionsPointsSlave.clear();
        return false;
    }

    rConditionsPointsSlave.clear();
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, bool TBelong, SizeType TNumNodesMaster>
inline bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::CheckCenterIsInside(
    const array_1d<double, 2>& rAuxiliarCenterLocalCoordinates,
    const SizeType NumNodes
    )
{
    if (NumNodes == 3) {
        if ( (rAuxiliarCenterLocalCoordinates[0] >= (0.0-ZeroTolerance)) && (rAuxiliarCenterLocalCoordinates[0] <= (1.0+ZeroTolerance)) ) {
            if ( (rAuxiliarCenterLocalCoordinates[1] >= (0.0-ZeroTolerance)) && (rAuxiliarCenterLocalCoordinates[1] <= (1.0+ZeroTolerance)) ) {
                if ( (rAuxiliarCenterLocalCoordinates[0] + rAuxiliarCenterLocalCoordinates[1]) <= (1.0+ZeroTolerance) ) {
                    return true;
                }
            }
        }
    } if (NumNodes == 4) {
        if ( std::abs(rAuxiliarCenterLocalCoordinates[0]) <= (1.0+ZeroTolerance) ) {
            if ( std::abs(rAuxiliarCenterLocalCoordinates[1]) <= (1.0+ZeroTolerance) ) {
                return true;
            }
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template class ExactMortarIntegrationUtility<2, 2, false>;
template class ExactMortarIntegrationUtility<3, 3, false>;
template class ExactMortarIntegrationUtility<3, 4, false>;
template class ExactMortarIntegrationUtility<3, 3, false, 4>;
template class ExactMortarIntegrationUtility<3, 4, false, 3>;

template class ExactMortarIntegrationUtility<2, 2, true>;
template class ExactMortarIntegrationUtility<3, 3, true>;
template class ExactMortarIntegrationUtility<3, 4, true>;
template class ExactMortarIntegrationUtility<3, 3, true, 4>;
template class ExactMortarIntegrationUtility<3, 4, true, 3>;

}  // namespace Kratos.
