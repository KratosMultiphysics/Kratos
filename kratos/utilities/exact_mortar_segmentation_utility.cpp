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
#include "utilities/exact_mortar_segmentation_utility.h"

namespace Kratos {
template <>
bool ExactMortarIntegrationUtility<2, 2, false>::GetExactIntegration(
    GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    ConditionArrayListType& ConditionsPointsSlave) {
    // We take the geometry GP from the core
    const double tolerance = 1.0e-6;  // std::numeric_limits<double>::epsilon();

    double total_weight = 0.0;
    array_1d<double, 2> auxiliar_coordinates = ZeroVector(2);

    // Declaring auxiliar values
    PointType projected_gp_global;
    GeometryNodeType::CoordinatesArrayType projected_gp_local;

    // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
    for (unsigned int i_slave = 0; i_slave < 2; ++i_slave) {
        const array_1d<double, 3>& normal = OriginalSlaveGeometry[i_slave].FastGetSolutionStepValue(NORMAL);
        
        const double distance = MortarUtilities::FastProjectDirection(OriginalMasterGeometry, OriginalSlaveGeometry[i_slave].Coordinates(), projected_gp_global, MasterNormal, -normal ); // The opposite direction
        
        if (distance > mDistanceThreshold) {
            ConditionsPointsSlave.clear();
            return false;
        }
        
        const bool is_inside = OriginalMasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local, tolerance );
        
        if (is_inside == true) {// The slave node belongs to the master
            if (i_slave == 0)
                auxiliar_coordinates[0] = -1.0;  // First node
            else
                auxiliar_coordinates[1] = 1.0;  // Second node
        }
    }

    // We check if the element is fully integrated
    if ((auxiliar_coordinates[0] == -1.0 && auxiliar_coordinates[1] == 1.0)) {
        total_weight = 2.0;
    } else  // If not then we proceed
    {
        std::vector<double> auxiliar_xi;
        for (unsigned int i_master = 0; i_master < 2; ++i_master) {
            projected_gp_local[0] = (i_master == 0) ? -1.0 : 1.0;
            double delta_xi = (i_master == 0) ? 0.5 : -0.5;
            const bool is_inside = MortarUtilities::ProjectIterativeLine2D(OriginalSlaveGeometry, OriginalMasterGeometry[i_master].Coordinates(), projected_gp_local, SlaveNormal, tolerance, delta_xi);
            
            if (is_inside == true)  auxiliar_xi.push_back(projected_gp_local[0]);
        }

        // In this case one edge of the slave belongs to the master and additionally one node of the master belongs to the slave
        if (auxiliar_xi.size() == 1 && ((auxiliar_coordinates[0] == - 1.0 || auxiliar_coordinates[1] == 1.0))) {
            if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) // NOTE: Equivalent to == -1.0
                auxiliar_coordinates[1] = auxiliar_xi[0];
            else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) // NOTE: Equivalent to == 1.0
                auxiliar_coordinates[0] = auxiliar_xi[0];
            else {
                KRATOS_WATCH(auxiliar_xi[0]);
                KRATOS_WATCH(auxiliar_coordinates[0]);
                KRATOS_WATCH(auxiliar_coordinates[1]);
                KRATOS_ERROR  << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!! (TYPE 0)" << std::endl;
            }
        } else if ( auxiliar_xi.size() == 2) { // Both nodes of the master belong to the slave (and none of the nodes of the slave belong to the master, the nodes can coincide, there is no other possibility)
            if (std::abs(auxiliar_coordinates[0] + 1.0) < tolerance) // NOTE: Equivalent to == -1.0. In this case the node in the left edge is already assigned 
                auxiliar_coordinates[1] = auxiliar_xi[0] < auxiliar_xi[1] ? auxiliar_xi[1] : auxiliar_xi[0]; // We set in the proper position
            else if (std::abs(auxiliar_coordinates[1] - 1.0) < tolerance) // NOTE: Equivalent to == 1.0. In this case the node in the right edge is already assigned 
                auxiliar_coordinates[0] = auxiliar_xi[0] < auxiliar_xi[1] ? auxiliar_xi[0] : auxiliar_xi[1]; // We set in the proper position
            else { // There isn't any coincidence with the edges
                if (auxiliar_xi[0] < auxiliar_xi[1]) { // We check that are in proper order
                    auxiliar_coordinates[0] = auxiliar_xi[0];
                    auxiliar_coordinates[1] = auxiliar_xi[1];
                } else {
                    auxiliar_coordinates[1] = auxiliar_xi[0];
                    auxiliar_coordinates[0] = auxiliar_xi[1];
                }
            }
        } else { // THIS IS NOT SUPPOSED TO HAPPEN
#ifdef KRATOS_DEBUG
            KRATOS_WATCH(OriginalSlaveGeometry);
            KRATOS_WATCH(OriginalMasterGeometry);
            KRATOS_ERROR
                << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!! (TYPE 1)"
                << std::endl;
#endif
            return false;  // NOTE: Giving problems
        }

        total_weight = auxiliar_coordinates[1] - auxiliar_coordinates[0];
    }
    
    KRATOS_ERROR_IF(total_weight < 0.0) << "WAAAAAAAAAAAAARNING!!!!!!!!, wrong order of the coordinates: "<< auxiliar_coordinates << std::endl;
    KRATOS_ERROR_IF(total_weight > 2.0) << "WAAAAAAAAAAAAARNING!!!!!!!!, impossible, Weight higher than 2: "<< auxiliar_coordinates << std::endl;

    // We do the final assignmen
    if (total_weight > std::numeric_limits<double>::epsilon()) {
        ConditionsPointsSlave.resize(1);
        array_1d<PointType, 2> list_points;
        list_points[0].Coordinates()[0] = auxiliar_coordinates[0];
        list_points[1].Coordinates()[0] = auxiliar_coordinates[1];
        ConditionsPointsSlave[0] = list_points;

        return true;
    } else {
        ConditionsPointsSlave.clear();
        return false;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 3, false>::GetExactIntegration(
    GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    ConditionArrayListType& ConditionsPointsSlave) {
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType& slave_center = OriginalSlaveGeometry.Center();

    // We define the condition tangents
    const array_1d<double, 3> slave_tangent_xi = (OriginalSlaveGeometry[1].Coordinates() - OriginalSlaveGeometry[0].Coordinates()) /
                                           norm_2(OriginalSlaveGeometry[1].Coordinates() -  OriginalSlaveGeometry[0].Coordinates());
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, SlaveNormal, slave_tangent_xi);

    // We define the tolerance
    const double tolerance = std::numeric_limits<double>::epsilon();

    // We define the auxiliar geometry
    std::vector<PointType::Pointer> points_array_slave(3);
    std::vector<PointType::Pointer> points_array_master(3);
    for (unsigned int i_node = 0; i_node < 3; ++i_node) {
        PointType aux_point;
        double distance;

        aux_point.Coordinates() = OriginalSlaveGeometry[i_node].Coordinates();  // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave[i_node] = boost::make_shared<PointType>(aux_point);
        
        aux_point = MortarUtilities::FastProject( slave_center, OriginalMasterGeometry[i_node], SlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master[i_node] = boost::make_shared<PointType>(aux_point);
        
        if (distance > mDistanceThreshold) {
            ConditionsPointsSlave.clear();
            return false;
        }
    }

    Triangle3D3<PointType> slave_geometry(points_array_slave);
    Triangle3D3<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 3> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, tolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside) == true) {
        ConditionsPointsSlave.resize(1);

        for (unsigned int i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            OriginalSlaveGeometry.PointLocalCoordinates( point, OriginalMasterGeometry[i_node]);
            ConditionsPointsSlave[0][i_node] = point;
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, tolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry);

        return TriangleIntersections<GeometryNodeType>(ConditionsPointsSlave, point_list, OriginalSlaveGeometry, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 4, false>::GetExactIntegration(
    GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    ConditionArrayListType& ConditionsPointsSlave) {
    // We define the tolerance
    const double tolerance = std::numeric_limits<double>::epsilon();

    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType& slave_center = OriginalSlaveGeometry.Center();

    // We define the condition tangents
    const array_1d<double, 3> slave_tangent_xi = (OriginalSlaveGeometry[2].Coordinates() - OriginalSlaveGeometry[0].Coordinates()) /
                                           norm_2(OriginalSlaveGeometry[2].Coordinates() - OriginalSlaveGeometry[0].Coordinates());
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, SlaveNormal, slave_tangent_xi);

    // We define the auxiliar geometry
    std::vector<PointType::Pointer> points_array_slave(4);
    std::vector<PointType::Pointer> points_array_slave_not_rotated(4);
    std::vector<PointType::Pointer> points_array_master(4);
    for (unsigned int i_node = 0; i_node < 4; ++i_node) {
        PointType aux_point;
        double distance_slave, distance_master;

        aux_point = MortarUtilities::FastProject(slave_center, OriginalSlaveGeometry[i_node], SlaveNormal, distance_slave);
        points_array_slave_not_rotated[i_node] = boost::make_shared<PointType>(aux_point);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave[i_node] = boost::make_shared<PointType>(aux_point);
        
        aux_point = MortarUtilities::FastProject( slave_center, OriginalMasterGeometry[i_node], SlaveNormal, distance_master);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master[i_node] = boost::make_shared<PointType>(aux_point);
        
        if (distance_slave > mDistanceThreshold || distance_master > mDistanceThreshold) {
            ConditionsPointsSlave.clear();
            return false;
        }
    }

    Quadrilateral3D4<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> slave_geometry_not_rotated(points_array_slave_not_rotated);
    Quadrilateral3D4<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 4> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, tolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside) == true) {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry);
        
        return TriangleIntersections<GeometryPointType>(ConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center, true);
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, tolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry);
        
        return TriangleIntersections<GeometryPointType>(ConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

// NOTE: The following are "hardcopies" of the previous ones, this is because C++ doesn't allow yet the partial specialization (CHECK THE ERROR TWICE!!!!!)

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<2, 2, true>::GetExactIntegration(
    GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    ConditionArrayListType& ConditionsPointsSlave) {
    // We take the geometry GP from the core
    const double tolerance = 1.0e-6;  // std::numeric_limits<double>::epsilon();

    double total_weight = 0.0;
    array_1d<double, 2> auxiliar_coordinates = ZeroVector(2);
    array_1d<PointBelongsLine2D2N, 2> auxiliar_belong;

    // Declaring auxiliar values
    PointType projected_gp_global;
    GeometryNodeType::CoordinatesArrayType projected_gp_local;

    // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
    for (unsigned int i_slave = 0; i_slave < 2; ++i_slave) {
        const array_1d<double, 3>& normal = OriginalSlaveGeometry[i_slave].FastGetSolutionStepValue(NORMAL);
        
        const double distance = MortarUtilities::FastProjectDirection(OriginalMasterGeometry, OriginalSlaveGeometry[i_slave].Coordinates(), projected_gp_global, MasterNormal, -normal ); // The opposite direction
        
        if (distance > mDistanceThreshold) {
            ConditionsPointsSlave.clear();
            return false;
        }
        
        const bool is_inside = OriginalMasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local, tolerance );
        
        if (is_inside == true) {// The slave node belongs to the master
            if (i_slave == 0) {  // First node
                auxiliar_coordinates[0] = -1.0;
                auxiliar_belong[0] = SlaveLine2D2N0;
            } else { // Second node
                auxiliar_coordinates[1] = 1.0;
                auxiliar_belong[1] = SlaveLine2D2N1;
            }
        }
    }

    // We check if the element is fully integrated
    if ((auxiliar_coordinates[0] == - 1.0 && auxiliar_coordinates[1] == 1.0) == true)
        total_weight = 2.0;
    else { // If not then we proceed
        std::vector<double> auxiliar_xi;
        std::vector<PointBelongsLine2D2N> auxiliar_master_belong;
        for (unsigned int i_master = 0; i_master < 2; ++i_master) {
            projected_gp_local[0] = (i_master == 0) ? -1.0 : 1.0;
            double delta_xi = (i_master == 0) ? 0.5 : -0.5;
            const bool is_inside = MortarUtilities::ProjectIterativeLine2D(OriginalSlaveGeometry,  OriginalMasterGeometry[i_master].Coordinates(), projected_gp_local, SlaveNormal, tolerance, delta_xi);

            if (is_inside == true) {
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
                KRATOS_WATCH(auxiliar_xi[0]);
                KRATOS_WATCH(auxiliar_coordinates[0]);
                KRATOS_WATCH(auxiliar_coordinates[1]);
                KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!! (TYPE 0)" << std::endl;
            }
        } else if (  auxiliar_xi.size() == 2) { // Both nodes of the master belong to the slave (and none of the nodes of the slave belong to the master, the nodes can coincide, there is no other possibility)
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
#ifdef KRATOS_DEBUG
            KRATOS_WATCH(OriginalSlaveGeometry);
            KRATOS_WATCH(OriginalMasterGeometry);
            KRATOS_ERROR << "WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!! (TYPE 1)" << std::endl;
#endif
            return false;  // NOTE: Giving problems
        }

        total_weight = auxiliar_coordinates[1] - auxiliar_coordinates[0];
    }
    
    KRATOS_ERROR_IF(total_weight < 0.0) << "WAAAAAAAAAAAAARNING!!!!!!!!, wrong order of the coordinates: "<< auxiliar_coordinates << std::endl;
    KRATOS_ERROR_IF(total_weight > 2.0) << "WAAAAAAAAAAAAARNING!!!!!!!!, impossible, Weight higher than 2: "<< auxiliar_coordinates << std::endl;

    // We do the final assignmen
    if (total_weight > std::numeric_limits<double>::epsilon()) {
        ConditionsPointsSlave.resize(1);
        array_1d<PointBelong<2>, 2> list_points;
        list_points[0].Coordinates()[0] = auxiliar_coordinates[0];
        list_points[0].SetBelong(auxiliar_belong[0]);
        list_points[1].Coordinates()[0] = auxiliar_coordinates[1];
        list_points[1].SetBelong(auxiliar_belong[1]);
        ConditionsPointsSlave[0] = list_points;

        return true;
    } else {
        ConditionsPointsSlave.clear();
        return false;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 3, true>::GetExactIntegration(
    GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    ConditionArrayListType& ConditionsPointsSlave) 
{
    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType& slave_center = OriginalSlaveGeometry.Center();

    // We define the condition tangents
    const array_1d<double, 3> slave_tangent_xi =  (OriginalSlaveGeometry[1].Coordinates() - OriginalSlaveGeometry[0].Coordinates()) /
                                            norm_2(OriginalSlaveGeometry[1].Coordinates() - OriginalSlaveGeometry[0].Coordinates());
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, SlaveNormal, slave_tangent_xi);

    // We define the tolerance
    const double tolerance = std::numeric_limits<double>::epsilon();

    // We define the auxiliar geometry
    std::vector<PointType::Pointer> points_array_slave(3);
    std::vector<PointType::Pointer> points_array_master(3);
    for (unsigned int i_node = 0; i_node < 3; ++i_node) {
        PointType aux_point;
        double distance;

        aux_point.Coordinates() = OriginalSlaveGeometry[i_node].Coordinates();  // NOTE: We are in a linear triangle, all the nodes belong already to the plane, so, the step one can be avoided, we directly project  the master nodes
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave[i_node] = boost::make_shared<PointType>(aux_point);

        aux_point = MortarUtilities::FastProject(slave_center, OriginalMasterGeometry[i_node], SlaveNormal, distance);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master[i_node] = boost::make_shared<PointType>(aux_point);
    }

    Triangle3D3<PointType> slave_geometry(points_array_slave);
    Triangle3D3<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 3> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, tolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside) == true) {
        ConditionsPointsSlave.resize(1);

        for (unsigned int i_node = 0; i_node < 3; ++i_node) {
            PointType point;
            OriginalSlaveGeometry.PointLocalCoordinates( point, OriginalMasterGeometry[i_node]);
            ConditionsPointsSlave[0][i_node] =  PointBelong<3>(point.Coordinates(), static_cast<PointBelongsTriangle3D3N>(i_node + 3));
        }

        return true;
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry, Master);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, tolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry, Slave);

        return TriangleIntersections<GeometryNodeType>(ConditionsPointsSlave, point_list, OriginalSlaveGeometry, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
bool ExactMortarIntegrationUtility<3, 4, true>::GetExactIntegration(
    GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    ConditionArrayListType& ConditionsPointsSlave) {
    // We define the tolerance
    const double tolerance = std::numeric_limits<double>::epsilon();

    // Firt we create an auxiliar plane based in the condition center and its normal
    const PointType& slave_center = OriginalSlaveGeometry.Center();

    // We define the condition tangents
    const array_1d<double, 3> slave_tangent_xi = (OriginalSlaveGeometry[2].Coordinates() - OriginalSlaveGeometry[0].Coordinates()) /
                                           norm_2(OriginalSlaveGeometry[2].Coordinates() - OriginalSlaveGeometry[0].Coordinates());
    array_1d<double, 3> slave_tangent_eta;
    MathUtils<double>::CrossProduct( slave_tangent_eta, SlaveNormal, slave_tangent_xi);

    // We define the auxiliar geometry
    std::vector<PointType::Pointer> points_array_slave(4);
    std::vector<PointType::Pointer> points_array_slave_not_rotated(4);
    std::vector<PointType::Pointer> points_array_master(4);
    for (unsigned int i_node = 0; i_node < 4; ++i_node) {
        PointType aux_point;
        double distance_slave, distance_master;

        aux_point = MortarUtilities::FastProject(slave_center, OriginalSlaveGeometry[i_node], SlaveNormal, distance_slave);
        points_array_slave_not_rotated[i_node] = boost::make_shared<PointType>(aux_point);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_slave[i_node] = boost::make_shared<PointType>(aux_point);

        aux_point = MortarUtilities::FastProject(slave_center, OriginalMasterGeometry[i_node], SlaveNormal, distance_master);
        MortarUtilities::RotatePoint(aux_point, slave_center, slave_tangent_xi, slave_tangent_eta, false);
        points_array_master[i_node] = boost::make_shared<PointType>(aux_point);
        
        if (distance_slave > mDistanceThreshold || distance_master > mDistanceThreshold) {
            ConditionsPointsSlave.clear();
            return false;
        }
    }

    Quadrilateral3D4<PointType> slave_geometry(points_array_slave);
    Quadrilateral3D4<PointType> slave_geometry_not_rotated(points_array_slave_not_rotated);
    Quadrilateral3D4<PointType> master_geometry(points_array_master);

    // No we project both nodes from the slave side and the master side
    array_1d<bool, 4> all_inside;

    // We check if the nodes are inside
    CheckInside(all_inside, slave_geometry, master_geometry, tolerance);

    // We create the pointlist
    PointListType point_list;

    // All the master points are inside the slave geometry
    if (CheckAllInside(all_inside) == true) {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry, Master);

        return TriangleIntersections<GeometryPointType>(ConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry,  master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center,  true);
    } else {
        // We add the internal nodes
        PushBackPoints(point_list, all_inside, master_geometry, Master);

        // We check if the nodes are inside
        CheckInside(all_inside, master_geometry, slave_geometry, tolerance);

        // We add the internal nodes
        PushBackPoints(point_list, all_inside, slave_geometry, Slave);

        return TriangleIntersections<GeometryPointType>(ConditionsPointsSlave, point_list, slave_geometry_not_rotated, slave_geometry, master_geometry, slave_tangent_xi, slave_tangent_eta, slave_center);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
bool ExactMortarIntegrationUtility<TDim, TNumNodes,
    TBelong>::GetExactIntegration(GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    IntegrationPointsType& IntegrationPointsSlave) {
    ConditionArrayListType conditions_points_slave;

    const bool is_inside = GetExactIntegration(OriginalSlaveGeometry, SlaveNormal, OriginalMasterGeometry, MasterNormal, conditions_points_slave);

    for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array(TDim);  // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (unsigned int i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            OriginalSlaveGeometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
            points_array[i_node] = boost::make_shared<PointType>(global_point);
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
            OriginalSlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);
            
            const double det_J = decomp_geom.DeterminantOfJacobian( local_point_decomp ) * (TDim == 2 ? 2.0 : 1.0);
            
            const array_1d<double, 3>& coordinates = local_point_parent.Coordinates();
            IntegrationPointsSlave.push_back( IntegrationPointType( coordinates[0], coordinates[1], weight * det_J )); // TODO: Change push_back for a fix operation
        }
    }

    return is_inside;
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
bool ExactMortarIntegrationUtility<TDim, TNumNodes,
    TBelong>::GetExactAreaIntegration(GeometryNodeType& OriginalSlaveGeometry,
    const array_1d<double, 3>& SlaveNormal,
    GeometryNodeType& OriginalMasterGeometry,
    const array_1d<double, 3>& MasterNormal, double& rArea) {
        
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = GetExactIntegration(OriginalSlaveGeometry, SlaveNormal, OriginalMasterGeometry, MasterNormal, conditions_points_slave);
    if (is_inside) GetTotalArea(OriginalSlaveGeometry, conditions_points_slave, rArea);
    
    return is_inside;
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::GetTotalArea(
    GeometryNodeType& OriginalSlaveGeometry,
    ConditionArrayListType& ConditionsPointsSlave, double& rArea) {
    
    rArea = 0.0;

    for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array( TDim);  // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (unsigned int i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            OriginalSlaveGeometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = boost::make_shared<PointType>(global_point);
        }
        
        DecompositionType decomp_geom( points_array );
        
        rArea += decomp_geom.Area();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::TestGetExactIntegration(
    Condition::Pointer& SlaveCond,
    Condition::Pointer& MasterCond, 
    Matrix& CustomSolution) {
    
    IntegrationPointsType integration_points_slave;

    const bool solution = GetExactIntegration(SlaveCond->GetGeometry(), SlaveCond->GetValue(NORMAL), MasterCond->GetGeometry(), MasterCond->GetValue(NORMAL), integration_points_slave);

    CustomSolution.resize(integration_points_slave.size(), TDim, false);
    
    for (unsigned int GP = 0; GP < integration_points_slave.size(); ++GP) {
        // Solution save:
        const array_1d<double, 3>& coordinates_gp = integration_points_slave[GP].Coordinates();
        CustomSolution(GP, 0) = coordinates_gp[0];
        if (TDim == 2)
            CustomSolution(GP, 1) = integration_points_slave[GP].Weight();
        else {
            CustomSolution(GP, 1) = coordinates_gp[1];
            CustomSolution(GP, 2) = integration_points_slave[GP].Weight();
        }
    }

    return solution;
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes, bool TBelong>
double ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::TestGetExactAreaIntegration(
    ModelPart& rMainModelPart,
    Condition::Pointer& SlaveCond
    )
{        
    // Initalize values
    double area = 0.0;
    IndexSet::Pointer indexes_set = SlaveCond->GetValue( INDEX_SET );
    
    for (auto it_pair = indexes_set->begin(); it_pair != indexes_set->end(); ++it_pair ) {
        double local_area = 0.0;
        Condition::Pointer p_master_cond = rMainModelPart.pGetCondition(*it_pair);
        const bool is_inside = GetExactAreaIntegration(SlaveCond->GetGeometry(), SlaveCond->GetValue(NORMAL), p_master_cond->GetGeometry(), p_master_cond->GetValue(NORMAL), local_area);
        if (is_inside) area += local_area;
    }
    
    return area;
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::GetIntegrationMethod() {
    // Setting the auxiliar integration points
    switch (mIntegrationOrder) {
        case 1: mAuxIntegrationMethod = GeometryData::GI_GAUSS_1;
        case 2: mAuxIntegrationMethod = GeometryData::GI_GAUSS_2;
        case 3: mAuxIntegrationMethod = GeometryData::GI_GAUSS_3;
        case 4: mAuxIntegrationMethod = GeometryData::GI_GAUSS_4;
        case 5: mAuxIntegrationMethod = GeometryData::GI_GAUSS_5;
        default: mAuxIntegrationMethod = GeometryData::GI_GAUSS_2;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
GeometryNodeType::IntegrationPointsArrayType ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::GetIntegrationTriangle() {
    // Setting the auxiliar integration points
    switch (mIntegrationOrder) {
        case 1: return Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 2: return Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 3: return Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 4: return Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        case 5: return Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        default: return Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
inline void
ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::PushBackPoints(
    VectorPoints& PointList, const array_1d<bool, TNumNodes>& AllInside,
    GeometryPointType& ThisGeometry) {
    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
        if (AllInside[i_node] == true) {
            // We check if the node already exists
            bool add_point = true;
            for (unsigned int iter = 0; iter < PointList.size(); ++iter)
                if (CheckPoints(ThisGeometry[i_node], PointList[iter])) add_point = false;
                    
            if (add_point == true) 
                PointList.push_back(ThisGeometry[i_node]);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
inline void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::PushBackPoints(
    VectorPointsBelong& PointList, 
    const array_1d<bool, TNumNodes>& AllInside,
    GeometryPointType& ThisGeometry, 
    const PointBelongs& ThisBelongs) {
    
    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
        if (AllInside[i_node] == true) {
            // We check if the node already exists
            bool add_point = true;
            for (unsigned int iter = 0; iter < PointList.size(); ++iter)
                if (CheckPoints(ThisGeometry[i_node], PointList[iter])) add_point = false;
                    
            if (add_point == true) {
                unsigned int initial_index = 0;
                if (ThisBelongs == Master)
                    initial_index = TNumNodes;
                PointList.push_back(
                    PointBelong<TNumNodes>(ThisGeometry[i_node].Coordinates(),
                        static_cast<BelongType>(initial_index + i_node)));
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
inline void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::CheckInside(
    array_1d<bool, TNumNodes>& AllInside, 
    GeometryPointType& Geometry1,
    GeometryPointType& Geometry2, 
    const double Tolerance) {
    
    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
        GeometryNodeType::CoordinatesArrayType projected_gp_local;
        AllInside[i_node] = Geometry1.IsInside( Geometry2[i_node].Coordinates( ), projected_gp_local, Tolerance) ;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
inline std::vector<std::size_t> ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::ComputeAnglesIndexes(PointListType& PointList) const {
    
    const unsigned int list_size = PointList.size();

    // We reorder the nodes according with the angle they form with the first node
    std::vector<double> angles(list_size - 1);
    array_1d<double, 3> v =
        PointList[1].Coordinates() - PointList[0].Coordinates();
    v /= norm_2(v);
    array_1d<double, 3> n = GetNormalVector2D(v);

    for (unsigned int elem = 1; elem < list_size; ++elem) {
        angles[elem - 1] = AnglePoints(PointList[0], PointList[elem], v, n);
        if (angles[elem - 1] < 0.0) {
            v = PointList[elem].Coordinates() - PointList[0].Coordinates();
            v /= norm_2(v);
            n = GetNormalVector2D(v);
            for (unsigned int aux_elem = 0; aux_elem <= (elem - 1); ++aux_elem)
                angles[aux_elem] -= angles[elem - 1];
        }
    }

    return MortarUtilities::SortIndexes<double>(angles);
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
inline void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::ComputeClippingIntersections(
    PointListType& PointList,
    GeometryPointType& Geometry1, 
    GeometryPointType& Geometry2,
    const PointType& RefCenter) {
    
    // We consider the Z coordinate constant
    const double z_ref = RefCenter.Coordinates()[2];

    // We find the intersection in each side
    for (unsigned int i_edge = 0; i_edge < TNumNodes; ++i_edge) {
        const unsigned int ip_edge = (i_edge == (TNumNodes - 1)) ? 0 : i_edge + 1;
        for (unsigned int j_edge = 0; j_edge < TNumNodes; ++j_edge) {
            const unsigned int jp_edge = (j_edge == (TNumNodes - 1)) ? 0 : j_edge + 1;

            PointType intersected_point;
            const bool intersected = Clipping2D( intersected_point, Geometry1[i_edge], Geometry1[ip_edge], Geometry2[j_edge], Geometry2[jp_edge] );
            
            if (intersected == true) {
                // Set the coordinate
                intersected_point.Coordinates()[2] = z_ref;

                // Ititialize the check
                bool add_point = true;
                for (unsigned int iter = 0; iter < PointList.size(); ++iter) {
                    if (CheckPoints(intersected_point, PointList[iter]) == true) {
                        add_point = false;
                        break;
                    }
                }

                if (add_point == true) {
                    if (TBelong == true) { // NOTE: We do some kind of strange hash to know the intersected edges
                        const unsigned int hash =  2 * TNumNodes + 10 * i_edge + 100 * ip_edge + 1000 * j_edge + 10000 * jp_edge;
                        PointList.push_back(PointBelong<TNumNodes>(intersected_point.Coordinates(), static_cast<BelongType>(hash)));
                    } else
                        PointList.push_back(intersected_point);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, bool TBelong>
template <class TGeometryType>
inline bool
ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong>::TriangleIntersections(
    ConditionArrayListType& ConditionsPointsSlave, PointListType& PointList,
    TGeometryType& OriginalSlaveGeometry, GeometryPointType& Geometry1,
    GeometryPointType& Geometry2, const array_1d<double, 3>& SlaveTangentXi,
    const array_1d<double, 3>& SlaveTangentEta, const PointType& RefCenter,
    const bool IsAllInside) {
    // We do the clipping
    if (IsAllInside == false)
        ComputeClippingIntersections(
            PointList, Geometry1, Geometry2, RefCenter);

    // We compose the triangles
    const unsigned int list_size = PointList.size();
    if (list_size >  2) { // Technically the minimum is three, just in case I consider 2
        const std::vector<std::size_t> index_vector = ComputeAnglesIndexes(PointList);

        ConditionsPointsSlave.resize((list_size - 2));

        // We recover this point to the triangle plane and compute the local coordinates
        for (unsigned int i_point_list = 0; i_point_list < PointList.size(); ++i_point_list) {
            MortarUtilities::RotatePoint(PointList[i_point_list], RefCenter, SlaveTangentXi, SlaveTangentEta, true);
            PointType local_point;
            OriginalSlaveGeometry.PointLocalCoordinates(local_point, PointList[i_point_list].Coordinates());
            PointList[i_point_list].Coordinates() = local_point.Coordinates();
        }

        for (unsigned int elem = 0; elem < list_size - 2; ++elem) { // NOTE: We always have two points less that the number of nodes
            ArrayTriangleType points_locals;

            const bool inverted_triangle = (FastTriagleCheck2D(PointList[0],  PointList[index_vector[elem] + 1], PointList[index_vector[elem + 1] + 1]) < 0.0);

            points_locals[(inverted_triangle == false) ? 0 : 2] = PointList[0];
            points_locals[1] = PointList[index_vector[elem + 0] + 1];
            points_locals[(inverted_triangle == true) ? 0 : 2] = PointList[index_vector[elem + 1] + 1];

            ConditionsPointsSlave[elem] = points_locals;
        }
        
        if (ConditionsPointsSlave.size() > 0)                
            return true;
        else
            return false;
    }
    else // No intersection
    {
        ConditionsPointsSlave.clear();
        return false;
    }

    ConditionsPointsSlave.clear();
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template class ExactMortarIntegrationUtility<2, 2, false>;
template class ExactMortarIntegrationUtility<3, 3, false>;
template class ExactMortarIntegrationUtility<3, 4, false>;

template class ExactMortarIntegrationUtility<2, 2, true>;
template class ExactMortarIntegrationUtility<3, 3, true>;
template class ExactMortarIntegrationUtility<3, 4, true>;

}  // namespace Kratos.
