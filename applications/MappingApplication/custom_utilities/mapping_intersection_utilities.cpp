//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

// System includes

// External includes

// Project includes
#include "mapping_intersection_utilities.h"

namespace Kratos
{

void MappingIntersectionUtilities::FindIntersection1DGeometries2D(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    KRATOS_ERROR_IF(rModelPartDomainA.ConditionsBegin()->GetGeometry().LocalSpaceDimension() != 1 &&
        rModelPartDomainA.ConditionsBegin()->GetGeometry().WorkingSpaceDimension() != 2)
        << "Can compare only line segments with other line segments." << std::endl;

    std::vector<array_1d<double, 3>> dummy;

    for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
        condition_a_itr != rModelPartDomainA.ConditionsEnd();
        ++condition_a_itr)
    {
        for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
            condition_b_itr != rModelPartDomainB.ConditionsEnd();
            ++condition_b_itr)
        {
            if (FindOverlapExtents1DGeometries2D(condition_a_itr->GetGeometry(), condition_b_itr->GetGeometry(), dummy))
            {
                rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                    condition_a_itr->pGetGeometry(), condition_b_itr->pGetGeometry()));
            }
        }
    }
}

void MappingIntersectionUtilities::CreateFEMFEMSurfaceCouplingGeometries(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    ModelPart& rModelPartResult)
{
    /*for (auto element_a_itr = rModelPartDomainA.ElementsBegin();
             element_a_itr != rModelPartDomainA.ElementsEnd();
             ++element_a_itr)
    {
        for (auto element_b_itr = rModelPartDomainB.ElementsBegin();
            element_b_itr != rModelPartDomainB.ElementsEnd();
            ++element_b_itr)
        {
            rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                element_a_itr->pGetGeometry(), element_b_itr->pGetGeometry()));
        }
    }*/

   for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
            condition_a_itr != rModelPartDomainA.ConditionsEnd();
             ++condition_a_itr)
    {
        for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
           condition_b_itr != rModelPartDomainB.ConditionsEnd();
            ++condition_b_itr)
        {
            rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
               condition_a_itr->pGetGeometry(),condition_b_itr->pGetGeometry()));
        }
    }

}

void MappingIntersectionUtilities::CreateQuadraturePointsCoupling1DGeometries2D(
    ModelPart& rModelPartCoupling,
    double Tolerance)
{
    const ModelPart& rParentModelPart = rModelPartCoupling.GetParentModelPart();

    KRATOS_ERROR_IF(rModelPartCoupling.GeometriesBegin()->LocalSpaceDimension() != 1 &&
        rModelPartCoupling.GeometriesBegin()->WorkingSpaceDimension() != 2)
        << "Can compare only line segments with other line segments." << std::endl;

    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
        geometry_itr != rModelPartCoupling.GeometriesEnd();
        ++geometry_itr)
    {
        auto& r_geom_master = geometry_itr->GetGeometryPart(0);
        auto& r_geom_slave = geometry_itr->GetGeometryPart(1);
        CoordinatesArrayType local_parameter_1 = ZeroVector(3);
        CoordinatesArrayType local_parameter_2 = ZeroVector(3);

        std::vector<array_1d<double, 3>> overlap_extents;
        KRATOS_ERROR_IF_NOT(MappingIntersectionUtilities::FindOverlapExtents1DGeometries2D(
            r_geom_master, r_geom_slave, overlap_extents, 1e-6))
            << "Lines do not intersect." << std::endl;
        r_geom_master.PointLocalCoordinates(local_parameter_1, overlap_extents[0]); // min of overlap
        r_geom_master.PointLocalCoordinates(local_parameter_2, overlap_extents[1]); // max of overlap


        const SizeType IntegrationPointsPerSpan = 2; // TODO this should depend on the basis order

        IntegrationPointsArrayType integration_points(IntegrationPointsPerSpan);

        typename IntegrationPointsArrayType::iterator integration_point_iterator = integration_points.begin();

        IntegrationPointUtilities::IntegrationPoints1D(
            integration_point_iterator,
            IntegrationPointsPerSpan,
            local_parameter_1[0], local_parameter_2[0]);

        // Determine quadrature point locations of span on master and then create quadrature point geometries
        GeometriesArrayType quadrature_point_geometries_master(IntegrationPointsPerSpan);
        CreateQuadraturePointsUtility<NodeType>::Create(
            r_geom_master, quadrature_point_geometries_master, integration_points, 1);

        // Transfer quadrature point locations to slave
        for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i)
        {
            CoordinatesArrayType local_parameter_slave = ZeroVector(3);
            r_geom_slave.PointLocalCoordinates(local_parameter_slave, quadrature_point_geometries_master[i].Center());

            integration_points[i].X() = local_parameter_slave[0];
        }

        // create slave quadrature point geometries
        GeometriesArrayType quadrature_point_geometries_slave(IntegrationPointsPerSpan);
        CreateQuadraturePointsUtility<NodeType>::Create(
            r_geom_slave, quadrature_point_geometries_slave, integration_points, 1);

        // add the quadrature point geometry conditions to the result model part
        const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
            ? 1
            : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
        for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i) {
            rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
        }
    }
}

void MappingIntersectionUtilities::CreateQuadraturePointsCoupling2DGeometries3D(
    ModelPart& rModelPartCoupling)
{
    const ModelPart& rParentModelPart = rModelPartCoupling.GetParentModelPart();


    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
        geometry_itr != rModelPartCoupling.GeometriesEnd();
        ++geometry_itr)
    {
        auto& r_geom_master = geometry_itr->GetGeometryPart(0);
        auto& r_geom_slave = geometry_itr->GetGeometryPart(1);
        
        IndexType number_of_nodes_master = r_geom_master.PointsNumber();
        IndexType number_of_nodes_slave = r_geom_master.PointsNumber();

        // Construct the clipper library objects
        Clipper2Lib::Paths64 master_element(1), slave_element(1), solution;
        const double factor = 1e-10;

        Clipper2Lib::Point64 int_point;
        int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
        int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());

        // Build the master element for clipper
        for (IndexType u = 0; u < number_of_nodes_master; ++u)
        {
            CoordinatesArrayType master_element_node = *(r_geom_master.pGetPoint(u));
            auto new_int_point = BrepTrimmingUtilities::ToIntPoint(master_element_node[0], master_element_node[1], factor);
            if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
            {
                master_element[0].push_back(new_int_point);
                int_point.x = new_int_point.x;
                int_point.y = new_int_point.y;
            }
        }

        // Build the slave element for clipper
        for (IndexType u = 0; u < number_of_nodes_slave; ++u)
        {
            CoordinatesArrayType slave_element_node = *(r_geom_slave.pGetPoint(u));
            auto new_int_point = BrepTrimmingUtilities::ToIntPoint(slave_element_node[0], slave_element_node[1], factor);
            if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
            {
                slave_element[0].push_back(new_int_point);
                int_point.x = new_int_point.x;
                int_point.y = new_int_point.y;
            }
        }

        // Intersect the triangle with the polygon resulting from the tessellation of the outer loop
        solution = Clipper2Lib::Intersect(master_element, slave_element, Clipper2Lib::FillRule::NonZero);

        if (solution.size() != 0)
        {
            std::vector<Matrix> triangles;
            BrepTrimmingUtilities::Triangulate_OPT(solution[0], triangles, factor);


            for (IndexType k = 0; k < triangles.size(); k++){
                // Find the position of the new triangles in the parameter space of the master element
                std::vector<CoordinatesArrayType> triangle_coordinates_master_parameter_space, triangle_master_physical_space;
                CoordinatesArrayType triangle_node_master_parameter_space = ZeroVector(3);
                CoordinatesArrayType triangle_node_master_physical_space = ZeroVector(3);


                 // Extract the triangles from the Triangulation_OPT method
                for (IndexType vertex_index = 0; vertex_index < 3; vertex_index++)
                {
                    triangle_node_master_physical_space[0] = triangles[k](vertex_index, 0);
                    triangle_node_master_physical_space[1] = triangles[k](vertex_index, 1);
                    triangle_node_master_physical_space[2] = 0.0;

                    triangle_master_physical_space.push_back(triangle_node_master_physical_space);
                }

                // Sort the triangle counter-clockwise
                std::vector<CoordinatesArrayType> sorted_triangle_master_physical_space;
                IgaMappingIntersectionUtilities::sortVerticesCounterClockwise(triangle_master_physical_space, sorted_triangle_master_physical_space);


                for (IndexType vertex_index = 0; vertex_index < 3; vertex_index++)
                {
                    r_geom_master.PointLocalCoordinates(triangle_node_master_parameter_space, sorted_triangle_master_physical_space[vertex_index]);
                    triangle_coordinates_master_parameter_space.push_back(triangle_node_master_parameter_space);
                }

                // Create integration points container for the master
                IntegrationPointsArrayType integration_points_master(3);
                typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();
                GeometriesArrayType quadrature_point_geometries_master(3);

                // Create integration points container for the slave
                IntegrationPointsArrayType integration_points_slave(3);
                GeometriesArrayType quadrature_point_geometries_slave(3);


                // Coordinates of the triangle vertex in the master parameter space
                double xi_0 = triangle_coordinates_master_parameter_space[0][0];
                double xi_1 = triangle_coordinates_master_parameter_space[1][0];
                double xi_2 = triangle_coordinates_master_parameter_space[2][0];
                double eta_0 = triangle_coordinates_master_parameter_space[0][1];
                double eta_1 = triangle_coordinates_master_parameter_space[1][1];
                double eta_2 = triangle_coordinates_master_parameter_space[2][1];

                // Creating the integration points
                IntegrationPointUtilities::IntegrationPointsTriangle2D(integration_point_master_iterator, 1, xi_0, xi_1, xi_2, eta_0, eta_1, eta_2);

                // Create master quadrature points
                CreateQuadraturePointsUtility<NodeType>::Create(
                        r_geom_master, quadrature_point_geometries_master, integration_points_master, 1);

                // Get the position of the master quadrature points in its physical space
                CoordinatesArrayType master_quadrature_point_xyz = ZeroVector(3);
                CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);
                integration_points_slave = integration_points_master;

                for (IndexType k = 0; k < quadrature_point_geometries_master.size(); k++)
                {
                    // Get the position of the master i_th quadrature point in its physical space
                    master_quadrature_point_xyz = quadrature_point_geometries_master[k].Center();

                    // Obtain the local coordinate of the i_th master quadrature point in the slave local space
                    r_geom_slave.PointLocalCoordinates(slave_quadrature_point_local_space, master_quadrature_point_xyz);

                    // Modify the slave integration points container
                    integration_points_slave[k].X() = slave_quadrature_point_local_space[0];
                    integration_points_slave[k].Y() = slave_quadrature_point_local_space[1];
                }

                // Create slave quadrature points
                CreateQuadraturePointsUtility<NodeType>::Create(
                        r_geom_slave, quadrature_point_geometries_slave, integration_points_slave, 1);

                // add the quadrature point geometry conditions to the result model part
                    const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                                             ? 1
                                             : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
        
                for (IndexType i = 0; i < quadrature_point_geometries_master.size(); ++i)
                {
                    rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                        id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
                }
            }
        } else
        {
            continue;
        }

    }
}

bool MappingIntersectionUtilities::FindOverlapExtents1DGeometries2D(
    const GeometryType& rMasterLine,
    const GeometryType& rSlaveLine,
    std::vector<array_1d<double, 3>>& rOverlapExtents,
    const double tolerance)
{
    // This finds the overlap extents, in global cartesian space
    // for intersection cases refer https://github.com/KratosMultiphysics/Kratos/pull/6785#issuecomment-614300457

    // Initial checks
    if (rOverlapExtents.size() != 2) rOverlapExtents.resize(2);
    KRATOS_ERROR_IF(rMasterLine.LocalSpaceDimension() != 1 || rSlaveLine.LocalSpaceDimension() != 1)
        << "Can only find overlap between two lines." << std::endl;

    // See if we have a point intersection
    const array_1d<double, 3> first_point = rMasterLine[0]; //p1
    const array_1d<double, 3> second_point = rMasterLine[1]; //p2
    const array_1d<double, 3> first_point_other = rSlaveLine[0]; //p3
    const array_1d<double, 3> second_point_other = rSlaveLine[1]; //p4
    const array_1d<double, 3> AB = second_point - first_point;
    const array_1d<double, 3> CD = second_point_other - first_point_other;
    const double master_length = inner_prod(AB, AB);
    const double slave_length = inner_prod(CD, CD);

    // parametric coordinate of intersection on current line
    const double numerator = ((first_point[0] - first_point_other[0]) * (first_point_other[1] - second_point_other[1]) - (first_point[1] - first_point_other[1]) * (first_point_other[0] - second_point_other[0]));
    const double denominator = ((first_point[0] - second_point[0]) * (first_point_other[1] - second_point_other[1]) - (first_point[1] - second_point[1]) * (first_point_other[0] - second_point_other[0]));

    // Filter out point overlap
    if (std::abs(denominator) > tolerance) {
        const double t = numerator / denominator;
        rOverlapExtents[0] = first_point + t * (second_point - first_point);
        rOverlapExtents[1] = rOverlapExtents[0];
        return false;
    }
    else
    {
        // We have parallel lines. Lines can only intersect if they are co-linear. Check this now.
        const double lhs = (second_point[1] - first_point[1]) * (first_point_other[0] - second_point[0]);
        const double rhs = (first_point_other[1] - second_point[1]) * (second_point[0] - first_point[0]);
        if (std::abs(lhs - rhs) < tolerance) { // Lines are parallel and co-linear, check that at least one point of the other line is within the local line
            const array_1d<double, 3> AC = first_point_other - first_point;
            const array_1d<double, 3> AD = second_point_other - first_point;
            const array_1d<double, 3> CA = first_point - first_point_other;
            const array_1d<double, 3> CB = second_point - first_point_other;
            const array_1d<double, 3> BD = second_point_other - second_point;

            if (inner_prod(AB, AC) >= (0.0 - tolerance) && inner_prod(AB, AC) <= (inner_prod(AB, AB) + tolerance)) {// Check if p3 is within the line
                if (std::abs(inner_prod(AB, AC)) <= tolerance) { // p3 coincides with p1. Check if the lines are in the same direction
                    if (inner_prod(AB, CD) > tolerance) { // case 4
                        rOverlapExtents[0] = first_point;
                        rOverlapExtents[1] = (master_length < slave_length)
                            ? second_point
                            : second_point_other;
                        return true;
                    } else return false;
                } else if (std::abs(inner_prod(AB, AC) - inner_prod(AB, AB)) < tolerance) { // p3 coincides with p2. Check if the lines are in the same direction
                    if (inner_prod(-1.0 * AB, CD) > tolerance) { // case 7
                        rOverlapExtents[1] = second_point;
                        rOverlapExtents[0] = (master_length < slave_length)
                            ? first_point
                            : second_point_other;
                        return true;
                    }
                    else return false;
                } else { // p3 lies within the line
                    if (inner_prod(CB,CD) > tolerance) { // points towards p2
                        rOverlapExtents[0] = first_point_other;
                        rOverlapExtents[1] = (inner_prod(CB,CB) < slave_length)
                            ? second_point
                            : second_point_other;
                    } else { // points towards p1
                        rOverlapExtents[1] = first_point_other;
                        rOverlapExtents[0] = (inner_prod(CA, CA) < slave_length)
                            ? first_point
                            : second_point_other;
                    }
                    return true;
                }
            }
            else if (inner_prod(AB, AD) >= (0.0 - tolerance) && inner_prod(AB, AD) <= (inner_prod(AB, AB) + tolerance)) { // Check if p4 is within the line
                if (std::abs(inner_prod(AB, AD)) <= tolerance) { // p4 coincides with Point 1. Check if the lines are in the same direction
                    if (inner_prod(AB, -1.0 * CD) > tolerance) {
                        rOverlapExtents[0] = first_point;
                        rOverlapExtents[1] = (master_length < slave_length)
                            ? second_point
                            : first_point_other;
                        return true;
                    } else return false;
                } else if (std::abs(inner_prod(AB, AD) - inner_prod(AB, AB)) < tolerance) { // p4 coincides with p2. Check if the lines are in the same direction
                    if (inner_prod(-1.0 * AB, -1.0 * CD) > tolerance) {
                        rOverlapExtents[1] = second_point;
                        rOverlapExtents[0] = (master_length < slave_length)
                            ? first_point
                            : first_point_other;
                        return true;
                    } else return false;
                } else {// p4 lies within the line
                    if (inner_prod(AD, CD) > tolerance) { // points towards p1
                        rOverlapExtents[1] = second_point_other;
                        rOverlapExtents[0] = (inner_prod(AD, AD) < slave_length)
                            ? first_point
                            : first_point_other;
                    }
                    else { // points towards p2
                        rOverlapExtents[0] = second_point_other;
                        rOverlapExtents[1] = (inner_prod(BD, BD) < slave_length)
                            ? second_point
                            : first_point_other;
                    }
                    return true;
                }
            }
            else if (inner_prod(CD, CA) > tolerance && inner_prod(-1.0 * CD, (second_point - second_point_other)) > tolerance) {// check if the line lies entirely within the other line
                rOverlapExtents[0] = first_point;
                rOverlapExtents[1] = second_point;
                return true; // the line lies entirely within the other line
            } else return false; // Lines are colinear, but do not overlap at all
        } else return false; // Lines are parallel but not colinear
    }
}

} // namespace Kratos.
