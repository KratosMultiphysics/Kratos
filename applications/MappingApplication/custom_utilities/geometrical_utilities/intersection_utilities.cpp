//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "intersection_utilities.h"

namespace Kratos
{

void IntersectionUtilities::FindIntersection1DGeometries2D(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    // Only perform the test for 1D lines in 2D space. Maybe this is too strict...
    KRATOS_ERROR_IF(rModelPartDomainA.ConditionsBegin()->GetGeometry().LocalSpaceDimension() != 1 &&
        rModelPartDomainA.ConditionsBegin()->GetGeometry().Dimension() != 2)
        << "Can compare only line segments with other line segments." << std::endl;

    for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
        condition_a_itr != rModelPartDomainA.ConditionsEnd();
        ++condition_a_itr)
    {
        for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
            condition_b_itr != rModelPartDomainB.ConditionsEnd();
            ++condition_b_itr)
        {
            if (condition_a_itr->GetGeometry().HasIntersection(condition_b_itr->GetGeometry()))
            {
                rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                    condition_a_itr->pGetGeometry(), condition_b_itr->pGetGeometry()));
            }
        }
    }
}

void IntersectionUtilities::CreateQuadraturePointsCoupling1DGeometries2D(
    ModelPart& rModelPartCoupling,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    const ModelPart& rParentModelPart = *(rModelPartCoupling.GetParentModelPart());

    // Only perform the test for 1D lines in 2D space. Maybe this is too strict...
    KRATOS_ERROR_IF(rModelPartCoupling.GeometriesBegin()->LocalSpaceDimension() != 1 &&
        rModelPartCoupling.GeometriesBegin()->Dimension() != 2)
        << "Can compare only line segments with other line segments." << std::endl;

    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
        geometry_itr != rModelPartCoupling.GeometriesEnd();
        ++geometry_itr)
    {
        auto& r_geom_master = geometry_itr->GetGeometryPart(0);
        auto& r_geom_slave = geometry_itr->GetGeometryPart(1);

        CoordinatesArrayType local_parameter_1 = ZeroVector(3);
        CoordinatesArrayType local_parameter_2 = ZeroVector(3);

        r_geom_master.PointLocalCoordinates(local_parameter_1, r_geom_slave[0]);
        r_geom_master.PointLocalCoordinates(local_parameter_1, r_geom_slave[1]);

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
        // r_geom_slave.CreateQuadraturePointGeometries(quadrature_point_geometries_slave, 1, integration_points); TODO check if the line above is an OK replacement for this line

        // add the quadrature point geometry conditions to the result model part
        const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
            ? 1
            : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
        for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i) {
            // TODO check this approach
            //rModelPartResult.AddCondition(Kratos::make_intrusive<Condition>(
            //    id + i, CouplingGeometry<Node<3>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
            rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                id + i, CouplingGeometry<Node<3>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
        }
    }
}

} // namespace Kratos.
