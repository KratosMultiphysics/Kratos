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
    KRATOS_ERROR_IF(rModelPartDomainA.ConditionsBegin()->GetGeometry().Dimension() != 1)
        << "Can compare only line segments with other line segments." << std::endl;

    for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
        condition_a_itr != rModelPartDomainA.ConditionsEnd();
        ++condition_a_itr)
    {
        for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
            condition_b_itr != rModelPartDomainB.ConditionsEnd();
            ++condition_b_itr)
        {
            // ERROR HasIntersection does not properly show parallel curves, which are the significant subset
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

        SizeType IntegrationPointsPerSpan = 2;

        IntegrationPointsArrayType integration_points(IntegrationPointsPerSpan);

        typename IntegrationPointsArrayType::iterator integration_point_iterator = integration_points.begin();

        IntegrationPointUtilities::IntegrationPoints1D(
            integration_point_iterator,
            IntegrationPointsPerSpan,
            local_parameter_1[0], local_parameter_2[0]);

        GeometriesArrayType quadrature_point_geometries_master(IntegrationPointsPerSpan);
        r_geom_master.CreateQuadraturePointGeometries(quadrature_point_geometries_master, 1, integration_points);

        for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i)
        {
            CoordinatesArrayType local_parameter_slave = ZeroVector(3);
            r_geom_slave.PointLocalCoordinates(local_parameter_slave, quadrature_point_geometries_master[i].Center());

            integration_points[i].X() = local_parameter_slave[0];
        }

        GeometriesArrayType quadrature_point_geometries_slave(IntegrationPointsPerSpan);
        r_geom_slave.CreateQuadraturePointGeometries(quadrature_point_geometries_slave, 1, integration_points);

        for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i)
        {
            IndexType id = 1;
            if (rModelPartResult.NumberOfConditions() > 0) {
                id = rModelPartResult.ConditionsEnd()->Id() + 1;
            }
            rModelPartResult.AddCondition(Kratos::make_intrusive<Condition>(
                id, CouplingGeometry<Node<3>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
        }
    }
}

} // namespace Kratos.
