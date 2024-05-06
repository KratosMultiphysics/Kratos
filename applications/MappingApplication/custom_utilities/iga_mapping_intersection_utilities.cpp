//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti
//                   Andrea Gorgi
//

// System includes

// External includes

// Project includes
#include "iga_mapping_intersection_utilities.h"

namespace Kratos
{

void IgaMappingIntersectionUtilities::IgaCreateBrepCurveOnSurfaceCouplingGeometries(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
        condition_a_itr != rModelPartDomainA.ConditionsEnd();
        ++condition_a_itr)
    {
        for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
            condition_b_itr != rModelPartDomainB.ConditionsEnd();
            ++condition_b_itr)
        {
            rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                condition_a_itr->pGetGeometry(), condition_b_itr->pGetGeometry()));
        }
    }
}

void IgaMappingIntersectionUtilities::IgaCreateQuadraturePointsCoupling1DGeometries2D(
    ModelPart& rModelPartCoupling,
    double Tolerance)
{
    const ModelPart& rParentModelPart = rModelPartCoupling.GetParentModelPart();

    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
        geometry_itr != rModelPartCoupling.GeometriesEnd();
        ++geometry_itr)
    {

        IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = geometry_itr->GetDefaultIntegrationInfo();

        geometry_itr->CreateIntegrationPoints(integration_points, integration_info);

        GeometriesArrayType master_and_slave_quadrature_points_geometries(integration_points.size()); // Vector of coupling geometries which stores the master and slave quedrature point geometries
        IndexType number_of_shape_functions_derivatives=3;

        if (integration_points.size() != 0) {
            geometry_itr->CreateQuadraturePointGeometries(master_and_slave_quadrature_points_geometries,number_of_shape_functions_derivatives,integration_points,integration_info);

            const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                ? 1
                : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
            
            const SizeType IntegrationPointsPerSpan = integration_info.GetNumberOfIntegrationPointsPerSpan(0);

            for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i) {
                rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                    id + i, master_and_slave_quadrature_points_geometries(i)));
            }
        }
        
    }
}

} // namespace Kratos.
