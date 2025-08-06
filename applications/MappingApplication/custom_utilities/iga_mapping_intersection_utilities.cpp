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
//

// System includes

// External includes

// Project includes
#include "iga_mapping_intersection_utilities.h"

namespace Kratos
{

void IgaMappingIntersectionUtilities::CreateIgaFEMCouplingGeometries(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    const bool& rIsOriginIga,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    const auto& r_iga_model_part = rIsOriginIga ? rModelPartDomainA : rModelPartDomainB;
    const auto& r_fem_model_part = rIsOriginIga ? rModelPartDomainB : rModelPartDomainA;

    for (auto& fem_cond : r_fem_model_part.Conditions()) {
        const auto& line_geom = fem_cond.GetGeometry();
        const Point line_center = line_geom.Center();
        
        double min_distance = std::numeric_limits<double>::max();
        Condition::Pointer p_closest_iga_condition = nullptr;

        // Get the closest condition on the IGA interface
        for (auto& iga_cond : r_iga_model_part.Conditions()) {
            const auto& brep_curve_on_surface_geom = iga_cond.GetGeometry();

            CoordinatesArrayType local_coords = ZeroVector(3);
            CoordinatesArrayType projected_point = ZeroVector(3);

            bool success = brep_curve_on_surface_geom.ProjectionPointGlobalToLocalSpace(line_center, local_coords, Tolerance);
            if (!success) continue;

            brep_curve_on_surface_geom.GlobalCoordinates(projected_point, local_coords);
            double distance = norm_2(line_center - projected_point);

            if (distance < min_distance) {
                min_distance = distance;
                p_closest_iga_condition = r_iga_model_part.pGetCondition(iga_cond.Id());
            }
        }

        // Create a coupling geometry relating both sides
        if (p_closest_iga_condition) {
            rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                p_closest_iga_condition->pGetGeometry(),
                fem_cond.pGetGeometry()
            ));
        }
    }
}

void IgaMappingIntersectionUtilities::CreateIgaFEMQuadraturePointsCouplingInterface(
    ModelPart& rModelPartCoupling,
    double Tolerance)
{
   const ModelPart &rParentModelPart = rModelPartCoupling.GetParentModelPart();

   // Loop over the coupling geometries and for each one create the integration points on the interface
   for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
             geometry_itr != rModelPartCoupling.GeometriesEnd();
             ++geometry_itr)
    {
        IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = geometry_itr->GetDefaultIntegrationInfo();

        geometry_itr->CreateIntegrationPoints(integration_points, integration_info);

        // Vector of coupling geometries which stores the origin and destination quadrature point geometries
        GeometriesArrayType origin_and_destination_quadrature_points_geometries(integration_points.size());
        IndexType number_of_shape_functions_derivatives = 3;

        if (integration_points.size() != 0){
            geometry_itr->CreateQuadraturePointGeometries(origin_and_destination_quadrature_points_geometries, number_of_shape_functions_derivatives, integration_points, integration_info);

            const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                                    ? 1
                                    : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

            for (IndexType i = 0; i < origin_and_destination_quadrature_points_geometries.size(); ++i)
            {
                rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                    id + i, origin_and_destination_quadrature_points_geometries(i)));
            }
        }
    }
}

} // namespace Kratos.
