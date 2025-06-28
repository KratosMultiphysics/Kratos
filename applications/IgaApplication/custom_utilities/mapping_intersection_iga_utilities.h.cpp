//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

// System includes

// External includes

// Project includes
#include "mapping_intersection_iga_utilities.h"

namespace Kratos
{

void MappingIntersectionIgaUtilities::FindIntersection1DGeometries2D(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    // KRATOS_ERROR_IF(rModelPartDomainA.ConditionsBegin()->GetGeometry().LocalSpaceDimension() != 1 &&
    //     rModelPartDomainA.ConditionsBegin()->GetGeometry().WorkingSpaceDimension() != 2)
    //     << "Can compare only line segments with other line segments." << std::endl;

    typename CouplingGeometry<NodeType>::GeometryPointerVector geometry_vector;

    KRATOS_WATCH(rModelPartDomainA)
    KRATOS_WATCH(rModelPartDomainB)
    KRATOS_WATCH(*(rModelPartDomainA.NodesBegin()))
    KRATOS_WATCH(*(rModelPartDomainB.NodesBegin()))

    // GeometryType& p_curve_on_surface_1 = ((rModelPartDomainA.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0);
    // GeometryType& p_curve_on_surface_2 = ((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0);

    // geometry_vector.push_back(p_curve_on_surface_1.pGetGeometryPart(std::numeric_limits<IndexType>::max() - 2));
    // geometry_vector.push_back(p_curve_on_surface_2.pGetGeometryPart(std::numeric_limits<IndexType>::max() - 2));

    geometry_vector.push_back(rModelPartDomainA.pGetGeometry(10));
    geometry_vector.push_back(rModelPartDomainB.pGetGeometry(10));

    auto p_coupling_geometry = Kratos::make_shared<CouplingGeometry<NodeType>>(
            geometry_vector);

    p_coupling_geometry->SetId(0);

    rModelPartResult.AddGeometry(p_coupling_geometry);
}

void MappingIntersectionIgaUtilities::CreateQuadraturePointsCoupling1DGeometries2D(
    ModelPart& rModelPartCoupling,
    double Tolerance)
{
    const ModelPart& rParentModelPart = rModelPartCoupling.GetParentModelPart();

    // KRATOS_ERROR_IF(rModelPartCoupling.GeometriesBegin()->LocalSpaceDimension() != 1 &&
    //     rModelPartCoupling.GeometriesBegin()->WorkingSpaceDimension() != 2)
    //     << "Can compare only line segments with other line segments." << std::endl;

    auto& geometry_itr = rModelPartCoupling.GetGeometry(0);

    IntegrationInfo integration_info = geometry_itr.GetDefaultIntegrationInfo();
    GeometriesArrayType quadrature_points;

    geometry_itr.CreateQuadraturePointGeometries(quadrature_points, 2, integration_info);

    // add the quadrature point geometry conditions to the result model part
    IndexType id = (rParentModelPart.NumberOfConditions() == 0)
        ? 1
        : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

    KRATOS_WATCH(quadrature_points.size())

    ModelPart::ConditionsContainerType new_condition_list;
    for (auto it = quadrature_points.ptr_begin(); it != quadrature_points.ptr_end(); ++it)
    {

        for (SizeType i = 0; i < (*it)->size(); ++i) {
            rModelPartCoupling.AddNode((*it)->pGetPoint(i));
        }

        rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(id, (*it)));
        id++;
    }

    KRATOS_WATCH(rModelPartCoupling)
}

void MappingIntersectionIgaUtilities::FindIntersection2DGeometries3D(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    // KRATOS_ERROR_IF(rModelPartDomainA.ConditionsBegin()->GetGeometry().LocalSpaceDimension() != 1 &&
    //     rModelPartDomainA.ConditionsBegin()->GetGeometry().WorkingSpaceDimension() != 2)
    //     << "Can compare only line segments with other line segments." << std::endl;

    typename CouplingGeometry<NodeType>::GeometryPointerVector geometry_vector;

    // GeometryType& p_curve_on_surface_1 = ((rModelPartDomainA.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0);
    // GeometryType& p_curve_on_surface_2 = ((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0);

    // geometry_vector.push_back(p_curve_on_surface_1.pGetGeometryPart(std::numeric_limits<IndexType>::max() - 2));
    // geometry_vector.push_back(p_curve_on_surface_2.pGetGeometryPart(std::numeric_limits<IndexType>::max() - 2));

    KRATOS_WATCH(rModelPartDomainA)
    KRATOS_WATCH(rModelPartDomainB)
    // KRATOS_WATCH(*(rModelPartDomainA.NodesBegin()))
    // KRATOS_WATCH(*(rModelPartDomainB.NodesBegin()))

    // KRATOS_WATCH(*((rModelPartDomainA.ConditionsBegin())->pGetGeometry())); //surface3D
    // KRATOS_WATCH(rModelPartDomainB.GetParentModelPart()) //coupling: why? // we need to assign the parent as well.
    // // KRATOS_WATCH(((rModelPartDomainB.GetParentModelPart().ConditionsBegin())->pGetGeometry())->GetGeometryParent(0)); //brep_surface
    // KRATOS_WATCH(*((rModelPartDomainB.ConditionsBegin())->pGetGeometry())); //quadrature point
    // KRATOS_WATCH(((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0)); //brep_surface
    // KRATOS_WATCH((((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0)).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)); //nurbs_surface
    // KRATOS_WATCH((((rModelPartDomainB.ConditionsEnd())->pGetGeometry())->GetGeometryParent(0)).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
    // exit(0);

    //hack for multipatch : TODO
    // IndexType temp = 0;
    // IndexType temp_surface_id = 0;
    // std::vector<int> ids;
    // for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
    //     condition_b_itr != rModelPartDomainB.ConditionsEnd();
    //     ++condition_b_itr)
    // {
    //     IndexType id = condition_b_itr->Id();

    //     if((id - temp) == 1)
    //     {
    //         temp = condition_b_itr->Id();
    //         if(temp_surface_id != ((condition_b_itr->pGetGeometry())->GetGeometryParent(0)).Id())
    //         {
    //             temp_surface_id = ((condition_b_itr->pGetGeometry())->GetGeometryParent(0)).Id();
    //             ids.push_back(condition_b_itr->Id());
    //         }
    //     }
    //     else
    //     {
    //         break;
    //     }
    // }
    // KRATOS_WATCH(ids)

    //find a surface
    // GeometryType& p_surface = ((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0); //brep_curve_on_surface 
    // auto p_surface_2 = rModelPartDomainB.GetRootModelPart().pGetGeometry(p_surface.Id());
    // auto p_brep_surface_2 = dynamic_pointer_cast<BrepSurfaceType>(p_surface_2);

    // KRATOS_WATCH(p_brep_surface_2)
    // KRATOS_WATCH(*p_brep_surface_2)

    for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
        condition_a_itr != rModelPartDomainA.ConditionsEnd();
        ++condition_a_itr)
    {
        // geometry_vector.push_back(condition_a_itr->pGetGeometry());
        // geometry_vector.push_back(((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0));

        //hack for multipatch : TODO
        // for(IndexType i = 0; i < ids.size(); i++)
        // {
            // rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
            //     condition_a_itr->pGetGeometry(), (((rModelPartDomainB.pGetCondition(ids[i]))->pGetGeometry())->GetGeometryParent(0)).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)));
        // }
        rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                condition_a_itr->pGetGeometry(), (((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0)).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)));
    }

   

    // auto p_coupling_geometry = Kratos::make_shared<CouplingGeometry<NodeType>>(
    //         geometry_vector);

    // p_coupling_geometry->SetId(0);
    // rModelPartResult.AddGeometry(p_coupling_geometry);

    KRATOS_WATCH(rModelPartResult)
    // exit(0);

    // CHECKPOINT!

    // const auto gp_canonical_tri = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3>>::GenerateIntegrationPoints();

    // //Master: FEM
    // IntegrationPointsArrayType integration_points_master(3*rModelPartDomainA.NumberOfConditions());
    // typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();

    // GeometriesArrayType quadrature_point_geometries_master(3*rModelPartDomainA.NumberOfConditions());

    // for (auto& r_cond : rModelPartDomainA.Conditions()) {

    //     std::vector<IndexType> Ids(r_cond.GetGeometry().size());
    //     std::vector<array_1d<double, 3>> rPointsCoordinates;

    //     // std::vector<std::size_t> equation_ids_vector;
    //     // const ProcessInfo& r_current_process_info = rModelPartDomainA.GetProcessInfo();
    //     // r_elem.EquationIdVector(equation_ids_vector,r_current_process_info);

    //     for (SizeType i = 0; i < r_cond.GetGeometry().size(); i++)
    //     {
    //         array_1d<double, 3> r_xyz = r_cond.GetGeometry().GetPoint( i ).GetInitialPosition();
    //         Ids.push_back(r_cond.GetGeometry().GetPoint( i ).Id());

    //         array_1d<double, 3> point_local_coords = ZeroVector(3);
    //         ((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0).ProjectionPointGlobalToLocalSpace(r_xyz,point_local_coords); //

    //         rPointsCoordinates.push_back(point_local_coords);
    //     }

    //     // const SizeType IntegrationPointsPerSpan = 2; // TODO this should depend on the basis order
    //     // IntegrationPointsArrayType integration_points(IntegrationPointsPerSpan);

    //     array_1d<double,3> r_local_coordinate = ZeroVector(3);
    //     for (IndexType gp_i = 0; gp_i < 3; ++gp_i) {
    //         // Calculate current Gauss point local coordinates
    //         const auto aux_gp = gp_canonical_tri[gp_i];

    //         r_local_coordinate[0] = rPointsCoordinates[0][0] * (1 - aux_gp[0] - aux_gp[1]) +
    //                               rPointsCoordinates[1][0] * aux_gp[0] +
    //                               rPointsCoordinates[2][0] * aux_gp[1];
    //         r_local_coordinate[1] = rPointsCoordinates[0][1] * (1 - aux_gp[0] - aux_gp[1]) +
    //                               rPointsCoordinates[1][1] * aux_gp[0] +
    //                               rPointsCoordinates[2][1] * aux_gp[1];

    //         const double area = (rPointsCoordinates[0][0] * (rPointsCoordinates[1][1] - rPointsCoordinates[2][1]) + rPointsCoordinates[1][0] * (rPointsCoordinates[2][1] - rPointsCoordinates[0][1]) + rPointsCoordinates[2][0] * (rPointsCoordinates[0][1] - rPointsCoordinates[1][1])) / 2;

    //         (*integration_point_master_iterator)[0] = r_local_coordinate[0];
    //         (*integration_point_master_iterator)[1] = r_local_coordinate[1];
    //         (*integration_point_master_iterator).Weight() = area / 3.0; //TODO

    //         integration_point_master_iterator++;
    //     }
    // }
    // KRATOS_WATCH(integration_points_master.size())

    // //up untl this point, we obtain the projection of gauss points from FEM mesh on the brep surface of IGA
    // //the general idea, to create coupling geometry based on these gauss point informations of both domains.

    // CreateQuadraturePointsUtility<NodeType>::Create(
    //         rModelPartDomainA.GetGeometry(0), quadrature_point_geometries_master, integration_points_master, 1);

    // //Slave: IGA :
    // GeometriesArrayType quadrature_point_geometries_slave(3*rModelPartDomainA.NumberOfConditions());
    // CreateQuadraturePointsUtility<NodeType>::Create(
    //         ((rModelPartDomainB.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0), quadrature_point_geometries_slave, integration_points_master, 1);


    // // add the quadrature point geometry conditions to the result model part
    // const ModelPart& rParentModelPart = rModelPartResult.GetParentModelPart();
    // const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
    //     ? 1
    //     : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
    // for (IndexType i = 0; i < 3*rModelPartDomainA.NumberOfConditions(); ++i) {
    //     rModelPartResult.AddCondition(Kratos::make_intrusive<Condition>(
    //         id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i)))); //different in IGA case
    // }

}

void MappingIntersectionIgaUtilities::CreateQuadraturePointsCoupling2DGeometries3D(
    ModelPart& rModelPartCoupling,
    double Tolerance)
{
    const auto gp_canonical_tri = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3>>::GenerateIntegrationPoints();

    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
        geometry_itr != rModelPartCoupling.GeometriesEnd();
        ++geometry_itr)
    {
        auto& r_geom_master = geometry_itr->GetGeometryPart(0);
        auto& r_geom_slave = geometry_itr->GetGeometryPart(1);

        // KRATOS_WATCH(r_geom_master) //2 dimensional triangle with three nodes in 3D space
        // KRATOS_WATCH(r_geom_slave) //3 dimensional nurbs surface.

        //Master: FEM
        IntegrationPointsArrayType integration_points_master(3);
        typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();
        GeometriesArrayType quadrature_point_geometries_master(3); 

        //Slave: IGA
        IntegrationPointsArrayType integration_points_slave(3);
        typename IntegrationPointsArrayType::iterator integration_point_slave_iterator = integration_points_slave.begin();
        GeometriesArrayType quadrature_point_geometries_slave(3);

        //Projection of the finite element mesh on the NURBS surface and to its parameter space
        std::vector<array_1d<double, 3>> rPointsCoordinatesMasters;
        std::vector<array_1d<double, 3>> rPointsCoordinatesSlaves;
        for (SizeType i = 0; i < r_geom_master.size(); i++)
        {
            array_1d<double, 3> r_xyz = r_geom_master.GetPoint( i ).GetInitialPosition();
            array_1d<double, 3> point_local_coords_master = ZeroVector(3);
            array_1d<double, 3> point_local_coords_slave = ZeroVector(3);
            r_geom_master.ProjectionPointGlobalToLocalSpace(r_xyz,point_local_coords_master); 
            r_geom_slave.ProjectionPointGlobalToLocalSpace(r_xyz,point_local_coords_slave); 
            rPointsCoordinatesMasters.push_back(point_local_coords_master);
            rPointsCoordinatesSlaves.push_back(point_local_coords_slave);
        }

        //Assign GPs on the NURBS surface (master GP)
        array_1d<double,3> r_local_coordinate_master = ZeroVector(3);
        for (IndexType gp_i = 0; gp_i < 3; ++gp_i) {
            // Calculate current Gauss point local coordinates
            const auto aux_gp = gp_canonical_tri[gp_i];

            r_local_coordinate_master[0] = rPointsCoordinatesMasters[0][0] * (1 - aux_gp[0] - aux_gp[1]) +
                                  rPointsCoordinatesMasters[1][0] * aux_gp[0] +
                                  rPointsCoordinatesMasters[2][0] * aux_gp[1];
            r_local_coordinate_master[1] = rPointsCoordinatesMasters[0][1] * (1 - aux_gp[0] - aux_gp[1]) +
                                  rPointsCoordinatesMasters[1][1] * aux_gp[0] +
                                  rPointsCoordinatesMasters[2][1] * aux_gp[1];

            const double area = (rPointsCoordinatesMasters[0][0] * (rPointsCoordinatesMasters[1][1] - rPointsCoordinatesMasters[2][1]) + rPointsCoordinatesMasters[1][0] * (rPointsCoordinatesMasters[2][1] - rPointsCoordinatesMasters[0][1]) + rPointsCoordinatesMasters[2][0] * (rPointsCoordinatesMasters[0][1] - rPointsCoordinatesMasters[1][1])) / 2;

            (*integration_point_master_iterator)[0] = r_local_coordinate_master[0];
            (*integration_point_master_iterator)[1] = r_local_coordinate_master[1];
            (*integration_point_master_iterator).Weight() = area / 3.0; 

            integration_point_master_iterator++;
        }

        //Assign GPs on the NURBS surface (slave GP)
        array_1d<double,3> r_local_coordinate_slave = ZeroVector(3);
        for (IndexType gp_i = 0; gp_i < 3; ++gp_i) {
            // Calculate current Gauss point local coordinates
            const auto aux_gp = gp_canonical_tri[gp_i];

            r_local_coordinate_slave[0] = rPointsCoordinatesSlaves[0][0] * (1 - aux_gp[0] - aux_gp[1]) +
                                  rPointsCoordinatesSlaves[1][0] * aux_gp[0] +
                                  rPointsCoordinatesSlaves[2][0] * aux_gp[1];
            r_local_coordinate_slave[1] = rPointsCoordinatesSlaves[0][1] * (1 - aux_gp[0] - aux_gp[1]) +
                                  rPointsCoordinatesSlaves[1][1] * aux_gp[0] +
                                  rPointsCoordinatesSlaves[2][1] * aux_gp[1];

            const double area = (rPointsCoordinatesSlaves[0][0] * (rPointsCoordinatesSlaves[1][1] - rPointsCoordinatesSlaves[2][1]) + rPointsCoordinatesSlaves[1][0] * (rPointsCoordinatesSlaves[2][1] - rPointsCoordinatesSlaves[0][1]) + rPointsCoordinatesSlaves[2][0] * (rPointsCoordinatesSlaves[0][1] - rPointsCoordinatesSlaves[1][1])) / 2;

            (*integration_point_slave_iterator)[0] = r_local_coordinate_slave[0];
            (*integration_point_slave_iterator)[1] = r_local_coordinate_slave[1];
            (*integration_point_slave_iterator).Weight() = area / 3.0; 

            integration_point_slave_iterator++;
        }

        //Slave: IGA :
        // CreateQuadraturePointsUtility<NodeType>::Create(
        //     r_geom_slave, quadrature_point_geometries_slave, integration_points_master, 1); 
        IntegrationInfo integration_info = r_geom_slave.GetDefaultIntegrationInfo();
        r_geom_slave.CreateQuadraturePointGeometries(quadrature_point_geometries_slave, 2, integration_points_slave, integration_info);

        //Master: FE mesh : (TO DO)
        CreateQuadraturePointsUtility<NodeType>::Create(
            r_geom_master, quadrature_point_geometries_master, integration_points_master, 1); //something is wrong here

        // KRATOS_WATCH(quadrature_point_geometries_master[1])
        // KRATOS_WATCH(quadrature_point_geometries_slave[1])

        // for (auto it = quadrature_point_geometries_slave.ptr_begin(); it != quadrature_point_geometries_slave.ptr_end(); ++it)
        // {
        //     for (SizeType i = 0; i < (*it)->size(); ++i) {
        //         KRATOS_WATCH((*it)->pGetPoint(i));
        //     }
        // }


        // add the quadrature point geometry conditions to the result model part
        const ModelPart& rParentModelPart = rModelPartCoupling.GetParentModelPart();
        const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
            ? 1
            : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
        for (IndexType i = 0; i < 3; ++i) {
            rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i)))); 
        }
    }

    KRATOS_WATCH(rModelPartCoupling)
}

} // namespace Kratos.
