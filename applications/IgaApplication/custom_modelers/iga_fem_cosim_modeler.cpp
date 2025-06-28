//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "iga_fem_cosim_modeler.h"

namespace Kratos
{
    ///@name Stages
    ///@{

    void IgaFemCosimModeler::SetupGeometryModel()
    {
        CheckParameters();

        // Create the coupling modelling part in the origin model
        ModelPart& coupling_model_part = (mpModels[0]->HasModelPart("coupling"))
            ? mpModels[0]->GetModelPart("coupling")
            : mpModels[0]->CreateModelPart("coupling");

        // Get the origin and destination interface names
        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;
        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
        }
        else
        {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
        }

        // create coupling conditions on interface depending on the dimension
        // const IndexType dim = 2;
        // if (dim == 2)
        // {
            // CreateInterfaceSurfaceCouplingConditions(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name), true); //FEM
            // CreateInterfaceSurfaceCouplingConditions(mpModels.back()->GetModelPart(destination_interface_sub_model_part_name), false); //IGA
        // }
        // else
        // {
        //     KRATOS_ERROR << "Not implemented yet" << std::endl;
        // }

        KRATOS_WATCH(origin_interface_sub_model_part_name)
        KRATOS_WATCH(destination_interface_sub_model_part_name)

        // KRATOS_WATCH(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name))
        // KRATOS_WATCH(*((mpModels[0]->GetModelPart(origin_interface_sub_model_part_name))).NodesBegin())

        // KRATOS_WATCH(mpModels[0]->GetModelPart(destination_interface_sub_model_part_name))
        // KRATOS_WATCH(*((mpModels[0]->GetModelPart(destination_interface_sub_model_part_name))).NodesBegin())

        

        // Transfer everything into the coupling modelpart
        //FEM
        ModelPart& coupling_interface_origin = (coupling_model_part.HasSubModelPart("interface_origin"))
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");
        CopySubModelPart(coupling_interface_origin,
            mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));


        //IGA
        ModelPart& coupling_interface_destination = (coupling_model_part.HasSubModelPart("interface_destination"))
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");
        CopySubModelPart(coupling_interface_destination,
            mpModels[1]->GetModelPart(destination_interface_sub_model_part_name));

        // KRATOS_WATCH((mpModels[1]->GetModelPart(destination_interface_sub_model_part_name)).GetParentModelPart()) 
        // KRATOS_WATCH(((mpModels[1]->GetModelPart(destination_interface_sub_model_part_name)).GetParentModelPart()).GetParentModelPart())
        // KRATOS_WATCH(*(((mpModels[1]->GetModelPart(destination_interface_sub_model_part_name)).ConditionsBegin())->pGetGeometry())); //quadrature point
        // KRATOS_WATCH(*((((mpModels[1]->GetModelPart(destination_interface_sub_model_part_name)).GetParentModelPart()).ConditionsBegin())->pGetGeometry())); 


        // KRATOS_WATCH(coupling_interface_destination)
        // exit(0);

        // KRATOS_ERROR_IF(coupling_interface_origin.NumberOfConditions() == 0)
        //     << "Coupling geometries are currently determined by conditions in the coupling sub model parts,"
        //     << " but there are currently not conditions in the coupling interface origin sub model part. Please specify some."
        //     << std::endl;
        // SizeType working_dim = coupling_interface_origin.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
        // SizeType local_dim = coupling_interface_origin.ConditionsBegin()->GetGeometry().LocalSpaceDimension();

        // if (working_dim == 3 && local_dim == 2) //IGA-IGA mapping
        // {
            // Divide overlapping interface into segments

            // Create quadrature coupling conditions from correponding integration points within segments
            MappingIntersectionIgaUtilities::FindIntersection2DGeometries3D(
                coupling_interface_origin,
                coupling_interface_destination,
                coupling_model_part, 1e-6);

            // add the quadrature point geometry conditions to the result model part
            MappingIntersectionIgaUtilities::CreateQuadraturePointsCoupling2DGeometries3D(
                coupling_model_part, 1e-6);

            //const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
            //    ? 1
            //    : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
            //for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i) {
            //    rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
            //        id + i, Kratos::make_shared<CouplingGeometry<Node>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
            //}
        // }
        // else
        // {
        //     KRATOS_ERROR << "Creation of coupling quadrature points not yet supported for requested"
        //         << " working space dimension = " << working_dim << " and local space dimension = "
        //         << local_dim << std::endl;
        // }

        KRATOS_WATCH(coupling_model_part)
    }


    void IgaFemCosimModeler::CheckParameters()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("is_interface_sub_model_parts_specified"))
            << "Missing \"is_interface_sub_model_parts_specified\" in MappingGeometriesModeler Parameters." << std::endl;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            KRATOS_ERROR_IF_NOT(mParameters.Has("origin_interface_sub_model_part_name"))
                << "Missing \"origin_interface_sub_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;

            KRATOS_ERROR_IF_NOT(mParameters.Has("destination_interface_sub_model_part_name"))
                << "Missing \"destination_interface_sub_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;
        }
    }

    void IgaFemCosimModeler::CreateInterfaceSurfaceCouplingConditions(ModelPart& rInterfaceModelPart, bool MasterToSlave)
    {
        rInterfaceModelPart.CreateSubModelPart("coupling_conditions");
        ModelPart& coupling_conditions = rInterfaceModelPart.GetSubModelPart("coupling_conditions");
        const ModelPart& root_mp = rInterfaceModelPart.GetRootModelPart();

        GeometryType& p_curve_on_surface = ((rInterfaceModelPart.ConditionsBegin())->pGetGeometry())->GetGeometryParent(0); //brep_curve_on_surface
        auto p_curve_on_surface_2 = root_mp.pGetGeometry(p_curve_on_surface.Id());
        auto p_brep_curve_on_surface_2 = dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_curve_on_surface_2);
        auto p_curve = p_brep_curve_on_surface_2->pGetCurveOnSurface()->pGetCurve();
        
        std::vector<double> spans;
        p_brep_curve_on_surface_2->pGetCurveOnSurface()->SpansLocalSpace(spans);

        // IndexType interface_node_id;
        rInterfaceModelPart.CreateSubModelPart("nodes_mp");
        ModelPart& nodes_mp = rInterfaceModelPart.GetSubModelPart("nodes_mp");

        PointerVector<Node> points;
        std::vector< GeometryPointerType> p_geom_vec;

        std::vector<int> nodes_id;

        for (auto& r_cond : rInterfaceModelPart.Conditions()) {
            
            auto& r_geom = r_cond.GetGeometry();
            auto& r_N = r_geom.ShapeFunctionsValues();

            for (IndexType i = 0; i<r_N.size2();++i)
            {
                if(r_N(0,i) > 1e-6)
                {
                    nodes_id.push_back(r_geom.pGetPoint(i)->Id());
                }
            }
        }

        std::sort(nodes_id.begin(), nodes_id.end());
        auto it = std::unique(nodes_id.begin(), nodes_id.end());
        nodes_id.erase(it, nodes_id.end());

        KRATOS_WATCH(nodes_id)

        for (size_t i = 0; i < nodes_id.size(); ++i)
        {
            points.push_back(rInterfaceModelPart.pGetNode(nodes_id[i]));
        }

        int p = points.size() + 1 - spans.size();
        Vector knot_vector = ZeroVector(points.size() + p - 1);
        for (IndexType i = p - 1; i < knot_vector.size(); ++i) {
            knot_vector(i) = spans[i - p + 1];
        }
        knot_vector(knot_vector.size()-1) = spans[spans.size()-1];

        auto curve = Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(points, p, knot_vector);
        curve->SetId(10);

        coupling_conditions.AddGeometry(curve);

        IntegrationInfo integration_info = curve->GetDefaultIntegrationInfo();
        GeometriesArrayType quadrature_points;

        curve->CreateQuadraturePointGeometries(quadrature_points, 2, integration_info);

        const ModelPart& rParentModelPart = rInterfaceModelPart.GetParentModelPart();
        IndexType id = (rParentModelPart.NumberOfConditions() == 0)
            ? 1
            : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

        for (auto it = quadrature_points.ptr_begin(); it != quadrature_points.ptr_end(); ++it)
        {
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                coupling_conditions.AssignNode((*it)->pGetPoint(i));
                // KRATOS_WATCH((*it)->pGetPoint(i))
            }
            KRATOS_WATCH(id)
            coupling_conditions.AddCondition(Kratos::make_intrusive<Condition>(id, (*it)));
            id++;
        }
        if(!MasterToSlave)
        {
            for (auto& node : coupling_conditions.Nodes()) {
                node.Fix(DISPLACEMENT_X);
                node.Fix(DISPLACEMENT_Y);
                node.Fix(DISPLACEMENT_Z);
            }
        }

        KRATOS_WATCH(*(coupling_conditions.NodesBegin()))
    }
}
