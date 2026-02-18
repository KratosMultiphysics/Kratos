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

// Project includes
#include "iga_fem_mapping_geometries_modeler.h"
#include "custom_utilities/iga_mapping_intersection_utilities.h"

namespace Kratos
{
    ///@name Stages
    ///@{

    void IgaFEMMappingGeometriesModeler::SetupGeometryModel()
    {
        CheckParameters();

        ModelPart& coupling_model_part = (mpModels[0]->HasModelPart("coupling"))
            ? mpModels[0]->GetModelPart("coupling")
            : mpModels[0]->CreateModelPart("coupling");

        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
        }

        // create coupling conditions on interface depending on the formulation
        const bool is_origin_iga = mParameters["is_origin_iga"].GetBool();
        const bool is_surface_mapping = mParameters["is_surface_mapping"].GetBool();

        KRATOS_ERROR_IF(is_surface_mapping && !is_origin_iga)
            << "Surface mapping with this modeler requires the origin ModelPart to be IGA.\n"
            << "Got is_origin_iga = false." << std::endl;

        double search_radius = 0.0;
        if (is_surface_mapping){
            KRATOS_ERROR_IF_NOT(mParameters.Has("search_radius")) << "'search_radius' was not specified in the CoSim parameters file\n";
            search_radius = mParameters["search_radius"].GetDouble();
        }
       
        if (is_origin_iga && !is_surface_mapping)
        {
            IgaFEMMappingGeometriesModeler::CreateIgaInterfaceBrepCurveOnSurfaceConditions(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));
            IgaFEMMappingGeometriesModeler::CreateFEMInterfaceNurbsCurveConditions(mpModels.back()->GetModelPart(destination_interface_sub_model_part_name));
        }
        else if (!is_origin_iga && !is_surface_mapping)
        {
            IgaFEMMappingGeometriesModeler::CreateIgaInterfaceBrepCurveOnSurfaceConditions(mpModels.back()->GetModelPart(destination_interface_sub_model_part_name));
            IgaFEMMappingGeometriesModeler::CreateFEMInterfaceNurbsCurveConditions(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));
        }
        
        // Transfer everything into the coupling modelpart
        ModelPart& coupling_interface_origin = (coupling_model_part.HasSubModelPart("interface_origin"))
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");
        ModelPart& coupling_interface_destination = (coupling_model_part.HasSubModelPart("interface_destination"))
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");

        if (is_origin_iga && !is_surface_mapping) {
            CopySubModelPartIgaInterface(coupling_interface_origin,
                mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));
            CopySubModelPartFEMInterface(coupling_interface_destination,
                mpModels[1]->GetModelPart(destination_interface_sub_model_part_name));
        } else if (!is_origin_iga && !is_surface_mapping){
            CopySubModelPartFEMInterface(coupling_interface_origin,
                mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));
            CopySubModelPartIgaInterface(coupling_interface_destination,
                mpModels[1]->GetModelPart(destination_interface_sub_model_part_name));
        } else if (is_origin_iga && is_surface_mapping){
            KRATOS_ERROR_IF(mpModels[1]->GetModelPart(destination_interface_sub_model_part_name).NumberOfConditions() == 0) 
                << "the destination model part has no conditions for the mapping. Please change the elements to conditions";
            
            for (const auto& r_cond : mpModels[1]->GetModelPart(destination_interface_sub_model_part_name).Conditions())
            {
                const auto& r_geom = r_cond.GetGeometry();

                KRATOS_ERROR_IF_NOT(r_geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3 || 
                                r_geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
                    << "Surface mapping currently supports only TRIANGULAR interface conditions.\n"
                    << "Condition Id: " << r_cond.Id() << "\n"
                    << "Geometry: " << r_geom.Info() << "\n"
                    << "PointsNumber(): " << r_geom.PointsNumber() << "\n";
            }

            CopySubModelPartSurfaceInterface(coupling_interface_origin,
                mpModels[0]->GetModelPart(origin_interface_sub_model_part_name), true);
            CopySubModelPartSurfaceInterface(coupling_interface_destination,
                mpModels[1]->GetModelPart(destination_interface_sub_model_part_name), false);
        }

        KRATOS_ERROR_IF(!is_surface_mapping && (coupling_interface_origin.NumberOfConditions() == 0 || coupling_interface_destination.NumberOfConditions() == 0))
            << "Coupling geometries are currently determined by conditions in the coupling sub model parts,"
            << " but there are currently not conditions in the coupling interface origin sub model part. Please specify some."
            << std::endl;

        if (!is_surface_mapping)
        {
            IgaMappingIntersectionUtilities::CreateIgaFEMCouplingGeometriesOnCurve(
                coupling_interface_origin,
                coupling_interface_destination,
                is_origin_iga,
                coupling_model_part, 1e-6);
            IgaMappingIntersectionUtilities::CreateIgaFEMQuadraturePointsOnCurve(
                coupling_model_part, 1e-6);
        } else {
            IgaMappingIntersectionUtilities::PatchCacheMap patch_cache;

            // Create coupling geometries connecting each finite element with the IGA surface 
            IgaMappingIntersectionUtilities::CreateIgaFEMCouplingGeometriesOnSurface(
                coupling_interface_origin,
                coupling_interface_destination,
                coupling_model_part, 
                is_origin_iga,
                search_radius,
                patch_cache);
            
             // Create quadrature point geometries in the origin and destination domain
            IgaMappingIntersectionUtilities::CreateIgaFEMQuadraturePointsOnSurface(
                coupling_model_part, 
                is_origin_iga,
                patch_cache,
                search_radius);
        }
    }

    void IgaFEMMappingGeometriesModeler::CheckParameters()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_model_part_name"))
            << "Missing \"origin_model_part_name\" in IgaFEMMappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_model_part_name"))
            << "Missing \"destination_model_part_name\" in IgaFEMMappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("is_interface_sub_model_parts_specified"))
            << "Missing \"is_interface_sub_model_parts_specified\" in IgaFEMMappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("is_origin_iga"))
            << "Missing \"is_origin_iga\" in IgaFEMMappingGeometriesModeler Parameters." << std::endl;
        
        KRATOS_ERROR_IF_NOT(mParameters.Has("is_surface_mapping"))
            << "Missing \"is_surface_mapping\" in IgaFEMMappingGeometriesModeler Parameters." << std::endl;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            KRATOS_ERROR_IF_NOT(mParameters.Has("origin_interface_sub_model_part_name"))
                << "Missing \"origin_interface_sub_model_part_name\" in IgaFEMMappingGeometriesModeler Parameters." << std::endl;

            KRATOS_ERROR_IF_NOT(mParameters.Has("destination_interface_sub_model_part_name"))
                << "Missing \"destination_interface_sub_model_part_name\" in IgaFEMMappingGeometriesModeler Parameters." << std::endl;
        }
    }

    // In this method, brep_curve_on_surface conditions are created on the iga interface. This is used later to compute the intersection with the FEM elements
    void IgaFEMMappingGeometriesModeler::CreateIgaInterfaceBrepCurveOnSurfaceConditions(ModelPart& rInterfaceModelPart)
    {   
        // Create the sub model part for coupling conditions
        rInterfaceModelPart.CreateSubModelPart("coupling_conditions");
        ModelPart& coupling_conditions = rInterfaceModelPart.GetSubModelPart("coupling_conditions");
        const ModelPart& root_mp = rInterfaceModelPart.GetRootModelPart();

        // This is the vector of brep curves on surface that will create the coupling conditions
        std::vector<GeometryType::Pointer> p_brep_curve_on_surface_vector;

        // Loop over all quadrature points in the interface model part to detect the different brep curves on surface
        IndexType old_brep_curve_on_surface_id = 0;
        for (auto& r_cond : rInterfaceModelPart.Conditions()){
            // Get the quadrature point geometry
            const GeometryPointerType p_geometry = r_cond.pGetGeometry();

            // Get the parent geometry of the quadrature point (brep curve on surface)
            const GeometryType& brep_curve_on_surface_geometry = p_geometry->GetGeometryParent(0);

            // Get the id of the new brep curve on surface 
            IndexType new_brep_curve_on_surface_id = brep_curve_on_surface_geometry.Id();

            if (new_brep_curve_on_surface_id != old_brep_curve_on_surface_id){
                for (auto geometry_itr = root_mp.GeometriesBegin(); geometry_itr != root_mp.GeometriesEnd(); geometry_itr ++){
                    if (geometry_itr->Id() == new_brep_curve_on_surface_id){
                        p_brep_curve_on_surface_vector.push_back(root_mp.pGetGeometry(geometry_itr->Id()));
                    }
                }
                old_brep_curve_on_surface_id = new_brep_curve_on_surface_id;
            }
        }

        // Determine next condition number
        IndexType condition_id = (root_mp.NumberOfConditions() == 0)
            ? 1 : (root_mp.ConditionsEnd() - 1)->Id() + 1;

        for (IndexType i = 0; i < p_brep_curve_on_surface_vector.size(); i++){
            coupling_conditions.CreateNewCondition("BrepCurveOnSurfaceCondition", condition_id, p_brep_curve_on_surface_vector[i], root_mp.ConditionsBegin()->pGetProperties(), 0);
            condition_id += 1;
        }
    }

    // In this method, nurbs_curve conditions (linear) are created on the fem interface (one curve per line element). This is used later to compute the intersection with the IGA boundary
    void IgaFEMMappingGeometriesModeler::CreateFEMInterfaceNurbsCurveConditions(ModelPart& rInterfaceModelPart)
    {
        rInterfaceModelPart.CreateSubModelPart("coupling_conditions");
        ModelPart& coupling_conditions = rInterfaceModelPart.GetSubModelPart("coupling_conditions");
        const ModelPart& root_mp = rInterfaceModelPart.GetRootModelPart();

        IndexType interface_node_id;
        IndexType trial_interface_node_id;
        IndexType trial_geom_node_id;

        // Determine next condition number
        IndexType condition_id = (root_mp.NumberOfConditions() == 0)
            ? 1 : (root_mp.ConditionsEnd() - 1)->Id() + 1;

        for (std::size_t node_index = 0; node_index < rInterfaceModelPart.NumberOfNodes() - 1; ++node_index)
        {
            interface_node_id = (rInterfaceModelPart.NodesBegin() + node_index)->Id();
            std::vector< GeometryPointerType> p_geom_vec;
            for (auto& r_elem: root_mp.Elements())
            {
                auto p_geom = r_elem.pGetGeometry();
                for (std::size_t i = 0; i < p_geom->size(); i++)
                {
                    if ((*p_geom)[i].Id() == interface_node_id)
                    {
                        p_geom_vec.push_back(p_geom);
                    }
                }
            }

            if (p_geom_vec.size() == 0) KRATOS_ERROR << "Interface node not found in modelpart geom\n";

            // Loop over all geometries that have nodes on the interface
            for (std::size_t interface_geom_index = 0; interface_geom_index < p_geom_vec.size(); interface_geom_index++)
            {
                GeometryType& r_interface_geom = *(p_geom_vec[interface_geom_index]);

                // Loop over remaining interface nodes, see if any of them are nodes in the interface geom
                for (std::size_t geom_node_index = 0; geom_node_index < r_interface_geom.size(); geom_node_index++)
                {
                    trial_geom_node_id = r_interface_geom[geom_node_index].Id();

                    for (std::size_t trial_index = node_index + 1; trial_index < rInterfaceModelPart.NumberOfNodes(); ++trial_index)
                    {
                        trial_interface_node_id = (rInterfaceModelPart.NodesBegin() + trial_index)->Id();
                        if (trial_geom_node_id == trial_interface_node_id)
                        {
                            Geometry<GeometricalObject::NodeType>::PointsArrayType brep_curve_points;
                            brep_curve_points.push_back(rInterfaceModelPart.pGetNode(interface_node_id));
                            brep_curve_points.push_back(rInterfaceModelPart.pGetNode(trial_interface_node_id));

                            // Creating the brep curve knot vector
                            Vector knot_vector = ZeroVector(4);
                            knot_vector[0] = 0.0;
                            knot_vector[1] = 0.0;
                            knot_vector[2] = 1.0;
                            knot_vector[3] = 1.0;

                            // Setting the polinomial degree of the b-rep curve 
                            int p = 1;

                            // Creating a pointer to the b-rep curve 
                            GeometryType::Pointer pBrep_curve=Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<NodeType>>>(brep_curve_points, p, knot_vector);

                            // Creating a new condition in the interface 
                            coupling_conditions.CreateNewCondition("NurbsCurveCondition", condition_id, pBrep_curve, rInterfaceModelPart.pGetProperties(0));
                            condition_id += 1;
                        }
                    }
                }
            }


        }
    }
}
