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

        if (is_origin_iga && is_surface_mapping == false)
        {
            CreateIgaInterfaceBrepCurveOnSurface(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));
            KRATOS_WATCH(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name))
            exit(0);
        }
        else
        {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
        }


        // Transfer everything into the coupling modelpart
        ModelPart& coupling_interface_origin = (coupling_model_part.HasSubModelPart("interface_origin"))
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");
        CopySubModelPart(coupling_interface_origin,
            mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));

        ModelPart& coupling_interface_destination = (coupling_model_part.HasSubModelPart("interface_destination"))
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");
        CopySubModelPart(coupling_interface_destination,
            mpModels[1]->GetModelPart(destination_interface_sub_model_part_name));

        KRATOS_ERROR_IF(coupling_interface_origin.NumberOfConditions() == 0)
            << "Coupling geometries are currently determined by conditions in the coupling sub model parts,"
            << " but there are currently not conditions in the coupling interface origin sub model part. Please specify some."
            << std::endl;
        SizeType working_dim = coupling_interface_origin.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
        SizeType local_dim = coupling_interface_origin.ConditionsBegin()->GetGeometry().LocalSpaceDimension();

        if (working_dim == 2 && local_dim == 1)
        {
            MappingIntersectionUtilities::FindIntersection1DGeometries2D(
                coupling_interface_origin,
                coupling_interface_destination,
                coupling_model_part, 1e-6);
            MappingIntersectionUtilities::CreateQuadraturePointsCoupling1DGeometries2D(
                coupling_model_part, 1e-6);
        }
        else
        {
            KRATOS_ERROR << "Creation of coupling quadrature points not yet supported for requested"
                << " working space dimension = " << working_dim << " and local space dimension = "
                << local_dim << std::endl;
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
    void IgaFEMMappingGeometriesModeler::CreateIgaInterfaceBrepCurveOnSurface(ModelPart& rInterfaceModelPart)
    {   
        // Create the sub model part for coupling conditions
        rInterfaceModelPart.CreateSubModelPart("coupling_conditions");
        ModelPart& coupling_conditions = rInterfaceModelPart.GetSubModelPart("coupling_conditions");
        const ModelPart& root_mp = rInterfaceModelPart.GetRootModelPart();

        // This is the vector of brep curves on surface that will create the coupling conditions
        std::vector<GeometryType::Pointer> p_brep_curve_on_surface_vector;

        // Loop over all quadrature points in the interface model part to detect the different brep curves on surface
        IndexType old_brep_curve_on_surface_id = 0;
        for (auto cond_it = rInterfaceModelPart.ConditionsBegin(); cond_it != rInterfaceModelPart.ConditionsEnd(); ++cond_it){
            // Get the quadrature point geometry
            const GeometryPointerType p_geometry = cond_it->pGetGeometry();

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
            coupling_conditions.CreateNewCondition("BrepCurveOnSurface", condition_id, p_brep_curve_on_surface_vector[i], root_mp.ElementsBegin()->pGetProperties(), 0);
            condition_id += 1;
        }
    }
}
