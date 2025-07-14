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

    // This implementation is to be changed after geometries allow the creation of GeometryType::Pointer from raw pointers (intrusive_ptr) -> no dynamic cast needed anymore
    // In the future: const GeometryPointerType copy_brep_curve_on_surface_geometry = Kratos::make_intrusive<BrepCurveOnSurfacePointer::element_type>(&brep_curve_on_surface_geometry);
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

            // If the brep curve on surface is not already in the vector, add it
            if (new_brep_curve_on_surface_id != old_brep_curve_on_surface_id)
            {
                // Get the nurbs curve on surface geometry from the brep curve on surface geometry
                const GeometryPointerType nurbs_curve_on_surface_geometry = brep_curve_on_surface_geometry.pGetGeometryPart(std::numeric_limits<IndexType>::max() - 2);

                // Downcast it from the base class (geometry) to the derived class (NurbsCurveOnSurfaceGeometry) 
                IgaFEMMappingGeometriesModeler::NurbsCurveOnSurfacePointer interface_nurbs_curve_on_surface_cast = std::dynamic_pointer_cast<IgaFEMMappingGeometriesModeler::NurbsCurveOnSurfacePointer::element_type>(nurbs_curve_on_surface_geometry); 

                // Get a pointer to the underlying surface (get a pointer to the base class)
                const GeometryPointerType interface_nurbs_surface = brep_curve_on_surface_geometry.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

                // Downcast the pointer to the derived class (NurbsSurfaceGeometry)
                IgaFEMMappingGeometriesModeler::NurbsSurfacePointer interface_nurbs_surface_cast = std::dynamic_pointer_cast<IgaFEMMappingGeometriesModeler::NurbsSurfacePointer::element_type>(interface_nurbs_surface); 
                KRATOS_ERROR_IF_NOT(interface_nurbs_surface_cast)
                    << "failed to downcast interface_nurbs_surface";

                // Get the Nurbs curve in the parameter space 
                auto p_nurbs_curve = interface_nurbs_curve_on_surface_cast->pGetCurve();

                // intersection between curve and surface parameter space in terms of the curve 1D parameter space.
                std::vector<double> curve_and_parameter_space_intersections;
                nurbs_curve_on_surface_geometry->SpansLocalSpace(curve_and_parameter_space_intersections);

                NurbsInterval active_range(curve_and_parameter_space_intersections.front(),curve_and_parameter_space_intersections.back());

                GeometryType::Pointer p_reconstructed_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfacePointer::element_type>(interface_nurbs_surface_cast, p_nurbs_curve, active_range);
                KRATOS_WATCH(p_reconstructed_brep_curve_on_surface->Center())
                p_brep_curve_on_surface_vector.push_back(p_reconstructed_brep_curve_on_surface);

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
