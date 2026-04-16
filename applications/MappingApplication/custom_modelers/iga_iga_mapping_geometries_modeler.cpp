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
#include "iga_iga_mapping_geometries_modeler.h"
#include "custom_utilities/iga_mapping_intersection_utilities.h"
#include <ratio>

namespace Kratos
{
    ///@name Stages
    ///@{

    void IgaIgaMappingGeometriesModeler::SetupGeometryModel()
    {
        CheckParameters();
        
        // Create the coupling model part
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
        const bool is_surface_mapping = mParameters["is_surface_mapping"].GetBool();

        if (is_surface_mapping){
            KRATOS_ERROR_IF_NOT(mParameters.Has("search_radius")) << "'search_radius' was not specified in the CoSim parameters file\n";
            mSearchRadius = mParameters["search_radius"].GetDouble();
        }

        ModelPart& r_origin_interface_mp = mpModels[0]->GetModelPart(origin_interface_sub_model_part_name);
        ModelPart& r_destination_interface_mp = mpModels[1]->GetModelPart(destination_interface_sub_model_part_name);

        if (mParameters.Has("origin_reconstruction_settings") &&
            mParameters["origin_reconstruction_settings"]["reconstruct_brep_geometry"].GetBool())
        {
            ReconstructBrepGeometryFromFile(
                r_origin_interface_mp,
                mParameters["origin_reconstruction_settings"]);
        }

        if (mParameters.Has("destination_reconstruction_settings") &&
            mParameters["destination_reconstruction_settings"]["reconstruct_brep_geometry"].GetBool())
        {
            ReconstructBrepGeometryFromFile(
                r_destination_interface_mp,
                mParameters["destination_reconstruction_settings"]);
        }

        // Transfer everything into the coupling modelpart
       ModelPart& coupling_interface_origin = (coupling_model_part.HasSubModelPart("interface_origin"))
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");
        CopySubModelPartBrepSurfaceInterface(coupling_interface_origin,
            mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));

        ModelPart& coupling_interface_destination = (coupling_model_part.HasSubModelPart("interface_destination"))
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");
        CopySubModelPartBrepSurfaceInterface(coupling_interface_destination,
            mpModels[1]->GetModelPart(destination_interface_sub_model_part_name));

        if (is_surface_mapping)
        {
            // Create a coupling geometry conecting both brep surfaces 
            IgaIgaMappingGeometriesModeler::CreateIgaIgaSurfaceCouplingGeometry(
                coupling_interface_origin,
                coupling_interface_destination,
                coupling_model_part
            );
             // Create quadrature point geometries in the origin and destination domain
            IgaMappingIntersectionUtilities::CreateIgaIgaQuadraturePointsCoupling2DGeometries3D(
                coupling_model_part);
        }
        else 
        {
        }
    }

    void IgaIgaMappingGeometriesModeler::CheckParameters()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_model_part_name"))
            << "Missing \"origin_model_part_name\" in IgaIgaMappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_model_part_name"))
            << "Missing \"destination_model_part_name\" in IgaIgaMappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("is_interface_sub_model_parts_specified"))
            << "Missing \"is_interface_sub_model_parts_specified\" in IgaIgaMappingGeometriesModeler Parameters." << std::endl;
        
        KRATOS_ERROR_IF_NOT(mParameters.Has("is_surface_mapping"))
            << "Missing \"is_surface_mapping\" in IgaIgaMappingGeometriesModeler Parameters." << std::endl;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            KRATOS_ERROR_IF_NOT(mParameters.Has("origin_interface_sub_model_part_name"))
                << "Missing \"origin_interface_sub_model_part_name\" in IgaIgaMappingGeometriesModeler Parameters." << std::endl;

            KRATOS_ERROR_IF_NOT(mParameters.Has("destination_interface_sub_model_part_name"))
                << "Missing \"destination_interface_sub_model_part_name\" in IgaIgaMappingGeometriesModeler Parameters." << std::endl;
        }
    }

    void IgaIgaMappingGeometriesModeler::CreateIgaIgaSurfaceCouplingGeometry(ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult)
    {
        std::vector<IndexType> patches_id_domain_A, patches_id_domain_B;

        // Looping over the geometries to get the brep surface ids 
        for (const auto& r_geom_A : rModelPartDomainA.Geometries())
        {
            if (r_geom_A.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Brep_Surface){
                patches_id_domain_A.push_back(r_geom_A.Id());
            }
        }

        for (const auto& r_geom_B : rModelPartDomainB.Geometries())
        {
            if (r_geom_B.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Brep_Surface){
                patches_id_domain_B.push_back(r_geom_B.Id());
            }
        }

        const IndexType patch_divisions = 100;
        IgaMappingIntersectionUtilities::PatchCacheMap patch_cache = IgaMappingIntersectionUtilities::BuildPatchCaches(
            patches_id_domain_A, rModelPartDomainA, patch_divisions);

        /* - Iterate over the patches in the domain B
           - For each patch in domain B, we should decide which patches in domain A are the ones that the element may have a projection on 
           - Create a new list called patchs_id_projection depending on the patch considered */
        for (auto patch_id_domain_B : patches_id_domain_B){
            // Get a pointer to the brep surface geometry of domain B
            auto p_brep_surface_domain_B = rModelPartDomainB.GetRootModelPart().pGetGeometry(patch_id_domain_B);
            // Cast the brep surface pointer  to the derived class
            auto brep_surface_cast_domain_B = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(p_brep_surface_domain_B);

            // Get the coordinates of the patch center
            CoordinatesArrayType patch_center = rModelPartDomainB.GetGeometry(patch_id_domain_B).Center();

            std::vector<IndexType> patches_with_probable_projection_id = IgaMappingIntersectionUtilities::GetPatchesWithProbableProjection(patches_id_domain_A, patch_cache, patch_center, mSearchRadius);

            for (IndexType i = 0; i < patches_with_probable_projection_id.size(); i++){
                // Get a pointer to the brep surface geometry of domain A
                auto p_brep_surface_domain_A = rModelPartDomainA.GetRootModelPart().pGetGeometry(patches_with_probable_projection_id[i]);

                // Cast the nurbs surface pointer and brep surface to the derived class
                auto brep_surface_cast_domain_A = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(p_brep_surface_domain_A);

                // Creating one coupling geometry for each combination of possible patches
                rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        p_brep_surface_domain_A, p_brep_surface_domain_B));
            }       
        }
    }

}
