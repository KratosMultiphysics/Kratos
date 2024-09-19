//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//                   Andrea Gorgi
//

// Project includes
#include "iga_iga_surface_mapping_geometries_modeler.h"
#include "custom_utilities/iga_mapping_intersection_utilities.h"

#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void IgaIgaSurfaceMappingGeometriesModeler::SetupGeometryModel()
    {
        CheckParameters();

        // Create the coupling model part inside the origin model
        ModelPart& coupling_model_part = (mpModels[0]->HasModelPart("coupling"))
            ? mpModels[0]->GetModelPart("coupling")
            : mpModels[0]->CreateModelPart("coupling");

        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;
        bool origin_is_iga = false;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
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

        // Create a coupling geometry conecting both brep surfaces 
        IgaIgaSurfaceMappingGeometriesModeler::CreateIgaIgaSurfaceCouplingGeometry(
            coupling_interface_origin,
            coupling_interface_destination,
            coupling_model_part
        );

        // Create quadrature point geometries in the origin and destination domain
        IgaMappingIntersectionUtilities::CreateIgaIgaQuadraturePointsCoupling2DGeometries3D(
            coupling_model_part);
    }

    void IgaIgaSurfaceMappingGeometriesModeler::CheckParameters()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_model_part_name"))
            << "Missing \"origin_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_model_part_name"))
            << "Missing \"destination_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;

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

     void IgaIgaSurfaceMappingGeometriesModeler::CreateIgaIgaSurfaceCouplingGeometry(ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult)
    {
        std::vector<IndexType> patchs_id_DomainA, patchs_id_DomainB;

        // Create a vector with the ids of the different brep surfaces for domain A
        for (auto condition_itr = rModelPartDomainA.ConditionsBegin(); condition_itr != rModelPartDomainA.ConditionsEnd(); condition_itr++){
            if (condition_itr == rModelPartDomainA.ConditionsBegin()){
                patchs_id_DomainA.push_back(condition_itr->pGetGeometry()->GetGeometryParent(0).Id());
                continue;
            }

            std::ostringstream condition_name_stream;
            condition_itr->pGetGeometry()->PrintInfo(condition_name_stream);
            std::string condition_name = condition_name_stream.str();

            // We should avoid calling the GetGeometryParent(0) method for the Point load conditions!!
            if (condition_name == "a point in 3D space"){
                continue;
            }

            IndexType current_patch_id = condition_itr->pGetGeometry()->GetGeometryParent(0).Id();
            
            
            auto find_patch_id_itr = std::find(patchs_id_DomainA.begin(), patchs_id_DomainA.end(), current_patch_id);

            if (find_patch_id_itr != patchs_id_DomainA.end()) {
                continue; // The patch id is already inside the vector
            } else {
                patchs_id_DomainA.push_back(current_patch_id);
            }
        }

        // Create a vector with the ids of the different brep surfaces for domain B
        for (auto condition_itr = rModelPartDomainB.ConditionsBegin(); condition_itr != rModelPartDomainB.ConditionsEnd(); condition_itr++){
            if (condition_itr == rModelPartDomainB.ConditionsBegin()){
                patchs_id_DomainB.push_back(condition_itr->pGetGeometry()->GetGeometryParent(0).Id());
                continue;
            }

            std::ostringstream condition_name_stream;
            condition_itr->pGetGeometry()->PrintInfo(condition_name_stream);
            std::string condition_name = condition_name_stream.str();

            // We should avoid calling the GetGeometryParent(0) method for the Point load conditions!!
            if (condition_name == "a point in 3D space"){
                continue;
            }

            IndexType current_patch_id = condition_itr->pGetGeometry()->GetGeometryParent(0).Id();
            
            
            auto find_patch_id_itr = std::find(patchs_id_DomainB.begin(), patchs_id_DomainB.end(), current_patch_id);

            if (find_patch_id_itr != patchs_id_DomainB.end()) {
                continue; // The patch id is already inside the vector
            } else {
                patchs_id_DomainB.push_back(current_patch_id);
            }
        }


        for (IndexType i = 0; i <  patchs_id_DomainA.size(); i++){
            for (IndexType j = 0; j <  patchs_id_DomainB.size(); j++){
                // Get a pointer to the brep surface geometry
                auto p_brep_surface_A = rModelPartDomainA.GetRootModelPart().pGetGeometry(patchs_id_DomainA[i]);
                auto p_brep_surface_B = rModelPartDomainB.GetRootModelPart().pGetGeometry(patchs_id_DomainB[j]);

                /*// Cast the  brep surface pointer to the derived class 
                auto brep_surface_cast_A = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_brep_surface_A);
                auto brep_surface_cast_B = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_brep_surface_B);*/

                // Creating one coupling geometry for each combination of brep surface 
                rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        p_brep_surface_A, p_brep_surface_B));
                
            }
        }
    }

}
