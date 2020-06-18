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
#include "mapping_geometries_modeler.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void MappingGeometriesModeler::SetupGeometryModel()
    {
        CheckParameters();

        const std::string origin_model_part_name = mParameters["origin_model_part_name"].GetString();
        ModelPart& origin_model_part = mpModel->HasModelPart(origin_model_part_name)
            ? mpModel->GetModelPart(origin_model_part_name)
            : mpModel->CreateModelPart(origin_model_part_name);

        const std::string destination_model_part_name = mParameters["destination_model_part_name"].GetString();
        ModelPart& destination_model_part = mpModel->HasModelPart(destination_model_part_name)
            ? mpModel->GetModelPart(destination_model_part_name)
            : mpModel->CreateModelPart(destination_model_part_name);

        ModelPart& coupling_model_part = mpModel->HasModelPart("coupling")
            ? mpModel->GetModelPart("coupling")
            : mpModel->CreateModelPart("coupling");

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
            // Some future functionality to automatically determine interfaces?
            // (Create interface sub model parts in origin and destination modelparts)
            // (set strings to correct values)
            //origin_interface_sub_model_part_name = ???
            //destination_interface_sub_model_part_name = ???
        }

        // Transfer everything into the coupling modelpart
        ModelPart& coupling_interface_origin = coupling_model_part.HasSubModelPart("interface_origin")
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");
        CopySubModelPart(coupling_interface_origin,
            origin_model_part.GetSubModelPart(origin_interface_sub_model_part_name));

        ModelPart& coupling_interface_destination = coupling_model_part.HasSubModelPart("interface_destination")
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");
        CopySubModelPart(coupling_interface_destination,
            destination_model_part.GetSubModelPart(destination_interface_sub_model_part_name));


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

    void MappingGeometriesModeler::CheckParameters()
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
    ///@}
}
