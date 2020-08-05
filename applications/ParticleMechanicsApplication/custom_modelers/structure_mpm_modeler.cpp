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
#include "structure_mpm_modeler.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void StructureMpmModeler::SetupGeometryModel()
    {
        CheckParameters();

        ModelPart& coupling_model_part = (mpModelStructure->HasModelPart("coupling"))
            ? mpModelStructure->GetModelPart("coupling")
            : mpModelStructure->CreateModelPart("coupling");

        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;
        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
        }
        else
        {
            KRATOS_ERROR << "Automatic interface creation is not implemented yet" << std::endl;
            // Some future functionality to automatically determine interfaces?
            // (Create interface sub model parts in origin and destination modelparts)
            // (set strings to correct values)
            //origin_interface_sub_model_part_name = ???
            //destination_interface_sub_model_part_name = ???
        }

        ModelPart& brackground_grid_model_part = (mpModelMpm->HasModelPart("BackgroundGrid"))
            ? mpModelMpm->GetModelPart("BackgroundGrid")
            : mpModelMpm->CreateModelPart("BackgroundGrid");

        // create coupling conditions on interface depending on the dimension
        // Todo

        std::vector<GeometryPointerType> quads_structure;
        CreateStructureQuadraturePointGeometries<
            typename ModelPart::GeometriesMapType, std::vector<GeometryPointerType>>(
                coupling_model_part.Geometries(), quads_structure);
        std::vector<GeometryPointerType> quads_mpm(quads_structure.size());
        CreateMpmQuadraturePointGeometries<2, std::vector<GeometryPointerType>>(
            quads_structure, quads_mpm, brackground_grid_model_part);

        IndexType id = 1;
        for (IndexType i = 0; i < quads_structure.size(); ++i) {
            coupling_model_part.AddCondition(Kratos::make_intrusive<Condition>(
                id + i, Kratos::make_shared<CouplingGeometry<Node<3>>>(quads_structure[i], quads_mpm[i])));
        }
    }

    void StructureMpmModeler::UpdateGeometryModel()
    {
        ModelPart& origin_model_part = (mpModelStructure->HasModelPart("origin"))
            ? mpModelStructure->GetModelPart("origin")
            : mpModelStructure->CreateModelPart("origin");

        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;
        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
        }
        else
        {
            KRATOS_ERROR << "Automatic interface creation is not implemented yet" << std::endl;
            // Some future functionality to automatically determine interfaces?
            // (Create interface sub model parts in origin and destination modelparts)
            // (set strings to correct values)
            //origin_interface_sub_model_part_name = ???
            //destination_interface_sub_model_part_name = ???
        }

        ModelPart& brackground_grid_model_part = (mpModelMpm->HasModelPart("BackgroundGrid"))
            ? mpModelMpm->GetModelPart("BackgroundGrid")
            : mpModelMpm->CreateModelPart("BackgroundGrid");

        // create coupling conditions on interface depending on the dimension
        // Todo

        std::vector<GeometryPointerType> quads_structure;
        UpdateMpmQuadraturePointGeometries<3,
            typename ModelPart::ConditionsContainerType>(
                origin_model_part.Conditions(), brackground_grid_model_part);
    }

    void StructureMpmModeler::CheckParameters()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_model_part_name"))
            << "Missing \"origin_model_part_name\" in StructureMpmModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_model_part_name"))
            << "Missing \"destination_model_part_name\" in StructureMpmModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("is_interface_sub_model_parts_specified"))
            << "Missing \"is_interface_sub_model_parts_specified\" in StructureMpmModeler Parameters." << std::endl;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            KRATOS_ERROR_IF_NOT(mParameters.Has("origin_interface_sub_model_part_name"))
                << "Missing \"origin_interface_sub_model_part_name\" in StructureMpmModeler Parameters." << std::endl;

            KRATOS_ERROR_IF_NOT(mParameters.Has("destination_interface_sub_model_part_name"))
                << "Missing \"destination_interface_sub_model_part_name\" in StructureMpmModeler Parameters." << std::endl;
        }
    }
}
