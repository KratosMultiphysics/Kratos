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
#include "iga_modeler.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void IgaModeler::GenerateModelPart(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Parameters PrepareGeometryParameters,
        IndexType EchoLevel) const
    {
        KRATOS_ERROR_IF_NOT(PrepareGeometryParameters.Has("element_condition_list"))
            << "Missing \"element_condition_list\" section" << std::endl;

        KRATOS_ERROR_IF_NOT(PrepareGeometryParameters["element_condition_list"].IsArray())
            << "\"element_condition_list\" needs to be an array." << std::endl;

        for (SizeType i = 0; i < PrepareGeometryParameters["element_condition_list"].size(); ++i)
        {
            CreateIntegrationDomain(
                rOriginModelPart,
                rDestinationModelPart,
                PrepareGeometryParameters["element_condition_list"][i],
                EchoLevel);
        }
    }

    ///@}

    void IgaModeler::CreateIntegrationDomain(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rParameters,
        IndexType EchoLevel) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("iga_model_part"))
            << "\"iga_model_part\" need to be specified." << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters.Has("parameters"))
            << "\"parameters\" need to be specified." << std::endl;

        std::string sub_model_part_name = rParameters["iga_model_part"].GetString();

        ModelPart& sub_model_part = rModelPart.HasSubModelPart(sub_model_part_name)
            ? rModelPart.GetSubModelPart(sub_model_part_name)
            : rModelPart.CreateSubModelPart(sub_model_part_name);

        // Generate the list of geometries, which are needed, here.
        GeometriesArrayType geometry_list;
        GetCadGeometryList(geometry_list, rModelPart, rParameters, EchoLevel);

        if (!rParameters.Has("geometry_type")) {
            CreateQuadraturePointGeometries(
                geometry_list, sub_model_part, rParameters["parameters"], EchoLevel);
        }
        else {
            std::string geometry_type = rParameters["geometry_type"].GetString();
            if (geometry_type == "GeometrySurfaceNodes"
                || geometry_type == "GeometrySurfaceVariationNodes"
                || geometry_type == "GeometryCurveNodes"
                || geometry_type == "GeometryCurveVariationNodes") {}
            else {
                CreateQuadraturePointGeometries(
                    geometry_list, sub_model_part, rParameters["parameters"], EchoLevel);
            }
        }
        KRATOS_INFO_IF("CreateIntegrationDomainElementCondition", EchoLevel > 3)
            << "Creation of elements/ conditions finished in: " << sub_model_part << std::endl;
    }

    void IgaModeler::CreateQuadraturePointGeometries(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        const Parameters& rParameters,
        IndexType EchoLevel) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("type"))
            << "\"type\" need to be specified." << std::endl;
        std::string type = rParameters["type"].GetString();
        KRATOS_ERROR_IF_NOT(rParameters.Has("name"))
            << "\"name\" need to be specified." << std::endl;
        std::string name = rParameters["name"].GetString();

        SizeType shape_function_derivatives_order = 1;
        if (rParameters.Has("shape_function_derivatives_order")) {
            shape_function_derivatives_order = rParameters["shape_function_derivatives_order"].GetInt();
        }
        else {
            KRATOS_INFO_IF("CreateQuadraturePointGeometries", EchoLevel > 4)
                << "shape_function_derivatives_order is not provided and thus being considered as 1. " << std::endl;
        }

        KRATOS_INFO_IF("CreateQuadraturePointGeometries", EchoLevel > 0)
            << "Creating " << name << "s of type: " << type
            << " for " << rQuadraturePointGeometryList.size() << " geometries"
            << " in " << rCadSubModelPart.Name() << "-SubModelPart." << std::endl;

        for (SizeType i = 0; i < rQuadraturePointGeometryList.size(); ++i)
        {
            GeometriesArrayType geometries;
            //rQuadraturePointGeometryList[i].CreateQuadraturePointGeometries(
            //    geometries, shape_function_derivatives_order);

            KRATOS_INFO_IF("CreateQuadraturePointGeometries", EchoLevel > 1)
                << geometries.size() << " quadrature point geometries have been created." << std::endl;

            if (type == "element" || type == "Element") {
                SizeType id = 0;
                if (rCadSubModelPart.GetRootModelPart().Elements().size() > 0)
                    id = rCadSubModelPart.GetRootModelPart().Elements().back().Id() + 1;

                Modeler::CreateElements<GeometriesArrayType>(
                    geometries.begin(), geometries.end(),
                    rCadSubModelPart, name, id, nullptr, EchoLevel);
            }
            else if (type == "condition" || type == "Condition") {
                SizeType id = 0;
                if (rCadSubModelPart.GetRootModelPart().Conditions().size() > 0)
                    id = rCadSubModelPart.GetRootModelPart().Conditions().back().Id() + 1;

                Modeler::CreateConditions<GeometriesArrayType>(
                    geometries.begin(), geometries.end(),
                    rCadSubModelPart, name, id, nullptr, EchoLevel);
            }
            else {
                KRATOS_ERROR << "\"type\" does not exist: " << type
                    << ". Possible types are \"element\" and \"condition\"." << std::endl;
            }
        }
    }

    void IgaModeler::GetCadGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters& rParameters,
        IndexType EchoLevel) const
    {
        if (rParameters.Has("brep_id")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_id"].GetInt()));
        }
        if (rParameters.Has("brep_ids")) {
            for (SizeType i = 0; i < rParameters["brep_ids"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_ids"][i].GetInt()));
            }
        }
        if (rParameters.Has("brep_name")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_name"].GetString()));
        }
        if (rParameters.Has("brep_names")) {
            for (SizeType i = 0; i < rParameters["brep_names"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_names"][i].GetString()));
            }
        }

        KRATOS_ERROR_IF(rGeometryList.size() == 0)
            << "Empty geometry list. Either \"brep_id\", \"brep_ids\", \"brep_name\" or \"brep_names\" are the possible options." << std::endl;
    }
}