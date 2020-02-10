//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Pooyan Dadvand
//
//


// Project includes
#include "cad_integration_domain.h"


namespace Kratos
{

    ///@}
    ///@name Life Cycle
    ///@{

    static void CadIntegrationDomain::CreateIntegrationDomain(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rPhysicsParameters,
        int EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("element_condition_list"))
            << "Missing \"element_condition_list\" section" << std::endl;

        Kratos::CadIntegrationDomain::CreateIntegrationDomainElementConditionList(
            rModelPart,
            rCadModelPart,
            rPhysicsParameters["element_condition_list"],
            EchoLevel);
    }

    static void CadIntegrationDomain::CreateIntegrationDomainElementConditionList(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rElementConditionListParameters,
        int EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rElementConditionListParameters.IsArray())
            << "\"element_condition_list\" needs to be an array." << std::endl;

        for (SizeType i = 0; i < rElementConditionListParameters.size(); ++i)
        {
            CadIntegrationDomain::CreateIntegrationDomainElementCondition(
                rModelPart,
                rCadModelPart,
                rElementConditionListParameters[i],
                EchoLevel);
        }
    }

    static void CadIntegrationDomain::CreateIntegrationDomainElementCondition(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("geometry_type"))
            << "\"geometry_type\" need to be specified." << std::endl;
        std::string geometry_type = rParameters["geometry_type"].GetString();

        KRATOS_ERROR_IF_NOT(rParameters.Has("iga_model_part"))
            << "\"iga_model_part\" need to be specified." << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters.Has("parameters"))
            << "\"parameters\" need to be specified." << std::endl;

        std::string sub_model_part_name = rParameters["iga_model_part"].GetString();

        ModelPart& sub_model_part = rModelPart.HasSubModelPart(sub_model_part_name)
            ? rModelPart.GetSubModelPart(sub_model_part_name)
            : rModelPart.CreateSubModelPart(sub_model_part_name);

        KRATOS_INFO_IF("CreateIntegrationDomainElementCondition", EchoLevel > 1)
            << "SubModelPart: " << sub_model_part_name << " prepared for geometrical integration of: "
            << geometry_type << "s" << std::endl;

        // Generate the list of geometries, which are needed, here.
        GeometriesArrayType geometry_list;
        CadIntegrationDomain::GetGeometryList(geometry_list, rModelPart, rParameters, EchoLevel);

        if (geometry_type == "GeometrySurfaceNodes"
            || geometry_type == "GeometrySurfaceVariationNodes"
            || geometry_type == "GeometryCurveNodes"
            || geometry_type == "GeometryCurveVariationNodes") {
            CadIntegrationDomain::GetGeometrySurfacePointsAt(
                geometry_list, sub_model_part, geometry_type, rParameters["parameters"], 0, EchoLevel);
        }
        else {
            CadIntegrationDomain::CreateQuadraturePointGeometries(
                geometry_list, sub_model_part, rParameters["parameters"], EchoLevel);
        }
        KRATOS_INFO_IF("CreateIntegrationDomainElementCondition", EchoLevel > 3)
            << "Creation of elements/ conditions finished in: " << sub_model_part << std::endl;
    }

    static void CadIntegrationDomain::CreateQuadraturePointGeometries(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0)
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
            rQuadraturePointGeometryList[i].CreateQuadraturePointGeometries(
                geometries, shape_function_derivatives_order);

            KRATOS_INFO_IF("CreateQuadraturePointGeometries", EchoLevel > 1)
                << geometries.size() << " quadrature point geometries have been created." << std::endl;

            if (type == "element" || type == "Element") {
                int id = 0;
                if (rCadSubModelPart.GetRootModelPart().Elements().size() > 0)
                    id = rCadSubModelPart.GetRootModelPart().Elements().back().Id() + 1;

                CadIntegrationDomain::CreateElements(geometries, rCadSubModelPart, name, id, EchoLevel);
            }
            else if (type == "condition" || type == "Condition") {
                int id = 0;
                if (rCadSubModelPart.GetRootModelPart().Conditions().size() > 0)
                    id = rCadSubModelPart.GetRootModelPart().Conditions().back().Id() + 1;

                CadIntegrationDomain::CreateConditions(geometries, rCadSubModelPart, name, id, EchoLevel);
            }
            else {
                KRATOS_ERROR << "\"type\" does not exist: " << type
                    << ". Possible types are \"element\" and \"condition\"." << std::endl;
            }
        }
    }

    static void CadIntegrationDomain::CreateElements(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        std::string& rElementName,
        int& rIdCounter,
        int EchoLevel = 0)
    {
        const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

        ElementsContainerType new_element_list;

        KRATOS_INFO_IF("CreateElements", EchoLevel > 2)
            << "Creating " << rQuadraturePointGeometryList.size()
            << " elements of type " << rElementName
            << " in " << rCadSubModelPart.Name() << "-SubModelPart." << std::endl;

        for (IndexType i = 0; i < rQuadraturePointGeometryList.size(); ++i)
        {
            auto p_element = rReferenceElement.Create(rIdCounter, rQuadraturePointGeometryList(i), nullptr);
            rIdCounter++;
            new_element_list.push_back(p_element);
        }

        rCadSubModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }

    static void CadIntegrationDomain::CreateConditions(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        std::string& rConditionName,
        int& rIdCounter,
        int EchoLevel = 0)
    {
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;
        new_condition_list.reserve(rQuadraturePointGeometryList.size());

        KRATOS_INFO_IF("CreateConditions", EchoLevel > 2)
            << "Creating " << rQuadraturePointGeometryList.size()
            << " conditions of type " << rConditionName
            << " in " << rCadSubModelPart.Name() << "-SubModelPart." << std::endl;

        for (IndexType i = 0; i < rQuadraturePointGeometryList.size(); ++i)
        {
            auto p_condition = rReferenceCondition.Create(rIdCounter, rQuadraturePointGeometryList[i], nullptr);
            rIdCounter++;
            new_condition_list.push_back(p_condition);
        }

        rCadSubModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
    }

    static void CadIntegrationDomain::GetGeometryPointsAt(
        GeometriesArrayType& rGeometryList,
        ModelPart& rCadSubModelPart,
        const std::string& rGeometryType,
        const Parameters& rParameters,
        IndexType SpecificationType,
        int EchoLevel = 0)
    {
        KRATOS_INFO_IF("GetGeometryPointsAt", EchoLevel > 1)
            << "In " << geometry_type  <<" collecting nodes." << std::endl;

        if (geometry_type == "GeometrySurfaceNodes") {
            CadIntegrationDomain::GetGeometrySurfacePointsAt(
                geometry_list, sub_model_part, rParameters["parameters"], 0, EchoLevel);
        }
        else if (geometry_type == "GeometrySurfaceVariationNodes") {
            CadIntegrationDomain::GetGeometrySurfacePointsAt(
                geometry_list, sub_model_part, rParameters["parameters"], 1, EchoLevel);
        }
        if (geometry_type == "GeometryCurveNodes") {
            CadIntegrationDomain::GetGeometryCurvePointsAt(
                geometry_list, sub_model_part, rParameters["parameters"], 0, EchoLevel);
        }
        else if (geometry_type == "GeometryCurveVariationNodes") {
            CadIntegrationDomain::GetGeometryCurvePointsAt(
                geometry_list, sub_model_part, rParameters["parameters"], 1, EchoLevel);
        }
    }

    static void CadIntegrationDomain::GetGeometrySurfacePointsAt(
        GeometriesArrayType& rGeometryList,
        ModelPart& rCadSubModelPart,
        const Parameters& rParameters,
        IndexType SpecificationType,
        int EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("local_parameters"))
            << "\"local_parameters\" need to be specified." << std::endl;
        KRATOS_ERROR_IF(rParameters["local_parameters"].size() > 3)
            << "\"local_parameters\" exceeds size. Maximum size is 3. Actual size is: "
            << rParameters["local_parameters"].size() << std::endl;

        CoordinatesArrayType local_parameters = ZeroVector(3);

        for (SizeType i = 0; i < rParameters["local_parameters"].size(); ++i) {
            local_parameters[i] = rParameters["local_parameters"][i].GetDouble();
        }

        KRATOS_INFO_IF("GetGeometrySurfacePointsAt", EchoLevel > 2)
            << "At parameters: " << local_parameters << "." << std::endl;

        const double tolerance;

        PointsArrayType points;
        for (SizeType i = 0; i < rGeometryList.size(); ++i)
        {
            // Edges
            if (std::abs(local_parameters[0] + 1) < tolerance && std::abs(local_parameters[1]) < tolerance) {
                GetPointsAtEdge(points, 0, SpecificationType);
            }
            else if (std::abs(local_parameters[0] - 1) < tolerance && std::abs(local_parameters[1] + 1) < tolerance) {
                GetPointsAtEdge(points, 1, SpecificationType);
            }
            else if (std::abs(local_parameters[0] + 1) < tolerance && std::abs(local_parameters[1] - 1) < tolerance) {
                GetPointsAtEdge(points, 2, SpecificationType);
            }
            else if (std::abs(local_parameters[0]) < tolerance && std::abs(local_parameters[1] + 1) < tolerance) {
                GetPointsAtEdge(points, 3, SpecificationType);
            }
            // Vertices
            else if (std::abs(local_parameters[0]) < tolerance && std::abs(local_parameters[1]) < tolerance) {
                GetPointsAtVertex(points, 0, SpecificationType);
            }
            else if (std::abs(local_parameters[0] - 1) < tolerance && std::abs(local_parameters[1]) < tolerance) {
                GetPointsAtVertex(points, 1, SpecificationType);
            }
            else if (std::abs(local_parameters[0] - 1) < tolerance && std::abs(local_parameters[1] - 1) < tolerance) {
                GetPointsAtVertex(points, 2, SpecificationType);
            }
            else if (std::abs(local_parameters[0]) < tolerance && std::abs(local_parameters[1] - 1) < tolerance) {
                GetPointsAtVertex(points, 3, SpecificationType);
            }
        }

        rCadSubModelPart.AddNodes(points.begin(), points.end());
    }

    static void CadIntegrationDomain::GetGeometryEdgePointsAt(
        GeometriesArrayType& rGeometryList,
        ModelPart& rCadSubModelPart,
        const Parameters& rParameters,
        IndexType SpecificationType,
        int EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("local_parameters"))
            << "\"local_parameters\" need to be specified." << std::endl;
        KRATOS_ERROR_IF(rParameters["local_parameters"].size() > 3)
            << "\"local_parameters\" exceeds size. Maximum size is 3. Actual size is: "
            << rParameters["local_parameters"].size() << std::endl;

        CoordinatesArrayType local_parameters = ZeroVector(3);

        for (SizeType i = 0; i < rParameters["local_parameters"].size(); ++i) {
            local_parameters[i] = rParameters["local_parameters"][i].GetDouble();
        }

        KRATOS_INFO_IF("GetGeometryEdgePointsAt", EchoLevel > 2)
            << "At parameters: " << local_parameters << "." << std::endl;

        const double tolerance;

        PointsArrayType points;
        for (SizeType i = 0; i < rGeometryList.size(); ++i) {
            // Vertices
            if (std::abs(local_parameters[0]) < tolerance) {
                GetPointsAtVertex(points, 0, SpecificationType);
            }
            else if (std::abs(local_parameters[0] - 1) < tolerance) {
                GetPointsAtVertex(points, 1, SpecificationType);
            }
        }

        rCadSubModelPart.AddNodes(points.begin(), points.end());
    }

    static void CadIntegrationDomain::GetGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0)
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

    ///@}

}; // namespace CadIntegrationDomain

}  // namespace Kratos.

#endif // KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED  defined