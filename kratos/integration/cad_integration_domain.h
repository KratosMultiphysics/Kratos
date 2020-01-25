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


#if !defined(KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED )
#define  KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// CadIntegrationDomain for the the CAD entities.
/* Provides functionalities and processes to
* create a numerical model from a CAD model.
*/
class CadIntegrationDomain
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::PointsArrayType PointsArrayType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    static void CreateIntegrationDomain(
        ModelPart& rModelPart,
        ModelPart &rCadModelPart,
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
private:
    static void CreateIntegrationDomainElementConditionList(
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

    static void CreateIntegrationDomainElementCondition(
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

        // Generate the list of geometries, which are needed, here.
        GeometriesArrayType geometry_list;
        CadIntegrationDomain::GetGeometryList(geometry_list, rModelPart, rParameters, EchoLevel);

        if (geometry_type == "GeometryCurveNodes" || geometry_type == "GeometrySurfaceNodes")
        {
            CadIntegrationDomain::GetGeometryPointsAt(geometry_list, sub_model_part, rParameters["parameters"], 0, EchoLevel);
        }
        if (geometry_type == "GeometryCurveVariationNodes" || geometry_type == "GeometrySurfaceVariationNodes")
        {
            CadIntegrationDomain::GetGeometryPointsAt(geometry_list, sub_model_part, rParameters["parameters"], 1, EchoLevel);
        }
        else
        {
            CadIntegrationDomain::CreateQuadraturePointGeometries(geometry_list, sub_model_part, rParameters, EchoLevel);
        }
    }

    static void CreateQuadraturePointGeometries(
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

        for (SizeType i = 0; i < rQuadraturePointGeometryList.size(); ++i)
        {
            GeometriesArrayType geometries;
            rQuadraturePointGeometryList[i].CreateQuadraturePointGeometries(
                geometries, 2);

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

    static void CreateElements(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        std::string& rElementName,
        int& rIdCounter,
        int EchoLevel = 0)
    {
        const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

        ElementsContainerType new_element_list;
        new_element_list.reserve(rQuadraturePointGeometryList.size());

        for (IndexType i = 0; i < rQuadraturePointGeometryList.size(); ++i)
        {
            auto p_element = rReferenceElement.Create(rIdCounter, rQuadraturePointGeometryList[i], nullptr);

            rIdCounter++;

            // Add nodes of new Condition to this ModelPart
            // Is that even necessary. This is a costlly operation
            rCadSubModelPart.AddNodes(p_element->pGetGeometry()->begin(), p_element->pGetGeometry()->end());

            new_element_list.push_back(p_element);
        }

        rCadSubModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }

    static void CreateConditions(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        std::string& rConditionName,
        int& rIdCounter,
        int EchoLevel = 0)
    {
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;
        new_condition_list.reserve(rQuadraturePointGeometryList.size());

        for (IndexType i = 0; i < rQuadraturePointGeometryList.size(); ++i)
        {
            auto p_condition = rReferenceCondition.Create(rIdCounter, rQuadraturePointGeometryList[i], nullptr);

            rIdCounter++;

            // Add nodes of new Condition to this ModelPart
            // Is that even necessary. This is a costlly operation
            rCadSubModelPart.AddNodes(p_condition->pGetGeometry()->begin(), p_condition->pGetGeometry()->end());

            new_condition_list.push_back(p_condition);
        }

        rCadSubModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
    }

    static void GetGeometryPointsAt(
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

        for (SizeType i = 0; i < rParameters["local_parameters"].size(); ++i)
        {
            local_parameters[i] = rParameters["local_parameters"][i].GetDouble();
        }

        PointsArrayType points;
        for (SizeType i = 0; i < rGeometryList.size(); ++i)
        {
            // Uncommented as the interface is not decided yet.
            //rGeometryList[i].GetPointsAt(
            //    points,
            //    local_parameters,
            //    SpecificationType);
        }

        rCadSubModelPart.AddNodes(points.begin(), points.end());
    }

    static void GetGeometryList(
        GeometriesArrayType rGeometryList,
        ModelPart& rModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0)
    {
        if (rParameters.Has("brep_id"))
        {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_id"].GetInt()));
        }
        if (rParameters.Has("brep_ids"))
        {
            for (SizeType i = 0; i < rParameters["brep_ids"].size(); ++i)
            {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_ids"][i].GetInt()));
            }
        }
        if (rParameters.Has("brep_name"))
        {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_name"].GetString()));
        }
        if (rParameters.Has("brep_names"))
        {
            for (SizeType i = 0; i < rParameters["brep_names"].size(); ++i)
            {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_names"][i].GetString()));
            }
        }

        KRATOS_ERROR_IF(rGeometryList.size() == 0)
            << "Empty geometry list. Either \"brep_id\", \"brep_ids\", \"brep_name\" or \"brep_names\" need to be specified." << std::endl;
    }

    ///@}

}; // namespace CadIntegrationDomain

}  // namespace Kratos.

#endif // KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED  defined