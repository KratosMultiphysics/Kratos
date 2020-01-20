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
namespace CadIntegrationDomain
{

public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef Geometry<Node<3>> GeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::PointsArrayType PointsArrayType;

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

        CreateIntegrationDomainElementConditionList(
            rModelPart,
            rCadModelPart,
            rPhysicsParameters["element_condition_list"],
            EchoLevel);
    }

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
            CreateIntegrationDomainElementCondition(
                rModelPart,
                rCadModelPart,
                rElementConditionListParameters[i],
                EchoLevel);
        }
    }

private:

    static void CreateIntegrationDomainElementCondition(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0)
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("geometry_type"))
            << "\"geometry_type\" need to be specified." << std::endl;
        std::string geometry_type = rParameters["geometry_type"];

        KRATOS_ERROR_IF_NOT(rParameters.Has("iga_model_part"))
            << "\"iga_model_part\" need to be specified." << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters.Has("parameters"))
            << "\"parameters\" need to be specified." << std::endl;

        std::string sub_model_part_name = rParameters["iga_model_part"].GetString();

        ModelPart& iga_model_part = rModelPart.HasSubModelPart(sub_model_part_name)
            ? rModelPart.GetSubModelPart(iga_model_part)
            : rModelPart.CreateSubModelPart(iga_model_part);

        std::vector<GeometryType> geometry_list;
        GetGeometryList(geometry_list, rModelPart, rParameters, EchoLevel);

        if (geometry_type == "GeometryCurveNodes" || geometry_type == "GeometrySurfaceNodes")
        {
            GetGeometryPointsAt(geometry_list, rModelPart, rCadModelPart, rParameters["parameters"], 0, EchoLevel);
        }
        if (geometry_type == "GeometryCurveVariationNodes" || geometry_type == "GeometrySurfaceVariationNodes")
        {
            GetGeometryPointsAt(geometry_list, rModelPart, rCadModelPart, rParameters["parameters"], 1, EchoLevel);
        }
    }

    static void GetGeometryPointsAt(
        std::vector<GeometryType>& rGeometryList,
        ModelPart& rModelPart,
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
            rGeometryList[i]->GetPointsAt(
                points,
                local_parameters,
                SpecificationType);
        }

        rCadSubModelPart.AddNodes(rGeometryList.Begin(), rGeometryList.End());
    }

    static void GetGeometryList(
        std::vector<GeometryType> rGeometryList,
        ModelPart& rModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0)
    {
        if (rParameters.Has("brep_id"))
        {
            rGeometryList.push_back(rModelPart.pGetGeometryPart(rParameters["brep_id"].GetInt()));
        }
        if (rParameters.Has("brep_ids"))
        {
            for (SizeType i = 0; i < rParameters["brep_ids"].size(); ++i)
            {
                rGeometryList.push_back(rModelPart.pGetGeometryPart(rParameters["brep_ids"][i].GetInt()));
            }
        }
        if (rParameters.Has("brep_name"))
        {
            rGeometryList.push_back(rModelPart.pGetGeometryPart(rParameters["brep_name"].GetString()))
        }
        if (rParameters.Has("brep_names"))
        {
            for (SizeType i = 0; i < rParameters["brep_names"].size(); ++i)
            {
                rGeometryList.push_back(rModelPart.pGetGeometryPart(rParameters["brep_names"][i].GetString()));
            }
        }

        KRATOS_ERROR_IF(rGeometryList.size() == 0)
            << "Empty geometry list. Either \"brep_id\", \"brep_ids\", \"brep_name\" or \"brep_names\" need to be specified." << std::endl;
    }

    ///@}

}; // namespace CadIntegrationDomain

}  // namespace Kratos.

#endif // KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED  defined