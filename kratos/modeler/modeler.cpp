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
#include "modeler.h"


namespace Kratos
{
    ///@name Operations
    ///@{

    void Modeler::GenerateModelPart(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for mesh generation." << std::endl;
    }

    void Modeler::GenerateMesh(
        ModelPart& ThisModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for mesh generation." << std::endl;
    }

    void Modeler::GenerateNodes(ModelPart& ThisModelPart)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for node generation." << std::endl;
    }

    ///@}
    ///@name Generate Elements and Conditions
    ///@{

    void CadIntegrationDomain::CreateElements(
        GeometriesArrayType& rGeometries,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        int& rIdCounter,
        SizeType EchoLevel)
    {
        const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

        ElementsContainerType new_element_list;

        KRATOS_INFO_IF("CreateElements", EchoLevel > 2)
            << "Creating " << rGeometries.size()
            << " elements of type " << rElementName
            << " in " << rDestinationModelPart.Name() << "-SubModelPart." << std::endl;

        for (IndexType i = 0; i < rGeometries.size(); ++i)
        {
            auto p_element = rReferenceElement.Create(rIdCounter, rGeometries(i), nullptr);
            rIdCounter++;
            new_element_list.push_back(p_element);
        }

        rDestinationModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }

    void CadIntegrationDomain::CreateConditions(
        GeometriesArrayType& rGeometries,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        int& rIdCounter,
        SizeType EchoLevel)
    {
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;
        new_condition_list.reserve(rQuadraturePointGeometryList.size());

        KRATOS_INFO_IF("CreateConditions", rGeometries > 2)
            << "Creating " << rGeometries.size()
            << " conditions of type " << rConditionName
            << " in " << rDestinationModelPart.Name() << "-SubModelPart." << std::endl;

        for (IndexType i = 0; i < rQuadraturePointGeometryList.size(); ++i)
        {
            auto p_condition = rReferenceCondition.Create(rIdCounter, rGeometries(i), nullptr);
            rIdCounter++;
            new_condition_list.push_back(p_condition);
        }

        rDestinationModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
    }

    ///@}
}