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
    ///@name Stages
    ///@{

    void Modeler::ImportGeometryModel(
        ModelPart& rModelPart,
        const std::string& rFileName,
        IndexType EchoLevel = 0) const
    {
        KRATOS_ERROR << "Calling ImportGeometryModel from the base modeler." << std::endl;
    }

    void Modeler::PrepareGeometryModel(
        ModelPart& rModelPart,
        const Parameters PrepareGeometryParameters,
        IndexType EchoLevel = 0) const
    {
        KRATOS_ERROR << "Calling ImportGeometryModel from the base modeler." << std::endl;
    }

    void Modeler::GenerateModelPart(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Parameters PrepareGeometryParameters,
        IndexType EchoLevel = 0) const
    {
        KRATOS_ERROR << "Calling GenerateModelPart from the base modeler." << std::endl;
    }

    void Modeler::ImportModelPart(
        ModelPart& rModelPart,
        const std::string& rFileName,
        IndexType EchoLevel = 0) const
    {
        KRATOS_ERROR << "Calling ImportModelPart from the base modeler." << std::endl;
    }

    void Modeler::PrepareModelPart(
        ModelPart& rModelPart,
        const Parameters PrepareModelPartParameters,
        IndexType EchoLevel = 0) const
    {
        KRATOS_ERROR << "Calling PrepareModelPart from the base modeler." << std::endl;
    }

    ///@}
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

    void Modeler::CreateElements(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        PropertiesPointerType pProperties,
        SizeType EchoLevel)
    {
        SizeType id = 1;
        if (rDestinationModelPart.GetRootModelPart().Elements().size() > 0)
            id = rDestinationModelPart.GetRootModelPart().Elements().back().Id() + 1;

        Modeler::CreateElements<ModelPart::GeometriesMapType>(
            rOriginModelPart.GeometriesBegin(),
            rOriginModelPart.GeometriesEnd(),
            rDestinationModelPart,
            rElementName, id, pProperties, EchoLevel);
    }

    template<class TContainerType>
    void Modeler::CreateElements(
        typename TContainerType::iterator rGeometriesBegin,
        typename TContainerType::iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        SizeType EchoLevel)
    {
        const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

        ElementsContainerType new_element_list;

        KRATOS_INFO_IF("CreateElements", EchoLevel > 2)
            << "Creating elements of type " << rElementName
            << " in " << rDestinationModelPart.Name() << "-SubModelPart." << std::endl;

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
        {
            new_element_list.push_back(
                rReferenceElement.Create(rIdCounter, *it, pProperties));
            rIdCounter++;
        }

        rDestinationModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }

    void Modeler::CreateConditions(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        PropertiesPointerType pProperties,
        SizeType EchoLevel)
    {
        SizeType id = 1;
        if (rDestinationModelPart.GetRootModelPart().Conditions().size() > 0)
            id = rDestinationModelPart.GetRootModelPart().Conditions().back().Id() + 1;

        Modeler::CreateConditions<ModelPart::GeometriesMapType>(
            rOriginModelPart.GeometriesBegin(),
            rOriginModelPart.GeometriesEnd(),
            rDestinationModelPart,
            rConditionName, id, pProperties, EchoLevel);
    }

    template<class TContainerType>
    void Modeler::CreateConditions(
        typename TContainerType::iterator rGeometriesBegin,
        typename TContainerType::iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        SizeType EchoLevel)
    {
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;

        KRATOS_INFO_IF("CreateConditions", EchoLevel > 2)
            << "Creating conditions of type " << rConditionName
            << " in " << rDestinationModelPart.Name() << "-SubModelPart." << std::endl;

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
        {
            new_condition_list.push_back(
                rReferenceCondition.Create(rIdCounter, *it, pProperties));
            rIdCounter++;
        }

        rDestinationModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
    }

    ///@}
}