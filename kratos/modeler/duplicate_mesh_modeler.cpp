//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "duplicate_mesh_modeler.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/timer.h"

namespace Kratos
{

void DuplicateMeshModeler::GenerateMesh(
    ModelPart& rThisModelPart,
    Element const& rReferenceElement,
    Condition const& rReferenceCondition)
{
    KRATOS_TRY;

    Timer::Start("Generating Mesh");

    Timer::Start("Generating Nodes");

    rThisModelPart.Nodes().clear();


    // Generate Copy of Nodes
    for (ModelPart::NodeIterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
    {
        rThisModelPart.CreateNewNode(i_node->Id(), *i_node);
    }

    Timer::Stop("Generating Nodes");

    rThisModelPart.PropertiesArray() = mrModelPart.PropertiesArray();

    Timer::Start("Generating Elements");


    //generating the elements
    const std::size_t element_size = rReferenceElement.GetGeometry().size();
    Element::NodesArrayType element_nodes_array(element_size);

    for(ModelPart::ElementsContainerType::iterator i_element = mrModelPart.ElementsBegin(); i_element != mrModelPart.ElementsEnd(); i_element++)
    {
        Element::GeometryType& geometry = i_element->GetGeometry();
        if(geometry.size() != element_size)
            KRATOS_THROW_ERROR(std::invalid_argument, "The given element is not compatible with the reference element", "");

        for(std::size_t i = 0 ; i < element_size ; i++)
            element_nodes_array(i) = rThisModelPart.pGetNode(geometry[i].Id());

        Element::Pointer p_element = rReferenceElement.Create(i_element->Id(), element_nodes_array, i_element->pGetProperties());
        rThisModelPart.Elements().push_back(p_element);
    }

    Timer::Stop("Generating Elements");

    Timer::Start("Generating Conditions");


    //generating the conditions
    const std::size_t condition_size = rReferenceCondition.GetGeometry().size();
    Condition::NodesArrayType condition_nodes_array(condition_size);

    for(ModelPart::ConditionsContainerType::iterator i_condition = mrModelPart.ConditionsBegin(); i_condition != mrModelPart.ConditionsEnd(); i_condition++)
    {
        Condition::GeometryType& geometry = i_condition->GetGeometry();
        if(geometry.size() != condition_size)
            KRATOS_THROW_ERROR(std::invalid_argument, "The given condition is not compatible with the reference condition", "");

        for(std::size_t i = 0 ; i < condition_size ; i++)
            condition_nodes_array(i) = rThisModelPart.pGetNode(geometry[i].Id());

        Condition::Pointer p_condition = rReferenceCondition.Create(i_condition->Id(), condition_nodes_array, i_condition->pGetProperties());
        rThisModelPart.Conditions().push_back(p_condition);
    }

    Timer::Stop("Generating Conditions");

    GenerateSubModelParts(mrModelPart, rThisModelPart);

    Timer::Stop("Generating Mesh");

    KRATOS_CATCH("");
}

void DuplicateMeshModeler::GenerateSubModelParts(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelpart)
{
    for (auto& r_ori_sub_model_part : rOriginModelPart.SubModelParts()) {
        const std::string sub_model_part_name = r_ori_sub_model_part.Name();
        auto& r_dest_sub_model_part = rDestinationModelpart.HasSubModelPart(sub_model_part_name) ?
            rDestinationModelpart.GetSubModelPart(sub_model_part_name) :
            rDestinationModelpart.CreateSubModelPart(sub_model_part_name);

        std::vector<std::size_t> new_nodes, new_elems, new_conds;

        new_nodes = block_for_each<AccumReduction<std::size_t>>(
            r_ori_sub_model_part.Nodes(),[](Node& rNode){
                return rNode.Id();
            });
        new_elems = block_for_each<AccumReduction<std::size_t>>(
            r_ori_sub_model_part.Elements(),[](Element& rElement){
                return rElement.Id();
            });
        new_conds = block_for_each<AccumReduction<std::size_t>>(
            r_ori_sub_model_part.Conditions(),[](Condition& rCondition){
                return rCondition.Id();
            });
        r_dest_sub_model_part.AddNodes(new_nodes);
        r_dest_sub_model_part.AddElements(new_elems);
        r_dest_sub_model_part.AddConditions(new_conds);

        GenerateSubModelParts(r_ori_sub_model_part, r_dest_sub_model_part);
    }

}

}
