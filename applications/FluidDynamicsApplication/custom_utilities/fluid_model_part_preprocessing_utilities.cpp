//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <string>
#include <vector>

// External includes

// Project includes

// Application includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Include base h
#include "fluid_model_part_preprocessing_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template<class TContainerType>
class FluidModelPartPreProcessingUtilitiesHelperClass
{
public:

    using IndexType = std::size_t;

    using TEntityType = typename TContainerType::data_type;

    static TContainerType& GetContainer(ModelPart& rModelPart);

    static void AddEntities(ModelPart& rModelPart, const std::vector<IndexType>& rListOfIndices);

    static IndexType FindAndAddCommonEntities(
        ModelPart& rMainModelPart,
        ModelPart& rOutputModelPart,
        const std::vector<std::string>& rListOfInterfaceModelPartNames)
    {
        KRATOS_TRY

        auto& r_container = GetContainer(rMainModelPart);

        block_for_each(r_container, [](TEntityType& rEntity) {
            rEntity.SetValue(DOMAIN_SIZE, 0);
        });

        for (const auto& r_adjacent_model_part_name : rListOfInterfaceModelPartNames) {
            auto& r_adjacent_model_part = rMainModelPart.GetSubModelPart(r_adjacent_model_part_name);
            block_for_each(GetContainer(r_adjacent_model_part), [](TEntityType& rEntity) {
                rEntity.GetValue(DOMAIN_SIZE) += 1;
            });
        }

        std::vector<IndexType> list_of_indices;
        for (auto& r_entity : r_container) {
            if (r_entity.GetValue(DOMAIN_SIZE) > 1) {
                list_of_indices.push_back(r_entity.Id());
            }
        }

        AddEntities(rOutputModelPart, list_of_indices);

        return list_of_indices.size();

        KRATOS_CATCH("");
    }
};

template<> ModelPart::NodesContainerType& FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::NodesContainerType>::GetContainer(ModelPart& rModelPart) { return rModelPart.Nodes(); }
template<> ModelPart::ConditionsContainerType& FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::ConditionsContainerType>::GetContainer(ModelPart& rModelPart) { return rModelPart.Conditions(); }
template<> ModelPart::ElementsContainerType& FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::ElementsContainerType>::GetContainer(ModelPart& rModelPart) { return rModelPart.Elements(); }

template<> void FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::NodesContainerType>::AddEntities(ModelPart& rModelPart, const std::vector<IndexType>& rListOfIndices) { rModelPart.AddNodes(rListOfIndices); }
template<> void FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::ConditionsContainerType>::AddEntities(ModelPart& rModelPart, const std::vector<IndexType>& rListOfIndices) { rModelPart.AddConditions(rListOfIndices); }
template<> void FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::ElementsContainerType>::AddEntities(ModelPart& rModelPart, const std::vector<IndexType>& rListOfIndices) { rModelPart.AddElements(rListOfIndices); }

IndexType FluidModelPartPreProcessingUtilities::BreakElement(
    ModelPart& rModelPart,
    ModelPart::ElementType& rElement,
    IndexType& NewNodeId,
    IndexType& NewElementId,
    const std::string& rNewElementName)
{
    KRATOS_TRY

    auto& r_geometry = rElement.GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();

    array_1d<double, 3> new_coordinates{0.0, 0.0, 0.0};
    for (const auto& r_node : r_geometry) {
        new_coordinates += (r_node.Coordinates() / number_of_nodes);
    }

    auto p_new_node = rModelPart.CreateNewNode(NewNodeId++, new_coordinates[0], new_coordinates[1], new_coordinates[2]);

    auto p_properties = rElement.pGetProperties();
    std::vector<IndexType> list_of_node_ids(number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        for (IndexType j = 0; j < number_of_nodes - 1; ++j) {
            list_of_node_ids[j] = r_geometry[(i + j) % number_of_nodes].Id();
        }
        list_of_node_ids[number_of_nodes - 1] = p_new_node->Id();
        rModelPart.CreateNewElement(rNewElementName, NewElementId++, list_of_node_ids, p_properties);
    }

    return number_of_nodes;

    KRATOS_CATCH("");
}

bool FluidModelPartPreProcessingUtilities::CreateModelPartForCommenInterface(
    ModelPart& rModelPart,
    const std::string& rCommonInterfaceModelPartName,
    const std::vector<std::string>& rListOfInterfaceModelPartNames)
{
    KRATOS_TRY

    std::stringstream model_part_names;
    for (const auto& r_model_part_name : rListOfInterfaceModelPartNames) {
        model_part_names << r_model_part_name;
    }

    auto& r_output_model_part = rModelPart.CreateSubModelPart(rCommonInterfaceModelPartName);

    const IndexType number_of_common_nodes = FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::NodesContainerType>::FindAndAddCommonEntities(rModelPart, r_output_model_part, rListOfInterfaceModelPartNames);
    KRATOS_INFO("FluidModelPartPreProcessingUtilities") << "--- Found " << number_of_common_nodes << " common node(s) in " << rModelPart.FullName() << " using " << model_part_names.str() << "as interface model parts.";

    const IndexType number_of_common_conditions = FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::ConditionsContainerType>::FindAndAddCommonEntities(rModelPart, r_output_model_part, rListOfInterfaceModelPartNames);
    KRATOS_INFO("FluidModelPartPreProcessingUtilities") << "--- Found " << number_of_common_conditions << " common condition(s) in " << rModelPart.FullName() << " using " << model_part_names.str() << "as interface model parts.";

    const IndexType number_of_common_elements = FluidModelPartPreProcessingUtilitiesHelperClass<ModelPart::ElementsContainerType>::FindAndAddCommonEntities(rModelPart, r_output_model_part, rListOfInterfaceModelPartNames);
    KRATOS_INFO("FluidModelPartPreProcessingUtilities") << "--- Found " << number_of_common_elements << " common element(s) in " << rModelPart.FullName() << " using " << model_part_names.str() << "as interface model parts.";

    return !(number_of_common_nodes == 0 && number_of_common_conditions == 0 && number_of_common_elements == 0);

    KRATOS_CATCH("");
}

std::vector<IndexType> FluidModelPartPreProcessingUtilities::GetElementIdsWithAllNodesOnBoundaries(
    ModelPart& rModelPart,
    const std::vector<std::string>& rListOfBoundaryModelPartNames)
{
    KRATOS_TRY

    block_for_each(rModelPart.Nodes(), [](ModelPart::NodeType& rNode) {
        rNode.Set(VISITED, false);
    });

    for (const auto& r_boundary_model_part_name : rListOfBoundaryModelPartNames) {
        auto& r_boundary_model_part = rModelPart.GetSubModelPart(r_boundary_model_part_name);
        block_for_each(r_boundary_model_part.Nodes(), [](ModelPart::NodeType& rNode) {
            rNode.Set(VISITED, true);
        });
    }

    auto problematic_element_ids = block_for_each<AccumReduction<IndexType>>(rModelPart.Elements(), [](ModelPart::ElementType& rElement) -> IndexType {
        bool is_all_nodes_on_boundary = true;
        for (const auto& r_node : rElement.GetGeometry()) {
            if (!r_node.Is(VISITED)) {
                is_all_nodes_on_boundary = false;
                break;
            }
        }

        if (is_all_nodes_on_boundary) {
            return rElement.Id();
        } else {
            return 0;
        }
    });

    std::sort(problematic_element_ids.begin(), problematic_element_ids.end());
    problematic_element_ids.erase(problematic_element_ids.begin(), std::find_if(problematic_element_ids.begin(), problematic_element_ids.end(), [](const IndexType Value) { return Value > 0;}));

    return problematic_element_ids;

    KRATOS_CATCH("");
}

void FluidModelPartPreProcessingUtilities::BreakElements(
    ModelPart& rModelPart,
    const std::string& rNewElementName,
    const std::vector<IndexType>& rElementIds)
{
    KRATOS_TRY

    const IndexType local_max_node_id = block_for_each<MaxReduction<IndexType>>(rModelPart.Nodes(), [](ModelPart::NodeType& rNode) -> IndexType {
        return rNode.Id();
    });
    IndexType new_node_id = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_max_node_id) + 1;

    const IndexType local_max_element_id  = block_for_each<MaxReduction<IndexType>>(rModelPart.Elements(), [&](ModelPart::ElementType& rElement) -> IndexType {
        rElement.Set(TO_ERASE, false);
        return rElement.Id();
    });
    IndexType new_element_id = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_max_element_id) + 1;

    IndexType number_of_created_elements = 0;
    for (const auto& element_id : rElementIds) {
        auto& r_element = rModelPart.GetElement(element_id);
        r_element.Set(TO_ERASE, true);
        number_of_created_elements += BreakElement(rModelPart, r_element, new_node_id, new_element_id, rNewElementName);
    }

    rModelPart.RemoveElementsFromAllLevels();

    KRATOS_INFO("FluidModelPartPreProcessingUtilities") << "--- Created " << number_of_created_elements << " element(s) in " << rModelPart.FullName() << " while breaking elements with given ids.";

    KRATOS_CATCH("");
}


} // namespace Kratos

