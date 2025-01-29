//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <algorithm>
#include <tuple>
#include <type_traits>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "expression/literal_flat_expression.h"
#include "utilities/brute_force_point_locator.h"

// Application includes

// Include base h
#include "control_utils.h"

namespace Kratos {

template<class TContainerType>
void ControlUtils::AssignEquivalentProperties(
        TContainerType& rSourceContainer,
        TContainerType& rDestinationContainer)
{
    KRATOS_TRY

    // TODO: To be removed once the PointerVectorSet mutable find is fixed.
    rSourceContainer.Sort();
    rDestinationContainer.Sort();

    const IndexType number_of_entities = rSourceContainer.size();

    KRATOS_ERROR_IF_NOT(number_of_entities == rDestinationContainer.size())
        << "Number of entities mismatch [ rSourceContainer.size() = " << number_of_entities
        << ", rDestinationContainer.size() = " << rDestinationContainer.size() << " ].\n";

    IndexPartition<IndexType>(number_of_entities).for_each([&rSourceContainer, &rDestinationContainer](const auto Index) {
        auto& r_destination_entity = *(rDestinationContainer.begin() + Index);
        auto p_itr = rSourceContainer.find(r_destination_entity.Id());

        KRATOS_ERROR_IF(p_itr == rSourceContainer.end()) <<
            "The entity with id = " << r_destination_entity.Id() << " not found in the source container.";

        r_destination_entity.SetProperties(p_itr->pGetProperties());
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
void ControlUtils::ClipContainerExpression(
    ContainerExpression<TContainerType>& rContainerExpression,
    const double Min,
    const double Max)
{
    KRATOS_TRY

    const auto& r_expression = rContainerExpression.GetExpression();
    const auto number_of_entities = r_expression.NumberOfEntities();
    const auto number_of_components = r_expression.GetItemComponentCount();
    auto p_output_expression = LiteralFlatExpression<double>::Create(number_of_entities, r_expression.GetItemShape());

    IndexPartition<IndexType>(number_of_entities).for_each([&p_output_expression, &r_expression, number_of_components, Min, Max](const auto Index) {
        const auto data_begin_index = Index * number_of_components;
        for (IndexType i = 0; i < number_of_components; ++i) {
            p_output_expression->SetData(data_begin_index, i, std::clamp(r_expression.Evaluate(Index, data_begin_index, i), Min, Max));
        }
    });

    rContainerExpression.SetExpression(p_output_expression);

    KRATOS_CATCH("");
}

template<class TContainerType>
void ControlUtils::GetIntegrationPoints(
    std::vector<Point>& rOutput,
    const TContainerType& rContainer)
{
    KRATOS_TRY

    rOutput.clear();
    rOutput.reserve(rContainer.size() * 5);

    for (const auto& r_entity : rContainer) {
        const auto& r_geometry = r_entity.GetGeometry();
        const Matrix& shape_function_values_list = r_geometry.ShapeFunctionsValues(r_entity.GetIntegrationMethod());
        const IndexType number_of_integration_points = shape_function_values_list.size1();

        for (IndexType i_gauss = 0; i_gauss < number_of_integration_points; ++i_gauss) {
            const Vector& shape_function_values = row(shape_function_values_list, i_gauss);

            double x{0}, y{0}, z{0};
            for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
                x += r_geometry[i_node].X() * shape_function_values[i_node];
                y += r_geometry[i_node].Y() * shape_function_values[i_node];
                z += r_geometry[i_node].Z() * shape_function_values[i_node];
            }

            rOutput.push_back(Point(x, y, z));
        }
    }

    KRATOS_CATCH("");
}

template<class EntityType, class TDataType>
void ControlUtils::EvaluateAtPoints(
    std::vector<TDataType>& rOutput,
    const Variable<TDataType>& rVariable,
    ModelPart& rModelPart,
    const std::vector<Point>& rCoordinates)
{
    KRATOS_TRY

    if (rOutput.size() != rCoordinates.size()) {
        rOutput.resize(rCoordinates.size());
    }

    BruteForcePointLocator point_locator(rModelPart);

    IndexPartition<IndexType>(rCoordinates.size()).for_each(Vector(), [&point_locator, &rCoordinates, &rModelPart, &rVariable, &rOutput](const auto Index, auto& rTLS) {
        if constexpr(std::is_same_v<EntityType, Node>) {
            const auto entity_id = point_locator.FindElement(*(rCoordinates.begin() + Index), rTLS);
            const auto& r_geometry = rModelPart.GetElement(entity_id).GetGeometry();

            auto value = rVariable.Zero();
            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                value += r_geometry[i].GetValue(rVariable) * rTLS[i];
            }

            rOutput[Index] = value;
        } else if constexpr(std::is_same_v<EntityType, Element>) {
            const auto entity_id = point_locator.FindElement(*(rCoordinates.begin() + Index), rTLS);
            rOutput[Index] = rModelPart.GetElement(entity_id).GetValue(rVariable);
        } else if constexpr(std::is_same_v<EntityType, Condition>) {
            const auto entity_id = point_locator.FindCondition(*(rCoordinates.begin() + Index), rTLS);
            rOutput[Index] = rModelPart.GetCondition(entity_id).GetValue(rVariable);
        }
    });


    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::AssignEquivalentProperties(ModelPart::ConditionsContainerType&, ModelPart::ConditionsContainerType&);
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::AssignEquivalentProperties(ModelPart::ElementsContainerType&, ModelPart::ElementsContainerType&);

template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::ClipContainerExpression(ContainerExpression<ModelPart::NodesContainerType>&, const double, const double);
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::ClipContainerExpression(ContainerExpression<ModelPart::ConditionsContainerType>&, const double, const double);
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::ClipContainerExpression(ContainerExpression<ModelPart::ElementsContainerType>&, const double, const double);

template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::GetIntegrationPoints(std::vector<Point>&, const ModelPart::ConditionsContainerType&);
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::GetIntegrationPoints(std::vector<Point>&, const ModelPart::ElementsContainerType&);

template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::EvaluateAtPoints<Node, double>(std::vector<double>&, const Variable<double>&, ModelPart&, const std::vector<Point>&);
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::EvaluateAtPoints<Element, double>(std::vector<double>&, const Variable<double>&, ModelPart&, const std::vector<Point>&);
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void ControlUtils::EvaluateAtPoints<Condition, double>(std::vector<double>&, const Variable<double>&, ModelPart&, const std::vector<Point>&);

} /* namespace Kratos.*/