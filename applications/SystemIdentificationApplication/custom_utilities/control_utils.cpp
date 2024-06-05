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

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "expression/literal_flat_expression.h"

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

// template instantiations
template void ControlUtils::AssignEquivalentProperties(ModelPart::ConditionsContainerType&, ModelPart::ConditionsContainerType&);
template void ControlUtils::AssignEquivalentProperties(ModelPart::ElementsContainerType&, ModelPart::ElementsContainerType&);

template void ControlUtils::ClipContainerExpression(ContainerExpression<ModelPart::NodesContainerType>&, const double, const double);
template void ControlUtils::ClipContainerExpression(ContainerExpression<ModelPart::ConditionsContainerType>&, const double, const double);
template void ControlUtils::ClipContainerExpression(ContainerExpression<ModelPart::ElementsContainerType>&, const double, const double);

} /* namespace Kratos.*/