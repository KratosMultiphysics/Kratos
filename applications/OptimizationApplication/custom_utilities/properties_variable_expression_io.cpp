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

// System incldues
#include <set>
#include <type_traits>

// Project includes
#include "expression/container_data_io.h"
#include "expression/expression_io_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/container_properties_data_io.h"

// Include base h
#include "properties_variable_expression_io.h"

namespace Kratos {

PropertiesVariableExpressionIO::Input::Input(
    const ModelPart& rModelPart,
    const VariableType& rVariable,
    Globals::DataLocation CurrentLocation)
    : mpModelPart(&rModelPart),
      mpVariable(rVariable),
      mDataLocation(CurrentLocation)
{
}

Expression::Pointer PropertiesVariableExpressionIO::Input::Execute() const
{
    switch (mDataLocation) {
        case Globals::DataLocation::Condition:
            return ExpressionIOUtils::ReadToExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::Properties>, const VariableType>(mpModelPart->Conditions(), mpVariable, mpModelPart->GetCommunicator().GetDataCommunicator());
        case Globals::DataLocation::Element:
            return ExpressionIOUtils::ReadToExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::Properties>, const VariableType>(mpModelPart->Elements(), mpVariable, mpModelPart->GetCommunicator().GetDataCommunicator());
        default:
            KRATOS_ERROR << "Invalid container type. Only supports Condition and Element.";
            break;
    }

    return nullptr;
}

PropertiesVariableExpressionIO::Output::Output(
    ModelPart& rModelPart,
    const VariableType& rVariable,
    Globals::DataLocation CurrentLocation)
    : mpModelPart(&rModelPart),
      mpVariable(rVariable),
      mDataLocation(CurrentLocation)
{
}

void PropertiesVariableExpressionIO::Output::Execute(const Expression& rExpression)
{
    KRATOS_TRY

    switch (mDataLocation) {
        case Globals::DataLocation::Condition:
            ExpressionIOUtils::WriteFromExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::Properties>, const VariableType>(mpModelPart->Conditions(), mpModelPart->GetCommunicator(), rExpression, mpVariable);
            break;
        case Globals::DataLocation::Element:
            ExpressionIOUtils::WriteFromExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::Properties>, const VariableType>(mpModelPart->Elements(), mpModelPart->GetCommunicator(), rExpression, mpVariable);
            break;
        default:
            KRATOS_ERROR << "Invalid container type. Only supports Condition and Element.";
            break;
    }

    KRATOS_CATCH("");
}

template<class TContainerType, MeshType TMeshType>
void PropertiesVariableExpressionIO::Read(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    static_assert(!std::is_same_v<TContainerType, ModelPart::NodesContainerType>,
                  "NodesContainerType expressions is not supported.\n");

    auto p_expression =
        Input(rContainerExpression.GetModelPart(), rVariable,
                                std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                    ? Globals::DataLocation::Condition
                                    : Globals::DataLocation::Element)
            .Execute();

    // p_expression is nullptr if there are no items in the ModelParts relevant container.
    // such as in ghost containers or interface containers.
    if (p_expression.get() != nullptr) {
        rContainerExpression.SetExpression(p_expression);
    }
}

template<class TContainerType, MeshType TMeshType>
void PropertiesVariableExpressionIO::Check(
    const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    static_assert(!std::is_same_v<TContainerType, ModelPart::NodesContainerType>,
                  "NodesContainerType expressions is not supported.\n");

    std::visit([&rContainerExpression](auto pVariable){
        using data_type = typename std::remove_const_t<std::remove_reference_t<decltype(*pVariable)>>::Type;

        const auto& value_ptrs = block_for_each<AccumReduction<data_type const*, std::set<data_type const*>>>(rContainerExpression.GetContainer(), [pVariable](const auto& rEntity) -> data_type const* {
            return &rEntity.GetProperties().GetValue(*pVariable);
        });

        const auto& r_data_communicator = rContainerExpression.GetModelPart().GetCommunicator().GetDataCommunicator();
        const auto number_of_values = r_data_communicator.SumAll(static_cast<unsigned int>(value_ptrs.size()));
        const auto number_of_entities = r_data_communicator.SumAll(static_cast<unsigned int>(rContainerExpression.GetContainer().size()));

        KRATOS_ERROR_IF_NOT(number_of_values == number_of_entities)
            << "Number of different values in properties for variable "
            << pVariable->Name() << " is different from number of entities in "
            << rContainerExpression.GetModelPart().FullName()
            << ". [ Number of different values in properties = " << number_of_values
            << ", number of entities = " << rContainerExpression.GetContainer().size()
            << " ].\n";

    }, rVariable);
}

template<class TContainerType, MeshType TMeshType>
void PropertiesVariableExpressionIO::Write(
    const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    static_assert(!std::is_same_v<TContainerType, ModelPart::NodesContainerType>,
                  "NodesContainerType expressions is not supported.\n");

    Output(*rContainerExpression.pGetModelPart(), rVariable,
                             std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                 ? Globals::DataLocation::Condition
                                 : Globals::DataLocation::Element)
        .Execute(rContainerExpression.GetExpression());
}

#define KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(CONTAINER_TYPE, MESH_TYPE)                                                                                      \
    template KRATOS_API(OPTIMIZATION_APPLICATION) void PropertiesVariableExpressionIO::Read(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const PropertiesVariableExpressionIO::VariableType&);          \
    template KRATOS_API(OPTIMIZATION_APPLICATION) void PropertiesVariableExpressionIO::Check(const ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const PropertiesVariableExpressionIO::VariableType&);   \
    template KRATOS_API(OPTIMIZATION_APPLICATION) void PropertiesVariableExpressionIO::Write(const ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const PropertiesVariableExpressionIO::VariableType&);

#define KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS_1(CONTAINER_TYPE)                \
    KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(CONTAINER_TYPE, MeshType::Local)

KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS_1(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS_1(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS
#undef KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS_1

} // namespace Kratos
