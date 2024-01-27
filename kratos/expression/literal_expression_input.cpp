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
#include <type_traits>

// Project includes
#include "expression/literal_expression.h"
#include "expression/expression_io_utils.h"

// Include base h
#include "literal_expression_input.h"

namespace Kratos {

LiteralExpressionIO::Input::Input(
    const ModelPart& rModelPart,
    const DataType& rValue,
    Globals::DataLocation CurrentLocation,
    MeshType CurrentMeshType)
    : mpModelPart(&rModelPart),
      mValue(rValue),
      mDataLocation(CurrentLocation),
      mMeshType(CurrentMeshType)
{
}
Expression::Pointer LiteralExpressionIO::Input::Execute() const
{
    return std::visit([&](auto& rValue){
        using data_type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;

        IndexType number_of_entities = 0;
        switch (mDataLocation) {
            case Globals::DataLocation::NodeHistorical:
            case Globals::DataLocation::NodeNonHistorical:
                number_of_entities = ExpressionIOUtils::GetMesh(mpModelPart->GetCommunicator(), mMeshType).NumberOfNodes();
                break;
            case Globals::DataLocation::Condition:
                number_of_entities = ExpressionIOUtils::GetMesh(mpModelPart->GetCommunicator(), mMeshType).NumberOfConditions();
                break;
            case Globals::DataLocation::Element:
                number_of_entities = ExpressionIOUtils::GetMesh(mpModelPart->GetCommunicator(), mMeshType).NumberOfElements();
                break;
            default:
                KRATOS_ERROR << "Invalid container type. Only supports NodeHistorical, NodeNonHistorical, Condition, Element.";
                break;
        }

        return LiteralExpression<data_type>::Create(rValue, number_of_entities);
    }, mValue);
}

template<class TContainerType, MeshType TMeshType>
void LiteralExpressionIO::SetData(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const DataType& rValue)
{
    rContainerExpression.SetExpression(
        Input(rContainerExpression.GetModelPart(), rValue,
                            std::is_same_v<TContainerType, ModelPart::NodesContainerType>
                                ? Globals::DataLocation::NodeNonHistorical
                                : std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                    ? Globals::DataLocation::Condition
                                    : Globals::DataLocation::Element,
                            TMeshType)
            .Execute());
}

template<class TContainerType, MeshType TMeshType>
void LiteralExpressionIO::SetDataToZero(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    std::visit([&rContainerExpression](const auto pVariable) {
        rContainerExpression.SetExpression(
            Input(rContainerExpression.GetModelPart(), pVariable->Zero(),
                                std::is_same_v<TContainerType, ModelPart::NodesContainerType>
                                    ? Globals::DataLocation::NodeNonHistorical
                                    : std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                        ? Globals::DataLocation::Condition
                                        : Globals::DataLocation::Element,
                                TMeshType)
                .Execute());
    }, rVariable);
}

// template instantiations

#define KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS(CONTAINER_TYPE, MESH_TYPE)                                                                                  \
    template KRATOS_API(KRATOS_CORE) void LiteralExpressionIO::SetData(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const LiteralExpressionIO::DataType&);          \
    template KRATOS_API(KRATOS_CORE) void LiteralExpressionIO::SetDataToZero(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const LiteralExpressionIO::VariableType&);

#define KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS_1(MESH_TYPE)                                  \
    KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS(ModelPart::NodesContainerType, MESH_TYPE)         \
    KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS(ModelPart::ConditionsContainerType, MESH_TYPE)    \
    KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS(ModelPart::ElementsContainerType, MESH_TYPE)

KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS_1(MeshType::Local)
KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS_1(MeshType::Interface)
KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS_1(MeshType::Ghost)

#undef KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS_1
#undef KRATOS_INSTANTIATE_DATA_EXPRESSION_IO_METHODS

} // namespace Kratos
