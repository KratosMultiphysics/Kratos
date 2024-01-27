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
#include "expression/container_data_io.h"
#include "expression_io_utils.h"

// Include base h
#include "variable_expression_io.h"

namespace Kratos {

VariableExpressionIO::Input::Input(
    const ModelPart& rModelPart,
    const VariableType& rVariable,
    Globals::DataLocation CurrentLocation,
    MeshType CurrentMeshType)
    : mpModelPart(&rModelPart),
      mpVariable(rVariable),
      mDataLocation(CurrentLocation),
      mMeshType(CurrentMeshType)
{
}

Expression::Pointer VariableExpressionIO::Input::Execute() const
{
    const auto& r_mesh = ExpressionIOUtils::GetMesh(mpModelPart->GetCommunicator(), mMeshType);
    const auto& r_data_communicator = mpModelPart->GetCommunicator().GetDataCommunicator();

    switch (mDataLocation) {
        case Globals::DataLocation::NodeHistorical:
            return ExpressionIOUtils::ReadToExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>, const VariableType>(r_mesh.Nodes(), mpVariable, r_data_communicator);
        case Globals::DataLocation::NodeNonHistorical:
            return ExpressionIOUtils::ReadToExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Nodes(), mpVariable, r_data_communicator);
        case Globals::DataLocation::Condition:
            return ExpressionIOUtils::ReadToExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Conditions(), mpVariable, r_data_communicator);
        case Globals::DataLocation::Element:
            return ExpressionIOUtils::ReadToExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Elements(), mpVariable, r_data_communicator);
        default:
            KRATOS_ERROR << "Invalid container type. Only supports NodeHistorical, NodeNonHistorical, Condition, Element.";
            break;
    }

    return nullptr;
}

VariableExpressionIO::Output::Output(
    ModelPart& rModelPart,
    const VariableType& rVariable,
    Globals::DataLocation CurrentLocation,
    MeshType CurrentMeshType)
    : mpModelPart(&rModelPart),
      mpVariable(rVariable),
      mDataLocation(CurrentLocation),
      mMeshType(CurrentMeshType)
{
}

void VariableExpressionIO::Output::Execute(const Expression& rExpression)
{
    KRATOS_TRY
    auto& r_communicator = mpModelPart->GetCommunicator();
    auto& r_mesh = ExpressionIOUtils::GetMesh(r_communicator, mMeshType);

    switch (mDataLocation) {
        case Globals::DataLocation::NodeHistorical:
            ExpressionIOUtils::WriteFromExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>, const VariableType>(r_mesh.Nodes(), r_communicator, rExpression, mpVariable);
            break;
        case Globals::DataLocation::NodeNonHistorical:
            ExpressionIOUtils::WriteFromExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Nodes(), r_communicator, rExpression, mpVariable);
            break;
        case Globals::DataLocation::Condition:
            ExpressionIOUtils::WriteFromExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Conditions(), r_communicator, rExpression, mpVariable);
            break;
        case Globals::DataLocation::Element:
            ExpressionIOUtils::WriteFromExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Elements(), r_communicator, rExpression, mpVariable);
            break;
        default:
            KRATOS_ERROR << "Invalid container type. Only supports NodeHistorical, NodeNonHistorical, Condition, Element.";
            break;
    }
    KRATOS_CATCH("");
}

template<MeshType TMeshType>
void VariableExpressionIO::Read(
    ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable,
    const bool IsHistorical)
{
    auto p_expression =
        Input(rContainerExpression.GetModelPart(), rVariable,
                                IsHistorical ? Globals::DataLocation::NodeHistorical
                                             : Globals::DataLocation::NodeNonHistorical,
                                TMeshType)
            .Execute();

    rContainerExpression.SetExpression(p_expression);
}

template<class TContainerType, MeshType TMeshType>
void VariableExpressionIO::Read(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    static_assert(!std::is_same_v<TContainerType, ModelPart::NodesContainerType>,
                  "NodesContainerType expressions should have the IsHistorical "
                  "stated.\n");

    auto p_expression =
        Input(rContainerExpression.GetModelPart(), rVariable,
                                std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                    ? Globals::DataLocation::Condition
                                    : Globals::DataLocation::Element,
                                TMeshType)
            .Execute();

    rContainerExpression.SetExpression(p_expression);
}

template<MeshType TMeshType>
void VariableExpressionIO::Write(
    const ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable,
    const bool IsHistorical)
{
    Output(*rContainerExpression.pGetModelPart(), rVariable,
                             IsHistorical ? Globals::DataLocation::NodeHistorical
                                          : Globals::DataLocation::NodeNonHistorical,
                             TMeshType)
        .Execute(rContainerExpression.GetExpression());
}

template<class TContainerType, MeshType TMeshType>
void VariableExpressionIO::Write(
    const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    static_assert(!std::is_same_v<TContainerType, ModelPart::NodesContainerType>,
                  "NodesContainerType expressions should have the IsHistorical "
                  "stated.\n");

    Output(*rContainerExpression.pGetModelPart(), rVariable,
                             std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                 ? Globals::DataLocation::Condition
                                 : Globals::DataLocation::Element,
                             TMeshType)
        .Execute(rContainerExpression.GetExpression());
}

#define KRATOS_INSTANTIATE_NODAL_CONTAINER_IO_METHODS(MESH_TYPE)                                                                                                \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionIO::Read(ContainerExpression<ModelPart::NodesContainerType, MESH_TYPE>&, const VariableExpressionIO::VariableType&, const bool); \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionIO::Write(const ContainerExpression<ModelPart::NodesContainerType, MESH_TYPE>&, const VariableExpressionIO::VariableType&, const bool);\

#define KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(CONTAINER_TYPE, MESH_TYPE)                                                    \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionIO::Read(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const VariableExpressionIO::VariableType&); \
    template KRATOS_API(KRATOS_CORE) void VariableExpressionIO::Write(const ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const VariableExpressionIO::VariableType&);\

#define KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MESH_TYPE)                              \
    KRATOS_INSTANTIATE_NODAL_CONTAINER_IO_METHODS(MESH_TYPE)                                        \
    KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(ModelPart::ConditionsContainerType, MESH_TYPE)   \
    KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(ModelPart::ElementsContainerType, MESH_TYPE)

KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MeshType::Local)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MeshType::Interface)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MeshType::Ghost)

#undef KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO
#undef KRATOS_INSTANTIATE_NODAL_CONTAINER_IO_METHODS
#undef KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS

} // namespace Kratos
