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
#include "variable_expression_io_utils.h"

// Include base h
#include "variable_expression_io.h"

namespace Kratos {

template<class TDataType>
VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const ContainerType& rContainerType,
    MeshType rMeshType)
    : mrModelPart(rModelPart),
      mpVariable(&rVariable),
      mContainerType(rContainerType),
      mMeshType(rMeshType)
{
}

VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(
    const ModelPart& rModelPart,
    const VariableType& rVariable,
    const ContainerType& rContainerType,
    MeshType rMeshType)
    : mrModelPart(rModelPart),
      mpVariable(rVariable),
      mContainerType(rContainerType),
      mMeshType(rMeshType)
{
}

template <class TDataType, MeshType TMeshType>
VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(
    const ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainer,
    const Variable<TDataType>& rVariable,
    const bool IsHistorical)
    : VariableExpressionInput(
        rContainer.GetModelPart(),
        rVariable,
        IsHistorical ? ContainerType::NodalHistorical : ContainerType::NodalNonHistorical,
        TMeshType == MeshType::Local ? MeshType::Local : TMeshType == MeshType::Interface ? MeshType::Interface : MeshType::Ghost)
{

}

template <class TContainerType, class TDataType, MeshType TMeshType>
VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(
    const ContainerExpression<TContainerType, TMeshType>& rContainer,
    const Variable<TDataType>& rVariable)
    : VariableExpressionInput(
        rContainer.GetModelPart(),
        rVariable,
        std::is_same_v<TContainerType, ModelPart::ConditionsContainerType> ? ContainerType::ConditionNonHistorical : ContainerType::ElementNonHistorical,
        TMeshType == MeshType::Local ? MeshType::Local : TMeshType == MeshType::Interface ? MeshType::Interface : MeshType::Ghost)
{
}

Expression::Pointer VariableExpressionIO::VariableExpressionInput::Execute() const
{
    const auto& r_mesh = GetMesh(mrModelPart.GetCommunicator(), mMeshType);

    switch (mContainerType) {
        case ContainerType::NodalHistorical:
            return VariableExpressionIOUtils::ReadToExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>, const VariableType>(r_mesh.Nodes(), mpVariable);
        case ContainerType::NodalNonHistorical:
            return VariableExpressionIOUtils::ReadToExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Nodes(), mpVariable);
        case ContainerType::ConditionNonHistorical:
            return VariableExpressionIOUtils::ReadToExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Conditions(), mpVariable);
        case ContainerType::ElementNonHistorical:
            return VariableExpressionIOUtils::ReadToExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Elements(), mpVariable);
    }

    return nullptr;
}

template<class TDataType>
VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const ContainerType& rContainerType,
    MeshType  rMeshType)
    : mrModelPart(rModelPart),
      mpVariable(&rVariable),
      mContainerType(rContainerType),
      mMeshType(rMeshType)
{
}

VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(
    ModelPart& rModelPart,
    const VariableType& rVariable,
    const ContainerType& rContainerType,
    MeshType  rMeshType)
    : mrModelPart(rModelPart),
      mpVariable(rVariable),
      mContainerType(rContainerType),
      mMeshType(rMeshType)
{
}

template <class TDataType, MeshType TMeshType>
VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(
    ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainer,
    const Variable<TDataType>& rVariable,
    const bool IsHistorical)
    : VariableExpressionOutput(
        rContainer.GetModelPart(),
        rVariable,
        IsHistorical ? ContainerType::NodalHistorical : ContainerType::NodalNonHistorical,
        TMeshType == MeshType::Local ? MeshType::Local : TMeshType == MeshType::Interface ? MeshType::Interface : MeshType::Ghost)
{

}

template <class TContainerType, class TDataType, MeshType TMeshType>
VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(
    ContainerExpression<TContainerType, TMeshType>& rContainer,
    const Variable<TDataType>& rVariable)
    : VariableExpressionOutput(
        rContainer.GetModelPart(),
        rVariable,
        std::is_same_v<TContainerType, ModelPart::ConditionsContainerType> ? ContainerType::ConditionNonHistorical : ContainerType::ElementNonHistorical,
        TMeshType == MeshType::Local ? MeshType::Local : TMeshType == MeshType::Interface ? MeshType::Interface : MeshType::Ghost)
{
}

void VariableExpressionIO::VariableExpressionOutput::Execute(const Expression& rExpression)
{
    KRATOS_TRY
    auto& r_communicator = mrModelPart.GetCommunicator();
    auto& r_mesh = GetMesh(r_communicator, mMeshType);

    switch (mContainerType) {
        case ContainerType::NodalHistorical:
            VariableExpressionIOUtils::WriteFromExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>, const VariableType>(r_mesh.Nodes(), r_communicator, rExpression, mpVariable);
            break;
        case ContainerType::NodalNonHistorical:
            VariableExpressionIOUtils::WriteFromExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Nodes(), r_communicator, rExpression, mpVariable);
            break;
        case ContainerType::ConditionNonHistorical:
            VariableExpressionIOUtils::WriteFromExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Conditions(), r_communicator, rExpression, mpVariable);
            break;
        case ContainerType::ElementNonHistorical:
            VariableExpressionIOUtils::WriteFromExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Elements(), r_communicator, rExpression, mpVariable);
            break;
    }
    KRATOS_CATCH("");
}

ModelPart::MeshType& VariableExpressionIO::GetMesh(
    Communicator& rCommunicator,
    MeshType  rMeshType)
{
    switch (rMeshType) {
        case MeshType::Local:
            return rCommunicator.LocalMesh();
        case MeshType::Interface:
            return rCommunicator.InterfaceMesh();
        case MeshType::Ghost:
            return rCommunicator.GhostMesh();
    }

    return rCommunicator.LocalMesh();
}

const ModelPart::MeshType& VariableExpressionIO::GetMesh(
    const Communicator& rCommunicator,
    MeshType  rMeshType)
{
    switch (rMeshType) {
        case MeshType::Local:
            return rCommunicator.LocalMesh();
        case MeshType::Interface:
            return rCommunicator.InterfaceMesh();
        case MeshType::Ghost:
            return rCommunicator.GhostMesh();
    }

    return rCommunicator.LocalMesh();
}

// template instantiations
#define KRATOS_VARIABLE_EXPRESSION_IO_NODES(...)                                                                                                                                                                        \
    template VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(const ContainerExpression<ModelPart::NodesContainerType, MeshType::Local>&, const Variable<__VA_ARGS__>&, const bool);      \
    template VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(const ContainerExpression<ModelPart::NodesContainerType, MeshType::Interface>&, const Variable<__VA_ARGS__>&, const bool);  \
    template VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(const ContainerExpression<ModelPart::NodesContainerType, MeshType::Ghost>&, const Variable<__VA_ARGS__>&, const bool);      \
    template VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(ContainerExpression<ModelPart::NodesContainerType, MeshType::Local>&, const Variable<__VA_ARGS__>&, const bool);          \
    template VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(ContainerExpression<ModelPart::NodesContainerType, MeshType::Interface>&, const Variable<__VA_ARGS__>&, const bool);      \
    template VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(ContainerExpression<ModelPart::NodesContainerType, MeshType::Ghost>&, const Variable<__VA_ARGS__>&, const bool);

#define KRATOS_VARIABLE_EXPRESSION_IO_ENTITIES(ENTITY_TYPE, ...)                                                                                                                            \
    template VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(const ContainerExpression<ENTITY_TYPE, MeshType::Local>&, const Variable<__VA_ARGS__>&);        \
    template VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(const ContainerExpression<ENTITY_TYPE, MeshType::Interface>&, const Variable<__VA_ARGS__>&);    \
    template VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(const ContainerExpression<ENTITY_TYPE, MeshType::Ghost>&, const Variable<__VA_ARGS__>&);        \
    template VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(ContainerExpression<ENTITY_TYPE, MeshType::Local>&, const Variable<__VA_ARGS__>&);            \
    template VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(ContainerExpression<ENTITY_TYPE, MeshType::Interface>&, const Variable<__VA_ARGS__>&);        \
    template VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(ContainerExpression<ENTITY_TYPE, MeshType::Ghost>&, const Variable<__VA_ARGS__>&);

#define KRATOS_VARIABLE_EXPRESSION_IO(...)                                                                                                                                                                              \
    template VariableExpressionIO::VariableExpressionInput::VariableExpressionInput(const ModelPart&, const Variable<__VA_ARGS__>&, const VariableExpressionIO::ContainerType&, MeshType); \
    template VariableExpressionIO::VariableExpressionOutput::VariableExpressionOutput(ModelPart&, const Variable<__VA_ARGS__>&, const VariableExpressionIO::ContainerType&, MeshType);     \
    KRATOS_VARIABLE_EXPRESSION_IO_NODES(__VA_ARGS__)                                                                                                                                                                    \
    KRATOS_VARIABLE_EXPRESSION_IO_ENTITIES(ModelPart::ConditionsContainerType, __VA_ARGS__)                                                                                                                             \
    KRATOS_VARIABLE_EXPRESSION_IO_ENTITIES(ModelPart::ElementsContainerType, __VA_ARGS__)

KRATOS_VARIABLE_EXPRESSION_IO(int)
KRATOS_VARIABLE_EXPRESSION_IO(double)
KRATOS_VARIABLE_EXPRESSION_IO(array_1d<double, 3>)
KRATOS_VARIABLE_EXPRESSION_IO(array_1d<double, 4>)
KRATOS_VARIABLE_EXPRESSION_IO(array_1d<double, 6>)
KRATOS_VARIABLE_EXPRESSION_IO(array_1d<double, 9>)
KRATOS_VARIABLE_EXPRESSION_IO(Vector)
KRATOS_VARIABLE_EXPRESSION_IO(Matrix)

#undef KRATOS_VARIABLE_EXPRESSION_IO
#undef KRATOS_VARIABLE_EXPRESSION_IO_ENTITIES
#undef KRATOS_VARIABLE_EXPRESSION_IO_NODES

} // namespace Kratos
