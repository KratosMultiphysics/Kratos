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

Expression::Pointer VariableExpressionIO::VariableExpressionInput::Execute() const
{
    const auto& r_mesh = GetMesh(mrModelPart.GetCommunicator(), mMeshType);

    switch (mContainerType) {
        case ContainerType::NodalHistorical: {
                return VariableExpressionIOUtils::ReadToExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>, const VariableType>(r_mesh.Nodes(), mpVariable);
                break;
            }
        case ContainerType::NodalNonHistorical: {
                return VariableExpressionIOUtils::ReadToExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Nodes(), mpVariable);
                break;
            }
        case ContainerType::ConditionNonHistorical: {
                return VariableExpressionIOUtils::ReadToExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Conditions(), mpVariable);
                break;
            }
        case ContainerType::ElementNonHistorical: {
                return VariableExpressionIOUtils::ReadToExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>, const VariableType>(r_mesh.Elements(), mpVariable);
                break;
            }
        default: {
                KRATOS_ERROR << "Invalid container type";
            break;
        }
    }

    return nullptr;
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

template<MeshType TMeshType>
void VariableExpressionIO::Read(
    ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable,
    const bool IsHistorical)
{
    auto p_expression =
        VariableExpressionInput(rContainerExpression.GetModelPart(), rVariable,
                                IsHistorical ? ContainerType::NodalHistorical
                                             : ContainerType::NodalNonHistorical,
                                TMeshType)
            .Execute();

    // p_expression is nullptr if there are no items in the ModelParts relevant container.
    // such as in ghost containers or interface containers.
    if (p_expression.get() != nullptr) {
        rContainerExpression.SetExpression(p_expression);
    }
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
        VariableExpressionInput(rContainerExpression.GetModelPart(), rVariable,
                                std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                    ? ContainerType::ConditionNonHistorical
                                    : ContainerType::ElementNonHistorical,
                                TMeshType)
            .Execute();

    // p_expression is nullptr if there are no items in the ModelParts relevant container.
    // such as in ghost containers or interface containers.
    if (p_expression.get() != nullptr) {
        rContainerExpression.SetExpression(p_expression);
    }
}

template<MeshType TMeshType>
void VariableExpressionIO::Write(
    const ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable,
    const bool IsHistorical)
{
    VariableExpressionOutput(*rContainerExpression.pGetModelPart(), rVariable,
                             IsHistorical ? ContainerType::NodalHistorical
                                          : ContainerType::NodalNonHistorical,
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

    VariableExpressionOutput(*rContainerExpression.pGetModelPart(), rVariable,
                             std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                 ? ContainerType::ConditionNonHistorical
                                 : ContainerType::ElementNonHistorical,
                             TMeshType)
        .Execute(rContainerExpression.GetExpression());
}

ModelPart::MeshType& VariableExpressionIO::GetMesh(
    Communicator& rCommunicator,
    MeshType  rMeshType)
{
    switch (rMeshType) {
        case MeshType::Local: {
                return rCommunicator.LocalMesh();
                break;
            }
        case MeshType::Interface: {
                return rCommunicator.InterfaceMesh();
                break;
            }
        case MeshType::Ghost: {
                return rCommunicator.GhostMesh();
                break;
            }
        default: {
            KRATOS_ERROR << "Invalid mesh type";
            break;
        }
    }
}

const ModelPart::MeshType& VariableExpressionIO::GetMesh(
    const Communicator& rCommunicator,
    MeshType  rMeshType)
{
    switch (rMeshType) {
        case MeshType::Local: {
                return rCommunicator.LocalMesh();
                break;
            }
        case MeshType::Interface: {
                return rCommunicator.InterfaceMesh();
                break;
            }
        case MeshType::Ghost: {
                return rCommunicator.GhostMesh();
                break;
            }
        default: {
            KRATOS_ERROR << "Invalid mesh type";
            break;
        }
    }
}

#define KRATOS_INSTANTIATE_NODAL_CONTAINER_IO_METHODS(MESH_TYPE)                                                                                                \
    template void VariableExpressionIO::Read(ContainerExpression<ModelPart::NodesContainerType, MESH_TYPE>&, const VariableExpressionIO::VariableType&, const bool); \
    template void VariableExpressionIO::Write(const ContainerExpression<ModelPart::NodesContainerType, MESH_TYPE>&, const VariableExpressionIO::VariableType&, const bool);\

#define KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(CONTAINER_TYPE, MESH_TYPE)                                                    \
    template void VariableExpressionIO::Read(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const VariableExpressionIO::VariableType&); \
    template void VariableExpressionIO::Write(const ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const VariableExpressionIO::VariableType&);\

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
