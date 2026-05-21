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

// Project includes
#include "expression/expression_io_utils.h"
#include "expression/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "nodal_position_expression_io.h"

namespace Kratos {

NodalPositionExpressionIO::Input::Input(
    const ModelPart& rModelPart,
    const Configuration& rConfiguration,
    const MeshType& rMeshType)
    : mrModelPart(rModelPart),
      mConfiguration(rConfiguration),
      mMeshType(rMeshType)
{
}

Expression::Pointer NodalPositionExpressionIO::Input::Execute() const
{
    const auto& r_mesh = ExpressionIOUtils::GetMesh(mrModelPart.GetCommunicator(), mMeshType);
    const auto number_of_nodes = r_mesh.Nodes().size();

    LiteralFlatExpression<double>::Pointer expression;

    switch (mConfiguration) {
        case Configuration::Initial:
            expression = LiteralFlatExpression<double>::Create(number_of_nodes, {3});
            IndexPartition<IndexType>(number_of_nodes).for_each([&r_mesh, &expression](const IndexType Index) {
                const auto& r_node = *(r_mesh.NodesBegin() + Index);
                double* current_data_start = expression->begin() + Index * 3;
                *current_data_start++ = r_node.X0();
                *current_data_start++ = r_node.Y0();
                *current_data_start = r_node.Z0();
            });
            break;
        case Configuration::Current:
            expression = LiteralFlatExpression<double>::Create(number_of_nodes, {3});
            IndexPartition<IndexType>(number_of_nodes).for_each([&r_mesh, &expression](const IndexType Index) {
                const auto& r_node = *(r_mesh.NodesBegin() + Index);
                double* current_data_start = expression->begin() + Index * 3;
                *current_data_start++ = r_node.X();
                *current_data_start++ = r_node.Y();
                *current_data_start = r_node.Z();
            });
            break;
    }

    return expression;
}

NodalPositionExpressionIO::Output::Output(
    ModelPart& rModelPart,
    const Configuration& rConfiguration,
    const MeshType& rMeshType)
    : mrModelPart(rModelPart),
      mConfiguration(rConfiguration),
      mMeshType(rMeshType)
{
}

void NodalPositionExpressionIO::Output::Execute(const Expression& rExpression)
{
    KRATOS_TRY

    auto& r_communicator = mrModelPart.GetCommunicator();
    auto& r_mesh = ExpressionIOUtils::GetMesh(r_communicator, mMeshType);
    auto& r_nodes = r_mesh.Nodes();

    KRATOS_ERROR_IF_NOT(rExpression.NumberOfEntities() == r_mesh.Nodes().size())
        << "Local nodes size and expression number of entities mismatch. [ local nodes size = "
        << r_nodes.size() << ", expression number of entities = " << rExpression.NumberOfEntities()
        << ", expression = " << rExpression << " ].\n";

    KRATOS_ERROR_IF_NOT(rExpression.GetItemComponentCount() == 3)
        << "Shape mismatch. [ Required shape of the expression = [3], given expression shape = "
        << rExpression.GetItemShape() << ", expression = " << rExpression << " ].\n";

    const auto number_of_nodes = r_mesh.Nodes().size();

    switch (mConfiguration) {
        case Configuration::Initial:
            IndexPartition<IndexType>(number_of_nodes).for_each([&r_mesh, &rExpression](const IndexType Index) {
                const IndexType entity_data_start_index = Index * 3;
                auto& r_node = *(r_mesh.NodesBegin() + Index);
                r_node.X0() = rExpression.Evaluate(Index, entity_data_start_index, 0);
                r_node.Y0() = rExpression.Evaluate(Index, entity_data_start_index, 1);
                r_node.Z0() = rExpression.Evaluate(Index, entity_data_start_index, 2);
            });

            // now apply for the ghost nodes
            ExpressionIOUtils::EvaluateExpressionOnGhostNodes(
                r_communicator,
                rExpression,
                [](auto& rNode, const auto& rValues){
                    rNode.X0() = rValues[0];
                    rNode.Y0() = rValues[1];
                    rNode.Z0() = rValues[2];
                });

            break;
        case Configuration::Current:
            IndexPartition<IndexType>(number_of_nodes).for_each([&r_mesh, &rExpression](const IndexType Index) {
                const IndexType entity_data_start_index = Index * 3;
                auto& r_node = *(r_mesh.NodesBegin() + Index);
                r_node.X() = rExpression.Evaluate(Index, entity_data_start_index, 0);
                r_node.Y() = rExpression.Evaluate(Index, entity_data_start_index, 1);
                r_node.Z() = rExpression.Evaluate(Index, entity_data_start_index, 2);
            });

            // now apply for the ghost nodes
            ExpressionIOUtils::EvaluateExpressionOnGhostNodes(
                r_communicator,
                rExpression,
                [](auto& rNode, const auto& rValues){
                    rNode.X() = rValues[0];
                    rNode.Y() = rValues[1];
                    rNode.Z() = rValues[2];
                });

            break;
    }

    KRATOS_CATCH("");
}

template<MeshType TMeshType>
void NodalPositionExpressionIO::Read(
    ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
    const Configuration& rConfiguration)
{
    auto p_expression = Input(rContainerExpression.GetModelPart(),
                                                     rConfiguration, TMeshType)
                            .Execute();

    rContainerExpression.SetExpression(p_expression);
}

template<MeshType TMeshType>
void NodalPositionExpressionIO::Write(
    const ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainerExpression,
    const Configuration& rConfiguration)
{
    Output(*rContainerExpression.pGetModelPart(),
                                  rConfiguration, TMeshType)
        .Execute(rContainerExpression.GetExpression());
}

#define KRATOS_INSTANTIATE_CONTAINER_NODAL_POSITION_EXPRESSION_IO(MESH_TYPE)                                                                                                            \
    template KRATOS_API(KRATOS_CORE) void NodalPositionExpressionIO::Read(ContainerExpression<ModelPart::NodesContainerType, MESH_TYPE>&, const Configuration&);        \
    template KRATOS_API(KRATOS_CORE) void NodalPositionExpressionIO::Write(const ContainerExpression<ModelPart::NodesContainerType, MESH_TYPE>&, const Configuration&); \

KRATOS_INSTANTIATE_CONTAINER_NODAL_POSITION_EXPRESSION_IO(MeshType::Local)
KRATOS_INSTANTIATE_CONTAINER_NODAL_POSITION_EXPRESSION_IO(MeshType::Interface)
KRATOS_INSTANTIATE_CONTAINER_NODAL_POSITION_EXPRESSION_IO(MeshType::Ghost)

#undef KRATOS_INSTANTIATE_CONTAINER_NODAL_POSITION_EXPRESSION_IO

} // namespace Kratos
