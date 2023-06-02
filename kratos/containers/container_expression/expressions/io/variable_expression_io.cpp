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
#include "containers/container_expression/variable_expression_data_io.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Include base h
#include "variable_expression_io.h"

namespace Kratos {

namespace VariableExpressionIOHelperUtilities
{

template <class TContainerType, class TContainerDataIO>
Expression::Pointer Read(
    const TContainerType& rContainer,
    const VariableExpressionIO::VariableType pVariable)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    if (number_of_entities != 0) {
        // initialize the shape with the first entity value
        return std::visit([&rContainer, number_of_entities](auto pVariable) {
            using data_type = typename std::remove_const_t<std::remove_pointer_t<decltype(pVariable)>>::Type;

            using raw_data_type = std::conditional_t<std::is_same_v<data_type, int>, int, double>;

            VariableExpressionDataIO<data_type> variable_flatten_data_io(TContainerDataIO::GetValue(*rContainer.begin(), *pVariable));

            auto p_expression = LiteralFlatExpression<raw_data_type>::Create(number_of_entities, variable_flatten_data_io.GetItemShape());
            auto& r_expression = *p_expression;

            IndexPartition<IndexType>(number_of_entities).for_each([&rContainer, &pVariable, &variable_flatten_data_io, &r_expression](const IndexType Index) {
                    const auto& values = TContainerDataIO::GetValue(*(rContainer.begin() + Index), *pVariable);
                    variable_flatten_data_io.Read(r_expression, Index, values);
                });

            return Kratos::intrusive_ptr<Expression>(&r_expression);

        }, pVariable);
    }

    return nullptr;

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO>
void Write(
    TContainerType& rContainer,
    Communicator& rCommunicator,
    const Expression& rExpression,
    const VariableExpressionIO::VariableType pVariable)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    if (number_of_entities > 0) {
        std::visit([&rContainer, number_of_entities, &rExpression, &rCommunicator](auto pVariable) {
            using data_type = typename std::remove_const_t<std::remove_pointer_t<decltype(pVariable)>>::Type;

            VariableExpressionDataIO<data_type> variable_flatten_data_io(rExpression.GetItemShape());

            // initialize the container variables first
            if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
                // initializes ghost nodes as for the later synchronization
                // only, the nodal non historical values needs to be set unless
                // they are properly initialized. Otherwise, in synchronization, the variables will
                // not be there in the ghost nodes hence seg faults.

                // the vectors and matrices needs to be initialized in historical and non-historical
                // data containers because they need to be initialized with the correct size for synchronization
                if constexpr(std::is_same_v<data_type, Vector> || std::is_same_v<data_type, Matrix>) {
                    data_type dummy_value{};
                    variable_flatten_data_io.Assign(dummy_value, rExpression, 0);
                    if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
                        VariableUtils().SetVariable(*pVariable, dummy_value, rCommunicator.GhostMesh().Nodes());
                    } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                        VariableUtils().SetNonHistoricalVariable(*pVariable, dummy_value, rCommunicator.GhostMesh().Nodes());
                    }
                } else {
                    // if it is a static type, then it only needs to be initialized in the non-historical container with zeros.
                    // historical container should be initialized with the default values when the container is created.
                    if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                        VariableUtils().SetNonHistoricalVariableToZero(*pVariable, rCommunicator.GhostMesh().Nodes());
                    }
                }
            }

            IndexPartition<IndexType>(number_of_entities).for_each(data_type{}, [&rContainer, &pVariable, &rExpression, &variable_flatten_data_io](const IndexType Index, data_type& rValue){
                variable_flatten_data_io.Assign(rValue, rExpression, Index);
                TContainerDataIO::SetValue(*(rContainer.begin() + Index), *pVariable, rValue);
            });

            if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
                // synchronize nodal values
                if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
                    rCommunicator.SynchronizeVariable(*pVariable);
                } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                    rCommunicator.SynchronizeNonHistoricalVariable(*pVariable);
                }
            }
        }, pVariable);
    }

    KRATOS_CATCH("");
}

} // namespace VariableExpressionIOHelperUtilities

template<class TDataType>
VariableExpressionIO::VariableExpressionIO(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const ContainerType& rContainerType,
    const MeshType& rMeshType)
    : mrModelPart(rModelPart),
      mpVariable(&rVariable),
      mContainerType(rContainerType),
      mMeshType(rMeshType)
{
}

template <class TDataType, class TMeshType>
VariableExpressionIO::VariableExpressionIO(
    ContainerExpression<ModelPart::NodesContainerType, TMeshType>& rContainer,
    const Variable<TDataType>& rVariable,
    const bool IsHistorical)
    : VariableExpressionIO(
        rContainer.GetModelPart(),
        rVariable,
        IsHistorical ? ContainerType::NodalHistorical : ContainerType::NodalNonHistorical,
        std::is_same_v<TMeshType, Kratos::MeshType::Local> ? MeshType::Local : std::is_same_v<TMeshType, Kratos::MeshType::Interface> ? MeshType::Interface : MeshType::Ghost)
{

}

template <class TContainerType, class TDataType, class TMeshType>
VariableExpressionIO::VariableExpressionIO(
    ContainerExpression<TContainerType, TMeshType>& rContainer,
    const Variable<TDataType>& rVariable)
    : VariableExpressionIO(
        rContainer.GetModelPart(),
        rVariable,
        std::is_same_v<TContainerType, ModelPart::ConditionsContainerType> ? ContainerType::ConditionNonHistorical : ContainerType::ElementNonHistorical,
        std::is_same_v<TMeshType, Kratos::MeshType::Local> ? MeshType::Local : std::is_same_v<TMeshType, Kratos::MeshType::Interface> ? MeshType::Interface : MeshType::Ghost)
{
}

Expression::Pointer VariableExpressionIO::Read()
{
    switch (mContainerType) {
        case ContainerType::NodalHistorical:
            return VariableExpressionIOHelperUtilities::Read<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>>(GetMesh().Nodes(), mpVariable);
        case ContainerType::NodalNonHistorical:
            return VariableExpressionIOHelperUtilities::Read<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>(GetMesh().Nodes(), mpVariable);
        case ContainerType::ConditionNonHistorical:
            return VariableExpressionIOHelperUtilities::Read<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>(GetMesh().Conditions(), mpVariable);
        case ContainerType::ElementNonHistorical:
            return VariableExpressionIOHelperUtilities::Read<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>(GetMesh().Elements(), mpVariable);
    }

    return nullptr;
}

void VariableExpressionIO::Write(const Expression& rExpression)
{
    auto& r_communicator = mrModelPart.GetCommunicator();

    switch (mContainerType) {
        case ContainerType::NodalHistorical:
            VariableExpressionIOHelperUtilities::Write<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>>(GetMesh().Nodes(), r_communicator, rExpression, mpVariable);
            break;
        case ContainerType::NodalNonHistorical:
            VariableExpressionIOHelperUtilities::Write<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>(GetMesh().Nodes(), r_communicator, rExpression, mpVariable);
            break;
        case ContainerType::ConditionNonHistorical:
            VariableExpressionIOHelperUtilities::Write<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>(GetMesh().Conditions(), r_communicator, rExpression, mpVariable);
            break;
        case ContainerType::ElementNonHistorical:
            VariableExpressionIOHelperUtilities::Write<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>(GetMesh().Elements(), r_communicator, rExpression, mpVariable);
            break;
    }
}

ModelPart::MeshType& VariableExpressionIO::GetMesh()
{
    switch (mMeshType) {
        case MeshType::Local:
            return mrModelPart.GetCommunicator().LocalMesh();
        case MeshType::Interface:
            return mrModelPart.GetCommunicator().InterfaceMesh();
        case MeshType::Ghost:
            return mrModelPart.GetCommunicator().GhostMesh();
    }

    return mrModelPart.GetCommunicator().LocalMesh();
}

// template instantiations
#define KRATOS_VARIABLE_EXPRESSION_IO_NODES(...)                                                                                                                                     \
    template VariableExpressionIO::VariableExpressionIO(ContainerExpression<ModelPart::NodesContainerType, Kratos::MeshType::Local>&, const Variable<__VA_ARGS__>&, const bool);     \
    template VariableExpressionIO::VariableExpressionIO(ContainerExpression<ModelPart::NodesContainerType, Kratos::MeshType::Interface>&, const Variable<__VA_ARGS__>&, const bool); \
    template VariableExpressionIO::VariableExpressionIO(ContainerExpression<ModelPart::NodesContainerType, Kratos::MeshType::Ghost>&, const Variable<__VA_ARGS__>&, const bool);

#define KRATOS_VARIABLE_EXPRESSION_IO_ENTITIES(ENTITY_TYPE, ...)                                                                                       \
    template VariableExpressionIO::VariableExpressionIO(ContainerExpression<ENTITY_TYPE, Kratos::MeshType::Local>&, const Variable<__VA_ARGS__>&);     \
    template VariableExpressionIO::VariableExpressionIO(ContainerExpression<ENTITY_TYPE, Kratos::MeshType::Interface>&, const Variable<__VA_ARGS__>&); \
    template VariableExpressionIO::VariableExpressionIO(ContainerExpression<ENTITY_TYPE, Kratos::MeshType::Ghost>&, const Variable<__VA_ARGS__>&);

#define KRATOS_VARIABLE_EXPRESSION_IO(...)                                                                                                                                              \
    template VariableExpressionIO::VariableExpressionIO(ModelPart&, const Variable<__VA_ARGS__>&, const VariableExpressionIO::ContainerType&, const VariableExpressionIO::MeshType&);   \
    KRATOS_VARIABLE_EXPRESSION_IO_NODES(__VA_ARGS__)                                                                                                                                    \
    KRATOS_VARIABLE_EXPRESSION_IO_ENTITIES(ModelPart::ConditionsContainerType, __VA_ARGS__)                                                                                             \
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
