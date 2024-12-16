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

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "expression/c_array_expression_io.h"
#include "expression/variable_expression_io.h"

// Application includes
#include "properties_variable_expression_io.h"

// Include base h
#include "collective_expression_io.h"

namespace Kratos {

template<class TRawDataType>
void CollectiveExpressionIO::Read(
    CollectiveExpression& rCollectiveExpression,
    TRawDataType const* pBegin,
    int const* NumberOfEntities,
    int const** pListShapeBegin,
    int const* ShapeSizes,
    const int NumberOfContainers)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(NumberOfContainers > 0 && static_cast<IndexType>(NumberOfContainers) == rCollectiveExpression.GetContainerExpressions().size())
        << "Number of containers mismatch. [ Input number of containers = " << NumberOfContainers
        << ", CollectiveExpression number of containers = "
        << rCollectiveExpression.GetContainerExpressions().size() << " ].\n";

    for (auto& p_container_variable_data_holder : rCollectiveExpression.GetContainerExpressions()) {
        std::visit([&pBegin, &pListShapeBegin, &ShapeSizes, &NumberOfEntities](auto& v) {
            CArrayExpressionIO::Read(*v, pBegin, *NumberOfEntities, *pListShapeBegin, *ShapeSizes);

            // now offset everything
            pBegin += v->GetContainer().size() * v->GetItemComponentCount();
            ++pListShapeBegin;
            ++ShapeSizes;
            ++NumberOfEntities;
        }, p_container_variable_data_holder);
    }

    KRATOS_CATCH("");
}

template<class TRawDataType>
void CollectiveExpressionIO::Move(
    CollectiveExpression& rCollectiveExpression,
    TRawDataType* pBegin,
    int const* NumberOfEntities,
    int const** pListShapeBegin,
    int const* ShapeSizes,
    const int NumberOfContainers)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(NumberOfContainers > 0 && static_cast<IndexType>(NumberOfContainers) == rCollectiveExpression.GetContainerExpressions().size())
        << "Number of containers mismatch. [ Input number of containers = " << NumberOfContainers
        << ", CollectiveExpression number of containers = "
        << rCollectiveExpression.GetContainerExpressions().size() << " ].\n";

    for (auto& p_container_variable_data_holder : rCollectiveExpression.GetContainerExpressions()) {
        std::visit([&pBegin, &pListShapeBegin, &ShapeSizes, &NumberOfEntities](const auto& v) {
            CArrayExpressionIO::Move(*v, pBegin, *NumberOfEntities, *pListShapeBegin, *ShapeSizes);

            // now offset everything
            pBegin += v->GetContainer().size() * v->GetItemComponentCount();
            ++pListShapeBegin;
            ++ShapeSizes;
            ++NumberOfEntities;
        }, p_container_variable_data_holder);
    }

    KRATOS_CATCH("");
}

template<class TRawDataType>
void CollectiveExpressionIO::Write(
    const CollectiveExpression& rCollectiveExpression,
    TRawDataType* pBegin,
    const int Size)
{
    KRATOS_ERROR_IF_NOT(Size > 0 && static_cast<IndexType>(Size) == rCollectiveExpression.GetCollectiveFlattenedDataSize())
        << "The size of the double vector does not match with the required "
           "collective expression size. [ "
           "Size = "
        << Size << ", collective expression data size = "
        << rCollectiveExpression.GetCollectiveFlattenedDataSize() << " ].\n";

    for (const auto& p_container_variable_data_holder : rCollectiveExpression.GetContainerExpressions()) {
        std::visit([&pBegin](const auto& v) {
            // get the shape of the container expression.
            const auto& r_shape = v->GetItemShape();

            // transform unsigned Index type shape to signed int shape.
            std::vector<int> shape(r_shape.size());
            std::transform(r_shape.begin(), r_shape.end(), shape.begin(), [](const IndexType Value) -> int { return Value; });

            // get the number of entities in the container.
            const auto number_of_entities = v->GetContainer().size();

            // evaluate the expression and put the result in a continuous array starting with pBegin.
            CArrayExpressionIO::Write(*v, pBegin, number_of_entities * v->GetItemComponentCount());

            // increase the offset to place the evaluated values of the next container expression correctly.
            pBegin += number_of_entities * v->GetItemComponentCount();
        }, p_container_variable_data_holder);
    }
}

void CollectiveExpressionIO::Read(
    CollectiveExpression& rCollectiveExpression,
    const std::vector<ContainerVariableType>& rContainerVariables)
{
    const auto& r_container_expressions = rCollectiveExpression.GetContainerExpressions();

    KRATOS_ERROR_IF_NOT(r_container_expressions.size() == rContainerVariables.size())
        << "Container expressions size and variables size mismatch. [ number of container expressions = "
        << r_container_expressions.size() << ", number of variable container types = " << rContainerVariables.size() << " ].\n";

    for (IndexType i = 0; i < r_container_expressions.size(); ++i) {
        std::visit([](auto& pContainer, auto& pContainerVariable) {
            using container_expression_type = std::decay_t<decltype(*pContainer)>;
            using container_variable_type = std::decay_t<decltype(*pContainerVariable)>;

            if constexpr(std::is_same_v<container_expression_type, ContainerExpression<ModelPart::NodesContainerType, MeshType::Local>> ||
                         std::is_same_v<container_expression_type, ContainerExpression<ModelPart::NodesContainerType, MeshType::Ghost>> ||
                         std::is_same_v<container_expression_type, ContainerExpression<ModelPart::NodesContainerType, MeshType::Interface>>) {
                if constexpr(std::is_same_v<container_variable_type, HistoricalVariable>) {
                    VariableExpressionIO::Read(*pContainer, pContainerVariable->mVariable, true);
                } else if constexpr(std::is_same_v<container_variable_type, NonHistoricalVariable>) {
                    VariableExpressionIO::Read(*pContainer, pContainerVariable->mVariable, false);
                } else {
                    KRATOS_ERROR
                        << "Nodal expressions only supports HistoricalVariable "
                           "and NonHistoricalVariable container types.\n";
                }
            } else {
                if constexpr(std::is_same_v<container_variable_type, NonHistoricalVariable>) {
                    VariableExpressionIO::Read(*pContainer, pContainerVariable->mVariable);
                } else if constexpr(std::is_same_v<container_variable_type, PropertiesVariable>) {
                    PropertiesVariableExpressionIO::Read(*pContainer, pContainerVariable->mVariable);
                } else {
                    KRATOS_ERROR << "Element/Condition expressions only "
                                    "supports NonHistoricalVariable and "
                                    "PropertiesVariable container types.\n";
                }
            }
        }, r_container_expressions[i], rContainerVariables[i]);
    }
}

void CollectiveExpressionIO::Read(
    CollectiveExpression& rCollectiveExpression,
    const ContainerVariableType& rContainerVariable)
{
    std::vector<ContainerVariableType> variables(rCollectiveExpression.GetContainerExpressions().size(), rContainerVariable);
    Read(rCollectiveExpression, variables);
}

void CollectiveExpressionIO::Write(
    const CollectiveExpression& rCollectiveExpression,
    const std::vector<ContainerVariableType>& rContainerVariables)
{
    const auto& r_container_expressions = rCollectiveExpression.GetContainerExpressions();

    KRATOS_ERROR_IF_NOT(r_container_expressions.size() == rContainerVariables.size())
        << "Container expressions size and variables size mismatch. [ number of container expressions = "
        << r_container_expressions.size() << ", number of variable container types = " << rContainerVariables.size() << " ].\n";

    for (IndexType i = 0; i < r_container_expressions.size(); ++i) {
        std::visit([](auto& pContainer, auto& pContainerVariable) {
            using container_expression_type = std::decay_t<decltype(*pContainer)>;
            using container_variable_type = std::decay_t<decltype(*pContainerVariable)>;

            if constexpr(std::is_same_v<container_expression_type, ContainerExpression<ModelPart::NodesContainerType, MeshType::Local>> ||
                         std::is_same_v<container_expression_type, ContainerExpression<ModelPart::NodesContainerType, MeshType::Ghost>> ||
                         std::is_same_v<container_expression_type, ContainerExpression<ModelPart::NodesContainerType, MeshType::Interface>>) {
                if constexpr(std::is_same_v<container_variable_type, HistoricalVariable>) {
                    VariableExpressionIO::Write(*pContainer, pContainerVariable->mVariable, true);
                } else if constexpr(std::is_same_v<container_variable_type, NonHistoricalVariable>) {
                    VariableExpressionIO::Write(*pContainer, pContainerVariable->mVariable, false);
                } else {
                    KRATOS_ERROR
                        << "Nodal expressions only supports HistoricalVariable "
                           "and NonHistoricalVariable container types.\n";
                }
            } else {
                if constexpr(std::is_same_v<container_variable_type, NonHistoricalVariable>) {
                    VariableExpressionIO::Write(*pContainer, pContainerVariable->mVariable);
                } else if constexpr(std::is_same_v<container_variable_type, PropertiesVariable>) {
                    PropertiesVariableExpressionIO::Write(*pContainer, pContainerVariable->mVariable);
                } else {
                    KRATOS_ERROR << "Element/Condition expressions only "
                                    "supports NonHistoricalVariable and "
                                    "PropertiesVariable container types.\n";
                }
            }
        }, r_container_expressions[i], rContainerVariables[i]);
    }
}

void CollectiveExpressionIO::Write(
    const CollectiveExpression& rCollectiveExpression,
    const ContainerVariableType& rContainerVariable)
{
    std::vector<ContainerVariableType> variables(rCollectiveExpression.GetContainerExpressions().size(), rContainerVariable);
    Write(rCollectiveExpression, variables);
}

// template instantiations
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::Read(CollectiveExpression&, int const*, int const*, int const**, int const*, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::Read(CollectiveExpression&, double const*, int const*, int const**, int const*, const int);

template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::Move(CollectiveExpression&, int*, int const*, int const**, int const*, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::Move(CollectiveExpression&, double*, int const*, int const**, int const*, const int);

template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::Write(const CollectiveExpression&, int*, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::Write(const CollectiveExpression&, double*, const int);

} // namespace Kratos
