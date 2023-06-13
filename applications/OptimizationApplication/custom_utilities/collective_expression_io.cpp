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

// Include base h
#include "collective_expression_io.h"

namespace Kratos {

template<class TRawDataType>
void CollectiveExpressionIO::ReadCArray(
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
void CollectiveExpressionIO::MoveCArray(
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
void CollectiveExpressionIO::WriteCArray(
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

// template instantiations
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::ReadCArray(CollectiveExpression&, int const*, int const*, int const**, int const*, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::ReadCArray(CollectiveExpression&, double const*, int const*, int const**, int const*, const int);

template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::MoveCArray(CollectiveExpression&, int*, int const*, int const**, int const*, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::MoveCArray(CollectiveExpression&, double*, int const*, int const**, int const*, const int);

template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::WriteCArray(const CollectiveExpression&, int*, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) void CollectiveExpressionIO::WriteCArray(const CollectiveExpression&, double*, const int);

} // namespace Kratos
