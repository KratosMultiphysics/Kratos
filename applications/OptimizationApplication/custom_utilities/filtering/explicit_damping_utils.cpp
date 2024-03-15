//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <numeric>

// Project includes
#include "utilities/model_part_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes

// Include base h
#include "explicit_damping_utils.h"

namespace Kratos {

template<class TContainerType>
void ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntityForComponent(
    LiteralFlatExpression<double>& rOutputExpression,
    const ContainerExpression<TContainerType>& rDampingRadiusExpression,
    const std::vector<const ModelPart*>& rDampedModelParts,
    const FilterFunction& rFilterFunction,
    const IndexType BucketSize,
    const IndexType DampedComponentIndex)
{
    KRATOS_TRY

    const auto& r_filtering_model_part = rDampingRadiusExpression.GetModelPart();

    // first construct the damped entities
    EntityPointsVector<TContainerType> damped_entities;
    damped_entities.resize(std::accumulate(
        rDampedModelParts.begin(), rDampedModelParts.end(), 0UL,
        [](const IndexType Value, const auto& pModelPart) {
            return Value + ModelPartUtils::GetContainer<TContainerType>(*pModelPart).size();
        }));

    IndexType local_start = 0;
    for (const auto& p_model_part : rDampedModelParts) {
        const auto& r_container = ModelPartUtils::GetContainer<TContainerType>(*p_model_part);
        IndexPartition<IndexType>(r_container.size()).for_each([&damped_entities, &r_container, local_start](const auto Index) {
            *(damped_entities.begin() + local_start + Index) =  Kratos::make_shared<EntityPointType<TContainerType>>(*(r_container.begin() + Index), Index);
        });
        local_start += r_container.size();
    }

    const auto stride = rOutputExpression.GetItemComponentCount();
    const auto& r_container = ModelPartUtils::GetContainer<TContainerType>(r_filtering_model_part);

    if (damped_entities.empty()) {
        IndexPartition<IndexType>(r_container.size()).for_each([&rOutputExpression, stride, DampedComponentIndex](const auto Index) {
            *(rOutputExpression.begin() + Index * stride + DampedComponentIndex) = 1.0;
        });
        return;
    }

    // now construct the kd tree
    auto p_search_tree =  Kratos::make_shared<KDTree<TContainerType>>(damped_entities.begin(), damped_entities.end(), BucketSize);

    // now calculate the damping for each entity
    IndexPartition<IndexType>(rOutputExpression.NumberOfEntities()).for_each([&rOutputExpression, &r_container, &p_search_tree, &rFilterFunction, &rDampingRadiusExpression, stride, DampedComponentIndex](const auto Index){
        EntityPointType<TContainerType> entity_point(*(r_container.begin() + Index), Index);
        const auto data_begin_index = Index * stride;
        const auto radius = rDampingRadiusExpression.GetExpression().Evaluate(Index, Index, 0);

        double dummy_distance{};
        auto p_nearest_damped_entity_point = p_search_tree->SearchNearestPoint(entity_point, dummy_distance);

        double& value = *(rOutputExpression.begin() + data_begin_index + DampedComponentIndex);
        value = 1.0 - rFilterFunction.ComputeWeight(entity_point.Coordinates(), p_nearest_damped_entity_point->Coordinates(), radius);
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity(
    const ContainerExpression<TContainerType>& rDampingRadiusExpression,
    const std::vector<std::vector<const ModelPart*>>& rDampedModelParts,
    const std::vector<IndexType>& rShape,
    const std::string& rFilterFunctionType,
    const IndexType BucketSize)
{
    KRATOS_TRY

    const auto number_of_entities = rDampingRadiusExpression.GetContainer().size();

    // now create the expression
    auto p_expression = LiteralFlatExpression<double>::Create(number_of_entities, rShape);
    auto& r_expression = *p_expression;

    const auto stride = r_expression.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == rDampedModelParts.size())
        << "Number of components in the shape and the number of model parts list mismatch. [ "
        << "number of components in shape = " << stride << ", number of model part lists = "
        << rDampedModelParts.size() << " ].\n";

    FilterFunction filter_function(rFilterFunctionType);

    // compute damping coefficients
    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        ComputeDampingCoefficientsBasedOnNearestEntityForComponent(
            r_expression, rDampingRadiusExpression, rDampedModelParts[i_comp],
            filter_function, BucketSize, i_comp);
    }

    auto result = rDampingRadiusExpression;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

// template instantiations

#define KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ContainerType)                                                                                  \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ContainerType> ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity( \
        const ContainerExpression<ContainerType>&,                                                                                                         \
        const std::vector<std::vector<const ModelPart*>>&,                                                                                                 \
        const std::vector<std::size_t>&, const std::string&, const std::size_t);                                                                           \
    template KRATOS_API(OPTIMIZATION_APPLICATION) void ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntityForComponent(                   \
        LiteralFlatExpression<double>&, const ContainerExpression<ContainerType>&,                                                                         \
        const std::vector<const ModelPart*>&, const FilterFunction&,                                                                                       \
        const IndexType, const IndexType);

KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType);
KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType);
KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType);

#undef KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE

} // namespace Kratos