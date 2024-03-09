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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "explicit_damping_utils.h"

namespace Kratos {

std::vector<std::vector<const ModelPart*>> ExplicitDampingUtils::GetComponentWiseDampedModelParts(
    Model& rModel,
    Parameters Settings,
    const IndexType NumberOfComponents)
{
    KRATOS_TRY

    std::vector<std::vector<const ModelPart*>> result;
    result.resize(NumberOfComponents);

    for (auto it = Settings.begin(); it != Settings.end(); ++it) {
        IndexType i_component = 0;
        for (const auto& r_value : *it) {
            if (r_value.GetBool()) {
                result[i_component].push_back(&rModel.GetModelPart(it.name()));
            }
            ++i_component;
        }

        KRATOS_ERROR_IF_NOT(i_component == NumberOfComponents)
            << "Number of damping component mismatch [ Number of components required = "
            << NumberOfComponents << ", number of components specified for damping model part "
            << it.name() << " = " << i_component << " ].\n";
    }

    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
void ExplicitDampingUtils::FillEntityIdEntityPointMap(
    std::unordered_map<IndexType, typename EntityPointType<TContainerType>::Pointer>& rMap,
    const TContainerType& rContainer)
{
    rMap.clear();
    rMap.reserve(rContainer.size());
    for (IndexType i = 0; i < rContainer.size(); ++i) {
        auto p_entity = Kratos::make_shared<EntityPointType<TContainerType>>(*(rContainer.begin() + i), i);
        rMap[p_entity->GetEntity().Id()] = p_entity;
    };
}

template<class TContainerType>
void ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntityForComponent(
    LiteralFlatExpression<double>& rOutputExpression,
    const ContainerExpression<TContainerType>& rDampingRadiusExpression,
    const std::unordered_map<IndexType, typename EntityPointType<TContainerType>::Pointer>& rFilterEntityIdEntityPointMap,
    const std::vector<const ModelPart*>& rDampedModelParts,
    const FilterFunction& rFilterFunction,
    const IndexType BucketSize,
    const IndexType DampedComponentIndex)
{
    KRATOS_TRY

    const auto& r_filtering_model_part = rDampingRadiusExpression.GetModelPart();

    // first construct the damped entities
    EntityPointsVector<TContainerType> damped_entities;
    damped_entities.reserve(std::accumulate(
        rDampedModelParts.begin(), rDampedModelParts.end(), 0UL,
        [](const IndexType Value, const auto& pModelPart) {
            return Value + OptimizationUtils::GetContainer<TContainerType>(*pModelPart).size();
        }));
    for (const auto& p_model_part : rDampedModelParts) {
        const auto& r_entities = block_for_each<AccumReduction<typename EntityPointType<TContainerType>::Pointer>>(OptimizationUtils::GetContainer<TContainerType>(*p_model_part), [&rFilterEntityIdEntityPointMap, &p_model_part, &r_filtering_model_part](const auto& rEntity) -> typename EntityPointType<TContainerType>::Pointer {
            auto p_itr = rFilterEntityIdEntityPointMap.find(rEntity.Id());

            KRATOS_ERROR_IF(p_itr == rFilterEntityIdEntityPointMap.end())
                << "Entity with id = " << rEntity.Id() << " in damping model part "
                << p_model_part->FullName() << " not found in the filtering model part "
                << r_filtering_model_part.FullName() << ".\n";

            return p_itr->second;
        });

        damped_entities.insert(damped_entities.end(), r_entities.begin(), r_entities.end());
    }

    const auto stride = rOutputExpression.GetItemComponentCount();
    const auto& r_container = OptimizationUtils::GetContainer<TContainerType>(r_filtering_model_part);

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

    const auto& r_container = rDampingRadiusExpression.GetContainer();
    const auto number_of_entities = rDampingRadiusExpression.GetContainer().size();

    // now create the expression
    auto p_expression = LiteralFlatExpression<double>::Create(number_of_entities, rShape);
    auto& r_expression = *p_expression;

    const auto stride = r_expression.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == rDampedModelParts.size())
        << "Number of components in the shape and the number of model parts list mismatch. [ "
        << "number of components in shape = " << stride << ", number of model part lists = "
        << rDampedModelParts.size() << " ].\n";

    // create the map of damped ids and entity indices
    std::unordered_map<IndexType, typename EntityPointType<TContainerType>::Pointer> id_entity_point_map;
    FillEntityIdEntityPointMap(id_entity_point_map, r_container);

    FilterFunction filter_function(rFilterFunctionType);

    // compute damping coefficients
    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        ComputeDampingCoefficientsBasedOnNearestEntityForComponent(
            r_expression, rDampingRadiusExpression, id_entity_point_map,
            rDampedModelParts[i_comp], filter_function, BucketSize, i_comp);
    }

    auto result = rDampingRadiusExpression;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

// template instantiations

#define KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ContainerType)                                                                         \
    template KRATOS_API(OPTIMIZATION_APPLICATION) void ExplicitDampingUtils::FillEntityIdEntityPointMap(                                                   \
        std::unordered_map<IndexType, typename EntityPointType<ContainerType>::Pointer>&,                                                         \
        const ContainerType&);                                                                                                                    \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ContainerType> ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity( \
        const ContainerExpression<ContainerType>&,                                                                                                \
        const std::vector<std::vector<const ModelPart*>>&,                                                                                        \
        const std::vector<std::size_t>&, const std::string&, const std::size_t);                                                                  \
    template KRATOS_API(OPTIMIZATION_APPLICATION) void ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntityForComponent(                   \
        LiteralFlatExpression<double>&, const ContainerExpression<ContainerType>&,                                                                \
        const std::unordered_map<IndexType, typename EntityPointType<ContainerType>::Pointer>&,                                                   \
        const std::vector<const ModelPart*>&, const FilterFunction&,                                                                              \
        const IndexType, const IndexType);

KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType);
KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType);
KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType);

#undef KRATOS_INSTANTIATE_FILTER_UTILS_FOR_CONTAINER_TYPE

} // namespace Kratos