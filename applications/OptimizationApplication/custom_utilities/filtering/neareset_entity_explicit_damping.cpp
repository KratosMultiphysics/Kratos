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
#include "includes/model_part.h"
#include "expression/literal_flat_expression.h"
#include "utilities/model_part_utils.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "neareset_entity_explicit_damping.h"

namespace Kratos {

template <class TContainerType>
NearestEntityExplicitDamping<TContainerType>::NearestEntityExplicitDamping(
    Model& rModel,
    Parameters Settings,
    const IndexType Stride)
    : mStride(Stride)
{
    KRATOS_TRY

    Parameters defaults = Parameters(R"(
    {
        "damping_type"              : "nearest_entity",
        "damping_function_type"     : "cosine",
        "damped_model_part_settings": {}
    })" );

    Settings.ValidateAndAssignDefaults(defaults);

    mpKernelFunction = Kratos::make_unique<FilterFunction>(Settings["damping_function_type"].GetString());
    mComponentWiseDampedModelParts = OptimizationUtils::GetComponentWiseModelParts(rModel, Settings["damped_model_part_settings"]);

    KRATOS_ERROR_IF_NOT(mStride == mComponentWiseDampedModelParts.size() || mComponentWiseDampedModelParts.size() == 0)
        << "Number of components found in the damped model parts mismatch with stride [ Stride = " << mStride
        << ", number of components in damped model parts = " << mComponentWiseDampedModelParts.size() << " ]. Settings = "
        << Settings;

    if (mComponentWiseDampedModelParts.empty()) {
        mComponentWiseDampedModelParts.resize(mStride);
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void NearestEntityExplicitDamping<TContainerType>::SetRadius(const ContainerExpression<TContainerType>& rDampingRadiusExpression)
{
    mpDampingRadius = rDampingRadiusExpression.Clone();
}

template <class TContainerType>
typename ContainerExpression<TContainerType>::Pointer NearestEntityExplicitDamping<TContainerType>::GetRadius() const
{
    return mpDampingRadius;
}

template <class TContainerType>
IndexType NearestEntityExplicitDamping<TContainerType>::GetStride() const
{
    return mStride;
}

template <class TContainerType>
std::vector<std::vector<ModelPart*>> NearestEntityExplicitDamping<TContainerType>::GetDampedModelParts() const
{
    return mComponentWiseDampedModelParts;
}

template <class TContainerType>
void NearestEntityExplicitDamping<TContainerType>::Update()
{
    KRATOS_TRY

    const auto stride = this->GetStride();
    auto& r_container = mpDampingRadius->GetContainer();
    const auto& r_damping_radius = mpDampingRadius->GetExpression();

    // first create the output expression with correct required shape
    LiteralFlatExpression<double>::Pointer p_output_expression;
    if (stride == 1) {
        p_output_expression = LiteralFlatExpression<double>::Create(r_container.size(), {});
    } else {
        p_output_expression = LiteralFlatExpression<double>::Create(r_container.size(), {stride});
    }

    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        auto& r_damped_model_parts = mComponentWiseDampedModelParts[i_comp];

        PositionAdapter adapter(r_damped_model_parts);

        if (adapter.kdtree_get_point_count() == 0) {
            IndexPartition<IndexType>(r_container.size()).for_each([&p_output_expression, stride, i_comp](const auto Index) {
                *(p_output_expression->begin() + Index * stride + i_comp) = 1.0;
            });
        } else {
            // now construct the kd tree
            KDTreeIndexType kd_tree_index(
                3, adapter,
                nanoflann::KDTreeSingleIndexAdaptorParams(
                    mLeafMaxSize, nanoflann::KDTreeSingleIndexAdaptorFlags::None, 0));
            kd_tree_index.buildIndex();

            const auto& kernel_function = *mpKernelFunction;

            // now calculate the damping for each entity
            IndexPartition<IndexType>(r_container.size()).for_each([&p_output_expression, &r_container, &kd_tree_index, &kernel_function, &r_damping_radius, stride, i_comp](const auto Index){
                const auto data_begin_index = Index * stride;
                const auto radius = r_damping_radius.Evaluate(Index, Index, 0);

                double squared_distance;
                unsigned int global_index;
                kd_tree_index.knnSearch(OptimizationUtils::GetEntityPosition(*(r_container.begin() + Index)).data().begin(), 1, &global_index, &squared_distance);

                double& value = *(p_output_expression->begin() + data_begin_index + i_comp);
                value = 1.0 - kernel_function.ComputeWeight(radius, std::sqrt(squared_distance));
            });
        }
    }

    mpDampingCoefficients = p_output_expression;

    KRATOS_CATCH("");
}

template <class TContainerType>
void NearestEntityExplicitDamping<TContainerType>::NearestEntityExplicitDamping::Apply(
    std::vector<std::vector<double>>& rDampedWeights,
    const std::vector<double>& rWeights,
    const IndexType Index,
    const IndexType NumberOfNeighbours,
    const EntityPointVector& rNeighbours) const
{
    KRATOS_TRY

    const IndexType stride = this->GetStride();

    rDampedWeights.resize(stride, std::vector<double>(rNeighbours.size()));

    const auto data_begin_index = Index * stride;

    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        auto& r_damped_weights = rDampedWeights[i_comp];
        r_damped_weights.resize(rNeighbours.size());
        for (IndexType i_neighbour = 0; i_neighbour < NumberOfNeighbours; ++i_neighbour) {
            r_damped_weights[i_neighbour] = rWeights[i_neighbour] * mpDampingCoefficients->Evaluate(Index, data_begin_index,  i_comp);
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void NearestEntityExplicitDamping<TContainerType>::NearestEntityExplicitDamping::CalculateMatrix(
    Matrix& rOutput,
    const IndexType ComponentIndex) const
{
    KRATOS_TRY

    const auto stride = this->GetStride();
    const auto& r_expression = *mpDampingCoefficients;
    const auto number_of_entities = r_expression.NumberOfEntities();

    KRATOS_ERROR_IF_NOT(ComponentIndex < stride)
        << "Invalid component index [ component index = " << ComponentIndex
        << ", stride = " << stride << " ].\n";

    if (rOutput.size1() != number_of_entities || rOutput.size2() != number_of_entities) {
        rOutput.resize(number_of_entities, number_of_entities, false);
    }

    rOutput.clear();

    IndexPartition<IndexType>(number_of_entities).for_each([&rOutput, &r_expression, stride, ComponentIndex, number_of_entities](const auto Index) {
        *(rOutput.data().begin() + Index * number_of_entities + Index) = r_expression.Evaluate(Index, Index * stride, ComponentIndex);
    });

    KRATOS_CATCH("");
}

// template instantiations
template class NearestEntityExplicitDamping<ModelPart::NodesContainerType>;
template class NearestEntityExplicitDamping<ModelPart::ConditionsContainerType>;
template class NearestEntityExplicitDamping<ModelPart::ElementsContainerType>;

} // namespace Kratos.
