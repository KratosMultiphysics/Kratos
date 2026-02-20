//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Suneth Warnakulasuriya
//

// System includes
#include <numeric>

// External includes
#include "nanoflann.hpp"

// Project includes
#include "includes/model_part.h"
#include "expression/literal_flat_expression.h"
#include "utilities/model_part_utils.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "integrated_neareset_entity_explicit_damping.h"

namespace Kratos {

template <class TContainerType>
IntegratedNearestEntityExplicitDamping<TContainerType>::IntegratedNearestEntityExplicitDamping(
    Model& rModel,
    Parameters Settings,
    const IndexType Stride)
    : mStride(Stride)
{
    KRATOS_TRY

    Parameters defaults = Parameters(R"(
    {
        "damping_type"               : "integrated_nearest_entity",
        "damping_function_type"      : "cosine",
        "damping_distance_multiplier": 100.0,
        "damping_max_items_in_bucket": 10,
        "damped_model_part_settings" : {}
    })" );

    Settings.ValidateAndAssignDefaults(defaults);

    mLeafMaxSize = Settings["damping_max_items_in_bucket"].GetInt();

    mpKernelFunction = Kratos::make_unique<DampingFunction>(Settings["damping_function_type"].GetString(), Settings["damping_distance_multiplier"].GetDouble());
    mComponentWiseDampedModelParts = OptimizationUtils::GetComponentWiseModelParts(rModel, Settings["damped_model_part_settings"]);

    KRATOS_ERROR_IF_NOT(mStride == mComponentWiseDampedModelParts.size() || mComponentWiseDampedModelParts.size() == 0)
        << "Number of components found in the damped model parts mismatch with stride [ Stride = " << mStride
        << ", number of components in damped model parts = " << mComponentWiseDampedModelParts.size() << " ]. Settings = "
        << Settings;

    if (mComponentWiseDampedModelParts.empty()) {
        mComponentWiseDampedModelParts.resize(mStride);
    }

    // initialize the KDTrees with nullptr
    mComponentWiseKDTrees.resize(mStride, nullptr);
    mComponentWisePositionAdapters.resize(mStride);

    KRATOS_CATCH("");
}

template <class TContainerType>
void IntegratedNearestEntityExplicitDamping<TContainerType>::SetRadius(const ContainerExpression<TContainerType>& rDampingRadiusExpression)
{
    mpDampingRadius = rDampingRadiusExpression.Clone();
}

template <class TContainerType>
typename ContainerExpression<TContainerType>::Pointer IntegratedNearestEntityExplicitDamping<TContainerType>::GetRadius() const
{
    return mpDampingRadius;
}

template <class TContainerType>
IndexType IntegratedNearestEntityExplicitDamping<TContainerType>::GetStride() const
{
    return mStride;
}

template <class TContainerType>
std::vector<std::vector<ModelPart*>> IntegratedNearestEntityExplicitDamping<TContainerType>::GetDampedModelParts() const
{
    return mComponentWiseDampedModelParts;
}

template <class TContainerType>
void IntegratedNearestEntityExplicitDamping<TContainerType>::Update()
{
    KRATOS_TRY

    const auto stride = this->GetStride();

    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        auto& r_damped_model_parts = mComponentWiseDampedModelParts[i_comp];

        mComponentWisePositionAdapters[i_comp] = Kratos::make_shared<PositionAdapter>(r_damped_model_parts);

        if (!mComponentWisePositionAdapters[i_comp]->kdtree_get_point_count()) {
            // now construct the kd tree
            mComponentWiseKDTrees[i_comp] =  Kratos::make_shared<KDTreeIndexType>(
                3, *mComponentWisePositionAdapters[i_comp],
                nanoflann::KDTreeSingleIndexAdaptorParams(mLeafMaxSize));
            mComponentWiseKDTrees[i_comp]->buildIndex();
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void IntegratedNearestEntityExplicitDamping<TContainerType>::IntegratedNearestEntityExplicitDamping::Apply(
    std::vector<std::vector<double>>& rDampedWeights,
    const std::vector<double>& rWeights,
    const IndexType Index,
    const IndexType NumberOfNeighbours,
    const EntityPointVector& rNeighbours) const
{
    KRATOS_TRY

    const IndexType stride = this->GetStride();

    rDampedWeights.resize(stride, std::vector<double>(rNeighbours.size()));

    const auto radius = mpDampingRadius->GetExpression().Evaluate(Index, Index, 0);

    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        auto& r_damped_weights = rDampedWeights[i_comp];
        r_damped_weights.resize(rNeighbours.size());
        if (mComponentWiseKDTrees[i_comp].get()) {
            for (IndexType i_neighbour = 0; i_neighbour < NumberOfNeighbours; ++i_neighbour) {

                double squared_distance;
                unsigned int global_index;
                mComponentWiseKDTrees[i_comp]->knnSearch(OptimizationUtils::GetEntityPosition(*rNeighbours[i_neighbour]).data().begin(), 1, &global_index, &squared_distance);

                r_damped_weights[i_neighbour] = rWeights[i_neighbour] * mpKernelFunction->ComputeWeight(radius, std::sqrt(squared_distance));
            }
        } else {
            std::copy(rWeights.begin(), rWeights.begin() + NumberOfNeighbours, r_damped_weights.begin());
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void IntegratedNearestEntityExplicitDamping<TContainerType>::IntegratedNearestEntityExplicitDamping::CalculateMatrix(
    Matrix& rOutput,
    const IndexType ComponentIndex) const
{
    KRATOS_TRY

    const auto stride = this->GetStride();
    const auto& r_radius_exp = mpDampingRadius->GetExpression();
    const auto& r_container = mpDampingRadius->GetContainer();
    const auto number_of_entities = r_container.size();

    KRATOS_ERROR_IF_NOT(ComponentIndex < stride)
        << "Invalid component index [ component index = " << ComponentIndex
        << ", stride = " << stride << " ].\n";

    if (rOutput.size1() != number_of_entities || rOutput.size2() != number_of_entities) {
        rOutput.resize(number_of_entities, number_of_entities, false);
    }

    rOutput.clear();

    if (!mComponentWiseKDTrees[ComponentIndex].get()) {
        // kd tree is not constructed. then there are no damped model parts
        // hence returning a Identity matrix.
        IndexPartition<IndexType>(number_of_entities).for_each([&rOutput, number_of_entities](const auto Index) {
            *(rOutput.data().begin() + Index * number_of_entities + Index) = 1.0;
        });
    } else {

        using search_adapter_type = NanoFlannSingleContainerPositionAdapter<TContainerType>;

        using search_distance_metric_type = typename nanoflann::metric_L2_Simple::traits<double, search_adapter_type>::distance_t;

        using search_tree_type = nanoflann::KDTreeSingleIndexAdaptor<search_distance_metric_type, search_adapter_type, 3>;

        // create the adapter
        search_adapter_type adapter(&r_container);

        // create domain search tree
        search_tree_type search_tree(3, adapter, nanoflann::KDTreeSingleIndexAdaptorParams(mLeafMaxSize));
        search_tree.buildIndex();

        auto& comp_search_tree = *mComponentWiseKDTrees[ComponentIndex];
        auto& kernel_function = *mpKernelFunction;

        IndexPartition<IndexType>(number_of_entities).for_each(NanoFlannKDTreeThreadLocalStorage<typename search_adapter_type::PointerVectorType>(), [&rOutput, &r_container, &r_radius_exp, &search_tree, &comp_search_tree, &kernel_function, &adapter, number_of_entities](const auto Index, auto& rTLS) {
            const auto radius = r_radius_exp.Evaluate(Index, Index, 0);

            // search for entities within radius
            search_tree.radiusSearch(OptimizationUtils::GetEntityPosition(*(r_container.begin() + Index)).data().begin(), radius * radius, rTLS.mNeighbourIndicesAndSquaredDistances, nanoflann::SearchParams());

            // update the neighbour entities from found indices
            const IndexType number_of_neighbors = rTLS.mNeighbourIndicesAndSquaredDistances.size();
            adapter.GetResultingEntityPointersVector(rTLS.mNeighbourEntityPoints, rTLS.mNeighbourIndicesAndSquaredDistances);

            for (IndexType i_neighbour = 0; i_neighbour < number_of_neighbors; ++i_neighbour) {
                const auto& r_neighbour_point = *rTLS.mNeighbourEntityPoints[i_neighbour];

                double squared_distance;
                unsigned int global_index;
                comp_search_tree.knnSearch(OptimizationUtils::GetEntityPosition(r_neighbour_point).data().begin(), 1, &global_index, &squared_distance);

                const double weight = kernel_function.ComputeWeight(radius, std::sqrt(squared_distance));

                *(rOutput.data().begin() + Index * number_of_entities + rTLS.mNeighbourIndicesAndSquaredDistances[i_neighbour].first) = weight;
            }

        });
    }

    KRATOS_CATCH("");
}

// template instantiations
template class IntegratedNearestEntityExplicitDamping<ModelPart::NodesContainerType>;
template class IntegratedNearestEntityExplicitDamping<ModelPart::ConditionsContainerType>;
template class IntegratedNearestEntityExplicitDamping<ModelPart::ElementsContainerType>;

} // namespace Kratos.