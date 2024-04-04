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
        "damped_model_part_settings" : {}
    })" );

    Settings.ValidateAndAssignDefaults(defaults);

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

    mComponentWiseEntityPoints.resize(mStride);

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
        auto& r_entity_points = mComponentWiseEntityPoints[i_comp];

        r_entity_points.resize(std::accumulate(
            r_damped_model_parts.begin(), r_damped_model_parts.end(), 0UL,
            [](const IndexType Value, const auto& pModelPart) {
                return Value + ModelPartUtils::GetContainer<TContainerType>(*pModelPart).size();
            }));

        IndexType local_start = 0;
        for (const auto& p_model_part : r_damped_model_parts) {
            const auto& r_damping_container = ModelPartUtils::GetContainer<TContainerType>(*p_model_part);
            IndexPartition<IndexType>(r_damping_container.size()).for_each([&r_entity_points, &r_damping_container, local_start](const auto Index) {
                *(r_entity_points.begin() + local_start + Index) =  Kratos::make_shared<EntityPointType>(*(r_damping_container.begin() + Index), Index);
            });
            local_start += r_damping_container.size();
        }

        if (!r_entity_points.empty()) {
            // now construct the kd tree
            mComponentWiseKDTrees[i_comp] =  Kratos::make_shared<KDTree>(r_entity_points.begin(), r_entity_points.end(), mBucketSize);
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

    const auto radius = mpDampingRadius->GetExpression().Evaluate(Index, Index, 0);

    for (IndexType i_comp = 0; i_comp < this->GetStride(); ++i_comp) {
        auto& r_damped_weights = rDampedWeights[i_comp];

        if (mComponentWiseKDTrees[i_comp]) {
            for (IndexType i_neighbour = 0; i_neighbour < NumberOfNeighbours; ++i_neighbour) {

                double squared_distance;
                auto p_nearest_damped_entity_point = mComponentWiseKDTrees[i_comp]->SearchNearestPoint(*rNeighbours[i_neighbour], squared_distance);

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

    using tls = OptimizationUtils::KDTreeThreadLocalStorage<typename EntityPointType::Pointer>;

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

    if (mComponentWiseKDTrees[ComponentIndex].get() == nullptr) {
        // kd tree is not constructed. then there are no damped model parts
        // hence returning a Identity matrix.
        IndexPartition<IndexType>(number_of_entities).for_each([&rOutput, number_of_entities](const auto Index) {
            *(rOutput.data().begin() + Index * number_of_entities + Index) = 1.0;
        });
    } else {
        // now fill the points vector
        EntityPointVector domain_entities(r_container.size());
        IndexPartition<IndexType>(r_container.size()).for_each([&domain_entities, &r_container](const IndexType Index) {
            domain_entities[Index] = Kratos::make_shared<EntityPoint<EntityType>>(*(r_container.begin() + Index), Index);
        });

        // create domain search tree
        KDTree search_tree(domain_entities.begin(), domain_entities.end(), mBucketSize);
        auto& comp_search_tree = *mComponentWiseKDTrees[ComponentIndex];
        auto& kernel_function = *mpKernelFunction;

        IndexPartition<IndexType>(number_of_entities).for_each(tls(1000, 1), [&rOutput, &r_container, &r_radius_exp, &search_tree, &comp_search_tree, &kernel_function, number_of_entities](const auto Index, auto& rTLS) {
            const auto radius = r_radius_exp.Evaluate(Index, Index, 0);
            EntityPoint<EntityType> entity_point(*(r_container.begin() + Index), Index);
            const auto number_of_neighbors = search_tree.SearchInRadius(
                                                entity_point,
                                                radius,
                                                rTLS.mNeighbourEntityPoints.begin(),
                                                rTLS.mResultingSquaredDistances.begin(),
                                                1000);

            KRATOS_ERROR_IF(number_of_neighbors >= 1000)
                << "Maximum number of allowed neighbours reached when "
                   "searching for neighbours with radii = "
                << radius << " [ max number of allowed neighbours = 1000 ].\n";

            for (IndexType i_neighbour = 0; i_neighbour < number_of_neighbors; ++i_neighbour) {
                const auto& r_neighbour_point = *rTLS.mNeighbourEntityPoints[i_neighbour];

                double squared_distance;
                auto p_nearest_damped_entity_point = comp_search_tree.SearchNearestPoint(r_neighbour_point, squared_distance);

                const double weight = kernel_function.ComputeWeight(radius, std::sqrt(squared_distance));

                *(rOutput.data().begin() + Index * number_of_entities + r_neighbour_point.Id()) = weight;
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
