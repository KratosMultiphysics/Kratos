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
void NearestEntityExplicitDamping<TContainerType>::SetRadius(TensorAdaptor<double>::Pointer pDampingRadiusTensorAdaptor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(std::holds_alternative<typename TContainerType::Pointer>(pDampingRadiusTensorAdaptor->GetContainer()))
        << "Radius container type and the explicit damping type container mismatch [ "
        << "tensor adaptor = " << *pDampingRadiusTensorAdaptor << " ].\n";

    mpDampingRadius = pDampingRadiusTensorAdaptor;

    KRATOS_CATCH("");
}

template <class TContainerType>
TensorAdaptor<double>::Pointer NearestEntityExplicitDamping<TContainerType>::GetRadius() const
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
    auto& r_container = *(std::get<typename TContainerType::Pointer>(mpDampingRadius->GetContainer()));
    const auto radius_view = mpDampingRadius->ViewData();

    mDampingCoefficients.resize(radius_view.size(), stride);

    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        auto& r_damped_model_parts = mComponentWiseDampedModelParts[i_comp];

        EntityPointVector damped_entities;
        damped_entities.resize(std::accumulate(
            r_damped_model_parts.begin(), r_damped_model_parts.end(), 0UL,
            [](const IndexType Value, const auto& pModelPart) {
                return Value + ModelPartUtils::GetContainer<TContainerType>(*pModelPart).size();
            }));

        IndexType local_start = 0;
        for (const auto& p_model_part : r_damped_model_parts) {
            const auto& r_damping_container = ModelPartUtils::GetContainer<TContainerType>(*p_model_part);
            IndexPartition<IndexType>(r_damping_container.size()).for_each([&damped_entities, &r_damping_container, local_start](const auto Index) {
                *(damped_entities.begin() + local_start + Index) =  Kratos::make_shared<EntityPointType>(*(r_damping_container.begin() + Index), Index);
            });
            local_start += r_damping_container.size();
        }

        if (damped_entities.empty()) {
            IndexPartition<IndexType>(r_container.size()).for_each([this, i_comp](const auto Index) {
                this->mDampingCoefficients(Index, i_comp) = 1.0;
            });
        } else {
            // now construct the kd tree
            auto p_search_tree =  Kratos::make_shared<KDTree>(damped_entities.begin(), damped_entities.end(), mBucketSize);
            const auto& kernel_function = *mpKernelFunction;

            // now calculate the damping for each entity
            IndexPartition<IndexType>(r_container.size()).for_each([this, &r_container, &p_search_tree, &kernel_function, &radius_view, stride, i_comp](const auto Index){
                EntityPointType entity_point(*(r_container.begin() + Index), Index);
                const auto radius = radius_view[Index];

                double squared_distance;
                auto p_nearest_damped_entity_point = p_search_tree->SearchNearestPoint(entity_point, squared_distance);

                this->mDampingCoefficients(Index, i_comp) = 1.0 - kernel_function.ComputeWeight(radius, std::sqrt(squared_distance));
            });
        }
    }

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

    for (IndexType i_comp = 0; i_comp < this->GetStride(); ++i_comp) {
        auto& r_damped_weights = rDampedWeights[i_comp];
        for (IndexType i_neighbour = 0; i_neighbour < NumberOfNeighbours; ++i_neighbour) {
            r_damped_weights[i_neighbour] = rWeights[i_neighbour] * mDampingCoefficients(Index, i_comp);
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
    const auto number_of_entities = mDampingCoefficients.size1();

    KRATOS_ERROR_IF_NOT(ComponentIndex < stride)
        << "Invalid component index [ component index = " << ComponentIndex
        << ", stride = " << stride << " ].\n";

    if (rOutput.size1() != number_of_entities || rOutput.size2() != number_of_entities) {
        rOutput.resize(number_of_entities, number_of_entities, false);
    }

    rOutput.clear();

    IndexPartition<IndexType>(number_of_entities).for_each([this, &rOutput, ComponentIndex, number_of_entities](const auto Index) {
        *(rOutput.data().begin() + Index * number_of_entities + Index) = this->mDampingCoefficients(Index, ComponentIndex);
    });

    KRATOS_CATCH("");
}

// template instantiations
template class NearestEntityExplicitDamping<ModelPart::NodesContainerType>;
template class NearestEntityExplicitDamping<ModelPart::ConditionsContainerType>;
template class NearestEntityExplicitDamping<ModelPart::ElementsContainerType>;

} // namespace Kratos.
