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
        "damping_type"               : "nearest_entity",
        "damping_function_type"      : "cosine",
        "damping_max_items_in_bucket": 10,
        "damped_model_part_settings" : {}
    })" );

    Settings.ValidateAndAssignDefaults(defaults);

    mLeafMaxSize = Settings["damping_max_items_in_bucket"].GetInt();

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

    #if defined(_WIN32) || defined(_WIN64)
        // The MSVC compiler does not support std::unique_lock<std::mutex> lock(mutex);
        // in the case the mutex is passed by reference. [See https://stackoverflow.com/questions/78598141/first-stdmutexlock-crashes-in-application-built-with-latest-visual-studio ]
        const unsigned int number_of_threads = 1;
        KRATOS_WARNING("NearestEntityExplicitDamping") << "The nanoflann will construct the KD tree in serial mode in Windows.";
    #else
        const unsigned int number_of_threads = 0;
    #endif

    const auto stride = this->GetStride();
    auto& r_container = *(std::get<typename TContainerType::Pointer>(mpDampingRadius->GetContainer()));
    const auto radius_view = mpDampingRadius->ViewData();

    mDampingCoefficients.resize(radius_view.size(), stride);
    auto& damping_coefficient_data = mDampingCoefficients.data();

    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        auto& r_damped_model_parts = mComponentWiseDampedModelParts[i_comp];

        PositionAdapter adapter(r_damped_model_parts);

        if (adapter.kdtree_get_point_count() == 0) {
            std::fill(damping_coefficient_data.begin(), damping_coefficient_data.end(), 1.0);
        } else {
            // now construct the kd tree
            KDTreeIndexType kd_tree_index(
                3, adapter,
                nanoflann::KDTreeSingleIndexAdaptorParams(
                    mLeafMaxSize, nanoflann::KDTreeSingleIndexAdaptorFlags::None, number_of_threads));

            const auto& kernel_function = *mpKernelFunction;

            // now calculate the damping for each entity
            IndexPartition<IndexType>(r_container.size()).for_each([&damping_coefficient_data, &r_container, &kd_tree_index, &kernel_function, &radius_view, stride, i_comp](const auto Index){
                const auto data_begin_index = Index * stride;
                const auto radius = radius_view[Index];

                double squared_distance;
                unsigned int global_index;
                kd_tree_index.knnSearch(&OptimizationUtils::GetEntityPosition(*(r_container.begin() + Index))[0], 1, &global_index, &squared_distance);

                double& value = *(damping_coefficient_data.begin() + data_begin_index + i_comp);
                value = 1.0 - kernel_function.ComputeWeight(radius, std::sqrt(squared_distance));
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

    const IndexType stride = this->GetStride();

    rDampedWeights.resize(stride);

    for (IndexType i_comp = 0; i_comp < stride; ++i_comp) {
        auto& r_damped_weights = rDampedWeights[i_comp];
        r_damped_weights.resize(rNeighbours.size());
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
