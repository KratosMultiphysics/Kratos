//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/data_containers.h"
#include "custom_utilities/generic_reduction_utilities.h"
#include "custom_utilities/method_utilities.h"
#include "custom_utilities/norms.h"

// Include base h
#include "spatial_methods.h"

namespace Kratos
{

namespace SpatialMethodHelperUtilities {

template <class TDataType>
void ResizeAndInitialize(
    TDataType& rOutput,
    const TDataType& rValueForSize,
    const typename DataTypeTraits<TDataType>::RawDataType& rInitializationValue)
{
    DataTypeTraits<TDataType>::Resize(rOutput, rValueForSize);
    DataTypeTraits<TDataType>::Initialize(rOutput, rInitializationValue);
}

template <class TDataType, class TIteratorType>
void ResizeAndFill(
    TDataType& rOutput,
    const TDataType& rValueForSize,
    TIteratorType Begin)
{
    const auto number_of_components = DataTypeTraits<TDataType>::Size(rValueForSize);
    DataTypeTraits<TDataType>::Resize(rOutput, rValueForSize);
    DataTypeTraits<TDataType>::FillFromVector(rOutput, Begin, Begin + number_of_components);
}

} // namespace SpatialMethodHelperUtilities

template <class TDataType>
SpatialMethods::DistributionInfoType SpatialMethods::Distribution(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const std::string& rNormType,
    const DataLocation& rLocation,
    Parameters Params)
{
    KRATOS_TRY

    const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, rLocation);

    const auto r_norm_type = Norms::GetNorm<TDataType>(rNormType);

    return std::visit(
        [&](auto& rDataContainer, auto& rNorm) -> DistributionInfoType {
            using data_container_type = std::decay_t<decltype(rDataContainer)>;
            using norm_type = std::decay_t<decltype(rNorm)>;
            using norm_return_type = typename norm_type::ResultantValueType<TDataType>;
            using data_type_traits = DataTypeTraits<norm_return_type>;

            Parameters default_parameters = Parameters(R"(
            {
                "number_of_value_groups" : 10,
                "min_value"              : "min",
                "max_value"              : "max"
            })");

            // first fix the min value
            if (Params.Has("min_value") && Params["min_value"].IsDouble()) {
                default_parameters["min_value"].SetDouble(0.0);
            }
            if (Params.Has("max_value") && Params["max_value"].IsDouble()) {
                default_parameters["max_value"].SetDouble(0.0);
            }
            Params.RecursivelyValidateAndAssignDefaults(default_parameters);

            DistributionInfo<norm_return_type> distribution_info;

            auto& min_value = distribution_info.mMin;
            if (Params["min_value"].IsDouble()) {
                data_type_traits::Initialize(min_value, Params["min_value"].GetDouble());
            } else if (Params["min_value"].IsString() && Params["min_value"].GetString() == "min") {
                min_value = std::get<0>(GenericReductionUtilities::GenericReduction<data_container_type, norm_type, MinOperation, true>(rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm).GetValue());
            } else {
                KRATOS_ERROR
                    << "Unknown min_value. Allowed only double or \"min\" "
                       "string as a value. [ min_value = "
                    << Params["min_value"] << " ]\n.";
            }

            // first fix the max value
            auto& max_value = distribution_info.mMax;
            if (Params["max_value"].IsDouble()) {
                data_type_traits::Initialize(max_value, Params["max_value"].GetDouble());
            }
            else if (Params["max_value"].IsString() && Params["max_value"].GetString() == "max") {
                max_value = std::get<0>(GenericReductionUtilities::GenericReduction<data_container_type, norm_type, MaxOperation, true>(rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm).GetValue());
            } else {
                KRATOS_ERROR
                    << "Unknown max_value. Allowed only double or \"max\" "
                       "string as a value. [ max_value = "
                    << Params["max_value"] << " ]\n.";
            }

            auto& group_limits = distribution_info.mGroupUpperValues;
            const IndexType number_of_groups = Params["number_of_value_groups"].GetInt();

            // we need additional two groups to store values below the specified
            // minimum and values above the specified maximum.
            const IndexType number_of_all_groups = number_of_groups + 2;
            group_limits.resize(number_of_all_groups);
            for (IndexType i = 0; i < number_of_groups + 1; ++i) {
                group_limits[i] = min_value + (max_value - min_value) * static_cast<double>(i) / static_cast<double>(number_of_groups);
            }

            const IndexType number_of_components = DataTypeTraits<norm_return_type>::Size(max_value);

            // final group limit is extended by a small amount. epsilon in numeric
            // limits cannot be used since testing also need to have the same
            // extending value in python. Therefore hard coded value is used
            auto& last_group_limit = group_limits[group_limits.size() - 2];
            for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                DataTypeTraits<norm_return_type>::GetComponent(last_group_limit, i_comp) += 1e-16;
            }
            norm_return_type additional_max_value;
            SpatialMethodHelperUtilities::ResizeAndInitialize(additional_max_value, max_value, std::numeric_limits<typename DataTypeTraits<norm_return_type>::RawDataType>::max());
            group_limits.back() = additional_max_value;

            /// reduction class
            class DistributionReduction {
            public:
                using data_type = std::tuple<IndexType, IndicesType, norm_return_type>;

                using return_type = std::tuple<
                                        std::vector<IndicesType>, // group value counts
                                        std::vector<norm_return_type>, // group means
                                        std::vector<norm_return_type>  // group variances
                                    >;

                return_type mValue;

                /// access to reduced value
                return_type GetValue() const
                {
                    return mValue;
                }

                /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
                void LocalReduce(const data_type& rValue)
                {
                    ResizeAndInitialize(rValue);

                    auto& r_current_group_counts = std::get<0>(mValue);
                    auto& r_current_group_means = std::get<1>(mValue);
                    auto& r_current_group_variances = std::get<2>(mValue);

                    const auto& r_group_indices = std::get<1>(rValue);
                    const auto& r_value = std::get<2>(rValue);
                    const IndexType number_of_components = DataTypeTraits<norm_return_type>::Size(r_value);

                    for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                        const auto group_index = r_group_indices[i_comp];
                        const auto raw_value = DataTypeTraits<norm_return_type>::GetComponent(r_value, i_comp);
                        ++r_current_group_counts[group_index][i_comp];
                        DataTypeTraits<norm_return_type>::GetComponent(r_current_group_means[group_index], i_comp) += raw_value;
                        DataTypeTraits<norm_return_type>::GetComponent(r_current_group_variances[group_index], i_comp) += std::pow(raw_value, 2);
                    }
                }

                /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
                void ThreadSafeReduce(const DistributionReduction& rOther)
                {
                    KRATOS_CRITICAL_SECTION

                    auto& r_group_counts = std::get<0>(mValue);
                    auto& r_group_means = std::get<1>(mValue);
                    auto& r_group_variances = std::get<2>(mValue);

                    const auto& r_other_group_counts = std::get<0>(rOther.mValue);
                    const auto& r_other_group_means = std::get<1>(rOther.mValue);
                    const auto& r_other_group_variances = std::get<2>(rOther.mValue);

                    const IndexType number_of_groups = r_other_group_counts.size();

                    // first check and resize everything
                    if (r_group_counts.size() != number_of_groups) {
                        r_group_counts.resize(number_of_groups);
                        r_group_means.resize(number_of_groups);
                        r_group_variances.resize(number_of_groups);
                    }

                    for (IndexType i_group = 0; i_group < number_of_groups; ++i_group) {
                        auto& r_group_count = r_group_counts[i_group];
                        auto& r_group_mean = r_group_means[i_group];
                        auto& r_group_variance = r_group_variances[i_group];

                        const auto& r_other_group_count = r_other_group_counts[i_group];
                        const auto& r_other_group_mean = r_other_group_means[i_group];
                        const auto& r_other_group_variance =
                            r_other_group_variances[i_group];

                        const IndexType number_of_components = r_other_group_count.size();

                        if (r_group_count.size() != number_of_components) {
                            r_group_count.resize(number_of_components);
                            std::fill(r_group_count.begin(), r_group_count.end(), 0);
                            SpatialMethodHelperUtilities::ResizeAndInitialize(r_group_mean, r_other_group_mean, 0.0);
                            SpatialMethodHelperUtilities::ResizeAndInitialize(r_group_variance, r_other_group_variance, 0.0);
                        }

                        for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                            r_group_count[i_comp] += r_other_group_count[i_comp];
                        }
                        r_group_mean += r_other_group_mean;
                        r_group_variance += r_other_group_variance;
                    }
                }

                void ResizeAndInitialize(const data_type& rValue)
                {
                    const auto number_of_groups = std::get<0>(rValue);
                    const auto& r_value = std::get<2>(rValue);
                    const auto number_of_components = DataTypeTraits<norm_return_type>::Size(r_value);

                    auto& r_current_group_counts = std::get<0>(mValue);
                    auto& r_current_group_means = std::get<1>(mValue);
                    auto& r_current_group_variances = std::get<2>(mValue);

                    if (r_current_group_counts.size() != number_of_groups) {
                        // resize
                        r_current_group_counts.resize(number_of_groups);
                        r_current_group_means.resize(number_of_groups);
                        r_current_group_variances.resize(number_of_groups);

                        // initialize
                        for (IndexType i_group = 0; i_group < number_of_groups; ++i_group) {
                            r_current_group_counts[i_group].resize(number_of_components);
                            std::fill(r_current_group_counts[i_group].begin(), r_current_group_counts[i_group].end(), 0);
                            SpatialMethodHelperUtilities::ResizeAndInitialize(r_current_group_means[i_group], r_value, 0.0);
                            SpatialMethodHelperUtilities::ResizeAndInitialize(r_current_group_variances[i_group], r_value, 0.0);
                        }
                    }
                }
            };

            const auto& reuduced_values =
                IndexPartition<IndexType>(rDataContainer.Size()).for_each<DistributionReduction>(TDataType{}, [&rDataContainer, &rNorm, &group_limits, number_of_components, number_of_groups](const IndexType Index, TDataType& rTLS) {
                    rDataContainer.GetValue(rTLS, Index);
                    auto norm_value = rNorm.Evaluate(rTLS);

                    IndicesType indices(number_of_components);
                    for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                        const auto& comp_value = DataTypeTraits<norm_return_type>::GetComponent(norm_value, i_comp);
                        for (IndexType i_group = 0; i_group < group_limits.size(); ++i_group) {
                            if (comp_value < DataTypeTraits<norm_return_type>::GetComponent(group_limits[i_group], i_comp)) {
                                indices[i_comp] = i_group;
                                break;
                            }
                        }
                    }

                    return std::make_tuple(group_limits.size(), indices, norm_value);
                });

            // now prepare data for mpi communication
            std::vector<typename DataTypeTraits<norm_return_type>::RawDataType> local_values;
            local_values.resize(number_of_components * number_of_all_groups * 2);
            auto values_begin = local_values.begin();
            IndicesType local_distribution;
            local_distribution.resize(number_of_components * number_of_all_groups);
            auto indices_begin = local_distribution.begin();
            for (IndexType i_group = 0; i_group < number_of_all_groups; ++i_group) {
                DataTypeTraits<IndicesType>::FillToVector(indices_begin, std::get<0>(reuduced_values)[i_group]); // get the group counts
                DataTypeTraits<norm_return_type>::FillToVector(values_begin, std::get<1>(reuduced_values)[i_group]); // put the means
                DataTypeTraits<norm_return_type>::FillToVector(values_begin + number_of_components, std::get<2>(reuduced_values)[i_group]); // put the variances
                values_begin += 2 * number_of_components;
                indices_begin += number_of_components;
            }

            // now do the mpi communicaton
            const auto& global_distribution = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_distribution);
            const auto& global_values = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_values);
            const double number_of_items = static_cast<double>(std::max(std::accumulate(global_distribution.begin(), global_distribution.end(), 0), 1)) / number_of_components;

            // now revert back to the group values
            distribution_info.mGroupNumberOfValues.resize(number_of_all_groups);
            distribution_info.mGroupValueDistributionPercentage.resize(number_of_all_groups);
            distribution_info.mGroupMean.resize(number_of_all_groups);
            distribution_info.mGroupVariance.resize(number_of_all_groups);
            auto global_indices_begin = global_distribution.begin();
            auto global_values_begin = global_values.begin();
            for (IndexType i_group = 0; i_group < number_of_all_groups; ++i_group) {
                auto& current_number_of_values = distribution_info.mGroupNumberOfValues[i_group];
                current_number_of_values.resize(number_of_components);
                DataTypeTraits<IndicesType>::FillFromVector(current_number_of_values, global_indices_begin, global_indices_begin + number_of_components);

                auto& current_distribution_percentage = distribution_info.mGroupValueDistributionPercentage[i_group];
                SpatialMethodHelperUtilities::ResizeAndFill(current_distribution_percentage, max_value, global_indices_begin);

                auto& current_mean = distribution_info.mGroupMean[i_group];
                SpatialMethodHelperUtilities::ResizeAndFill(current_mean, max_value, global_values_begin);

                auto& current_variance = distribution_info.mGroupVariance[i_group];
                SpatialMethodHelperUtilities::ResizeAndFill(current_variance, max_value, global_values_begin + number_of_components);

                global_indices_begin += number_of_components;
                global_values_begin += number_of_components * 2;

                // post processing of values
                current_distribution_percentage /= number_of_items;
                for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                    const auto n = current_number_of_values[i_comp];
                    if (n > 0) {
                        DataTypeTraits<norm_return_type>::GetComponent(current_mean, i_comp) /= n;
                        DataTypeTraits<norm_return_type>::GetComponent(current_variance, i_comp) /= n;
                    }
                    DataTypeTraits<norm_return_type>::GetComponent(current_variance, i_comp) -= std::pow(DataTypeTraits<norm_return_type>::GetComponent(current_mean, i_comp), 2);
                }
            }

            // reversing group limit extention
            for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                DataTypeTraits<norm_return_type>::GetComponent(group_limits[group_limits.size() - 2], i_comp) += 1e-16;
            }
            group_limits[group_limits.size() - 1] = max_value;

            return distribution_info;
        }, data_container, r_norm_type);

    KRATOS_CATCH("");
}

// template instantiations
template SpatialMethods::DistributionInfoType SpatialMethods::Distribution(const ModelPart&, const Variable<double>&, const std::string&, const DataLocation&, Parameters);
template SpatialMethods::DistributionInfoType SpatialMethods::Distribution(const ModelPart&, const Variable<array_1d<double, 3>>&, const std::string&, const DataLocation&, Parameters);
template SpatialMethods::DistributionInfoType SpatialMethods::Distribution(const ModelPart&, const Variable<array_1d<double, 4>>&, const std::string&, const DataLocation&, Parameters);
template SpatialMethods::DistributionInfoType SpatialMethods::Distribution(const ModelPart&, const Variable<array_1d<double, 6>>&, const std::string&, const DataLocation&, Parameters);
template SpatialMethods::DistributionInfoType SpatialMethods::Distribution(const ModelPart&, const Variable<array_1d<double, 9>>&, const std::string&, const DataLocation&, Parameters);
template SpatialMethods::DistributionInfoType SpatialMethods::Distribution(const ModelPart&, const Variable<Vector>&, const std::string&, const DataLocation&, Parameters);
template SpatialMethods::DistributionInfoType SpatialMethods::Distribution(const ModelPart&, const Variable<Matrix>&, const std::string&, const DataLocation&, Parameters);

} // namespace Kratos