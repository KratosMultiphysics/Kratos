//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <numeric>
#include <tuple>
#include <algorithm>

// External includes

// Project includes
#include "utilities/model_part_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "expression/literal_expression.h"
#include "expression/literal_flat_expression.h"

// Application includes
#include "custom_utilities/filtering/filter_utils.h"
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/smooth_clamper.h"

// Include base h
#include "mask_utils.h"

namespace Kratos {

template<class TContainerType>
void MaskUtils::CheckCompatibility(
    const ContainerExpression<TContainerType>& rMask1,
    const ContainerExpression<TContainerType>& rMask2)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rMask1.GetItemComponentCount() == 1)
        << "rMask1 should be a scalar expression. [ shape of the given expression = "
        << rMask1.GetItemShape() << " ].\n";

    KRATOS_ERROR_IF_NOT(rMask2.GetItemComponentCount() == 1)
        << "rMask2 should be a scalar expression. [ shape of the given expression = "
        << rMask2.GetItemShape() << " ].\n";

    KRATOS_ERROR_IF_NOT(rMask1.GetContainer().size() == rMask2.GetContainer().size())
        << "rMask1 and rMask2 entities size mismatch [ rMask1.size() = " << rMask1.GetContainer().size()
        << ", rMask2.size() = " << rMask2.GetContainer().size() << " ].\n";

    KRATOS_CATCH("");
}

template<class TContainerType>
std::size_t MaskUtils::GetMaskSize(
    const ContainerExpression<TContainerType>& rMask,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rMask.GetItemComponentCount() == 1)
        << "Mask should be a scalar expression. [ shape of the given expression = "
        << rMask.GetItemShape() << " ].\n";

    const auto& r_expression = rMask.GetExpression();

    return rMask.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(IndexPartition<IndexType>(r_expression.NumberOfEntities()).for_each<SumReduction<int>>([&r_expression, RequiredMinimumRedundancy](const auto Index) {
        return r_expression.Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy;
    }));

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> MaskUtils::GetMask(
    const ContainerExpression<TContainerType>& rScalarExpression)
{
    KRATOS_TRY

    const auto& r_input_expression = rScalarExpression.GetExpression();
    const auto number_of_entities = r_input_expression.NumberOfEntities();

    KRATOS_ERROR_IF_NOT(r_input_expression.GetItemComponentCount() == 1)
        << "rScalarExpression should be a scalar expression. [ shape of the given expression = "
        << r_input_expression.GetItemShape() << " ].\n";

    struct Data
    {
        IndexType mIndex;
        double mValue;
        bool operator<(const Data& rRight) const { return mValue < rRight.mValue; }
    };

    std::vector<Data> index_value_pairs_vector;
    index_value_pairs_vector.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&index_value_pairs_vector, &r_input_expression](const auto Index) {
        index_value_pairs_vector[Index].mIndex = Index;
        index_value_pairs_vector[Index].mValue = r_input_expression.Evaluate(Index, Index, 0);
    });

    // now sort expression values
    std::sort(index_value_pairs_vector.begin(), index_value_pairs_vector.end(), [](const auto& rV1, const auto& rV2){
        return rV1.mValue > rV2.mValue;
    });

    // now find the coverage
    const auto& r_data = IndexPartition<IndexType>(number_of_entities).for_each<MaxReduction<Data>>([&index_value_pairs_vector](const auto Index){
        const double current_sum = std::accumulate(index_value_pairs_vector.begin(), index_value_pairs_vector.begin() + Index + 1, 0.0, [](const auto& rValue, const auto& rIndexValuePair) {
                                                return rValue + rIndexValuePair.mValue;
                                            });
        return Data{Index, current_sum / std::sqrt(Index + 1)};
    });

    auto p_expression = LiteralFlatExpression<int>::Create(number_of_entities, {});
    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression](const auto Index) {
        *(p_expression->begin() + Index) = 0;
    });
    IndexPartition<IndexType>(r_data.mIndex + 1).for_each([&p_expression, &index_value_pairs_vector](const auto Index){
        *(p_expression->begin() + index_value_pairs_vector[Index].mIndex) = 1;
    });

    auto result = rScalarExpression;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
double MaskUtils::GetMaskThreshold(const ContainerExpression<TContainerType>& rScalarExpression)
{
    KRATOS_TRY

    const auto& r_input_expression = rScalarExpression.GetExpression();
    const auto number_of_entities = r_input_expression.NumberOfEntities();

    KRATOS_ERROR_IF_NOT(r_input_expression.GetItemComponentCount() == 1)
        << "rScalarExpression should be a scalar expression. [ shape of the given expression = "
        << r_input_expression.GetItemShape() << " ].\n";

    struct Data
    {
        IndexType mIndex;
        double mValue;
        bool operator<(const Data& rRight) const { return mValue < rRight.mValue; }
    };

    std::vector<Data> index_value_pairs_vector;
    index_value_pairs_vector.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&index_value_pairs_vector, &r_input_expression](const auto Index) {
        index_value_pairs_vector[Index].mIndex = Index;
        index_value_pairs_vector[Index].mValue = r_input_expression.Evaluate(Index, Index, 0);
    });

    // now sort expression values
    std::sort(index_value_pairs_vector.begin(), index_value_pairs_vector.end(), [](const auto& rV1, const auto& rV2){
        return rV1.mValue > rV2.mValue;
    });

    // now find the coverage
    const auto& r_data = IndexPartition<IndexType>(number_of_entities).for_each<MaxReduction<Data>>([&index_value_pairs_vector](const auto Index){
        const double current_sum = std::accumulate(index_value_pairs_vector.begin(), index_value_pairs_vector.begin() + Index + 1, 0.0, [](const auto& rValue, const auto& rIndexValuePair) {
                                                return rValue + rIndexValuePair.mValue;
                                            });
        return Data{Index, current_sum / std::sqrt(Index + 1)};
    });

    if (r_data.mIndex < number_of_entities - 1) {
        const auto threshold_value_index_1 = index_value_pairs_vector[r_data.mIndex].mIndex;
        const auto threshold_value_index_2 = index_value_pairs_vector[r_data.mIndex + 1].mIndex;
        return 0.5 * (r_input_expression.Evaluate(threshold_value_index_1, threshold_value_index_1, 0) + r_input_expression.Evaluate(threshold_value_index_2, threshold_value_index_2, 0));
    } else {
        const auto threshold_value_index = index_value_pairs_vector[r_data.mIndex].mIndex;
        return r_input_expression.Evaluate(threshold_value_index, threshold_value_index, 0);
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> MaskUtils::GetMask(
    const ContainerExpression<TContainerType>& rScalarExpression,
    const double Threshold)
{
    KRATOS_TRY

    const auto& r_input_expression = rScalarExpression.GetExpression();
    const auto number_of_entities = r_input_expression.NumberOfEntities();

    KRATOS_ERROR_IF_NOT(r_input_expression.GetItemComponentCount() == 1)
        << "rScalarExpression should be a scalar expression. [ shape of the given expression = "
        << r_input_expression.GetItemShape() << " ].\n";

    auto p_expression = LiteralFlatExpression<int>::Create(number_of_entities, {});
    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression, &r_input_expression, Threshold](const auto Index) {
        *(p_expression->begin() + Index) = r_input_expression.Evaluate(Index, Index, 0) > Threshold;
    });

    auto result = rScalarExpression;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> MaskUtils::Union(
    const ContainerExpression<TContainerType>& rMask1,
    const ContainerExpression<TContainerType>& rMask2,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rMask1, rMask2);

    const auto& r_mask_1_exp = rMask1.GetExpression();
    const auto& r_mask_2_exp = rMask2.GetExpression();
    const auto number_of_entities = r_mask_1_exp.NumberOfEntities();

    auto p_expression = LiteralFlatExpression<int>::Create(number_of_entities, {});
    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression, &r_mask_1_exp, &r_mask_2_exp, RequiredMinimumRedundancy](const auto Index) {
        *(p_expression->begin() + Index) = (r_mask_1_exp.Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy ||
                                            r_mask_2_exp.Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy)
                                                ? RequiredMinimumRedundancy
                                                : 0;
    });

    auto result = rMask1;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> MaskUtils::Intersect(
    const ContainerExpression<TContainerType>& rMask1,
    const ContainerExpression<TContainerType>& rMask2,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rMask1, rMask2);

    const auto& r_mask_1_exp = rMask1.GetExpression();
    const auto& r_mask_2_exp = rMask2.GetExpression();
    const auto number_of_entities = r_mask_1_exp.NumberOfEntities();

    auto p_expression = LiteralFlatExpression<int>::Create(number_of_entities, {});
    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression, &r_mask_1_exp, &r_mask_2_exp, RequiredMinimumRedundancy](const auto Index) {
        *(p_expression->begin() + Index) = (r_mask_1_exp.Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy &&
                                            r_mask_2_exp.Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy)
                                                ? RequiredMinimumRedundancy
                                                : 0;
    });

    auto result = rMask1;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> MaskUtils::Subtract(
    const ContainerExpression<TContainerType>& rMask1,
    const ContainerExpression<TContainerType>& rMask2,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rMask1, rMask2);

    const auto& r_mask_1_exp = rMask1.GetExpression();
    const auto& r_mask_2_exp = rMask2.GetExpression();
    const auto number_of_entities = r_mask_1_exp.NumberOfEntities();

    auto p_expression = LiteralFlatExpression<int>::Create(number_of_entities, {});
    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression, &r_mask_1_exp, &r_mask_2_exp, RequiredMinimumRedundancy](const auto Index) {
        *(p_expression->begin() + Index) = (r_mask_1_exp.Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy &&
                                            r_mask_2_exp.Evaluate(Index, Index, 0) < RequiredMinimumRedundancy)
                                                ? RequiredMinimumRedundancy
                                                : 0;
    });

    auto result = rMask1;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> MaskUtils::Scale(
    const ContainerExpression<TContainerType>& rScalarExpression,
    const ContainerExpression<TContainerType>& rMask,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rScalarExpression, rMask);

    const auto& r_scalar_expression = rScalarExpression.GetExpression();
    const auto& r_mask = rMask.GetExpression();
    const auto number_of_entities = r_scalar_expression.NumberOfEntities();

    auto p_expression = LiteralFlatExpression<double>::Create(number_of_entities, {});
    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression, &r_scalar_expression, &r_mask, RequiredMinimumRedundancy](const auto Index) {
        *(p_expression->begin() + Index) = (r_mask.Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy)
                                                ? r_scalar_expression.Evaluate(Index, Index, 0)
                                                : 0;
    });

    auto result = rScalarExpression;
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
std::vector<std::tuple<std::vector<IndexType>, typename ContainerExpression<TContainerType>::Pointer>> MaskUtils::ClusterMasks(
    const std::vector<ContainerExpression<TContainerType>>& rMasksList,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    std::vector<std::tuple<std::vector<IndexType>, typename ContainerExpression<TContainerType>::Pointer>> cluster_data;
    if (rMasksList.size() == 0) {
        return cluster_data;
    }

    const auto& r_front_cexp = rMasksList.front();
    const IndexType number_of_entities = r_front_cexp.GetExpression().NumberOfEntities();

    // Check if all the masks are compatible
    for (const auto& r_mask : rMasksList) {
        KRATOS_ERROR_IF_NOT(number_of_entities == r_mask.GetExpression().NumberOfEntities())
            << "Mismatch in mask size [required mask size = " << number_of_entities << ", "
            << " found one mask with size = " << r_mask.GetExpression().NumberOfEntities() << " ].\n";

        KRATOS_ERROR_IF_NOT(r_mask.GetExpression().GetItemComponentCount() == 1)
            << "Found a mask with a non scalar dimensionality which is not allowed [ requried mask shape = {}, found shape = "
            << r_mask.GetExpression().GetItemShape() << " ].\n";
    }

    std::vector<std::vector<IndexType>> domain_mask_indices(number_of_entities);
    IndexPartition<IndexType>(number_of_entities).for_each([&domain_mask_indices, &rMasksList, RequiredMinimumRedundancy](const auto Index) {
        auto& r_indices_list = domain_mask_indices[Index];
        for (IndexType i = 0; i < rMasksList.size(); ++i) {
            if (rMasksList[i].GetExpression().Evaluate(Index, Index, 0) >= RequiredMinimumRedundancy) {
                r_indices_list.push_back(i);
            }
        }
    });

    // now find unique list of cluster indices
    std::vector<LiteralFlatExpression<int>::Pointer> cluster_mask_exps;
    for (IndexType i = 0; i < number_of_entities; ++i) {
        auto& mask_indices = domain_mask_indices[i];
        std::sort(mask_indices.begin(), mask_indices.end());

        /// check against existing mask indices
        auto p_itr = std::find_if(cluster_data.begin(), cluster_data.end(), [&mask_indices](const auto& rData){ return std::get<0>(rData) == mask_indices; });
        if (p_itr == cluster_data.end()) {
            auto p_expression = LiteralFlatExpression<int>::Create(number_of_entities, {});
            IndexPartition<IndexType>(number_of_entities).for_each([&p_expression](const auto Index){
                *(p_expression->begin() + Index) = 0;
            });

            auto p_container_exp = Kratos::make_shared<ContainerExpression<TContainerType>>(*r_front_cexp.pGetModelPart());
            p_container_exp->SetExpression(p_expression);

            cluster_data.push_back(std::make_tuple(mask_indices, p_container_exp));
            cluster_mask_exps.push_back(p_expression);
        }
    }

    // now fill in the cluster masks
    IndexPartition<IndexType>(number_of_entities).for_each([&cluster_data, &domain_mask_indices, &cluster_mask_exps](const auto Index) {
        const auto& mask_indices = domain_mask_indices[Index];
        auto p_itr = std::find_if(cluster_data.begin(), cluster_data.end(), [&mask_indices](const auto& rData){ return std::get<0>(rData) == mask_indices; });
        const auto cluster_index = std::distance(cluster_data.begin(), p_itr);
        *(cluster_mask_exps[cluster_index]->begin() + Index) = 1;
    });

    return cluster_data;

    KRATOS_CATCH("");
}

template<class TContainerType>
std::vector<IndexType> MaskUtils::GetMasksDividingReferenceMask(
    const ContainerExpression<TContainerType>& rReferenceMask,
    const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rMasksList,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    const auto reference_mask_coverage = GetMaskSize(rReferenceMask, RequiredMinimumRedundancy);

    std::vector<IndexType> indices;
    for (IndexType i = 0; i < rMasksList.size(); ++i) {
        const auto& r_intersected_exp = Intersect(rReferenceMask, *rMasksList[i], RequiredMinimumRedundancy);
        const auto intesection_coverage = GetMaskSize(r_intersected_exp, RequiredMinimumRedundancy);

        if (intesection_coverage > 0 && intesection_coverage < reference_mask_coverage) {
            indices.push_back(i);
        }
    }

    return indices;

    KRATOS_CATCH("");
}

template<class TContainerType>
void MaskUtils::FillModelPartUsingClusterMask(
    ModelPart& rModelPart,
    const ContainerExpression<TContainerType>& rClusterMask,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rModelPart.NumberOfNodes() == 0 || rModelPart.NumberOfElements() == 0 || rModelPart.NumberOfConditions())
        << "The " << rModelPart.FullName() << " is not empty.\n";

    KRATOS_ERROR_IF_NOT(rClusterMask.GetItemComponentCount() == 1)
        << "Cluster mask can only be with a scalar value for every entity. [ rClusterMask = " << rClusterMask << " ].\n";

    const auto number_of_entities = rClusterMask.GetContainer().size();

    TContainerType temp;
    temp.reserve(number_of_entities);
    for (IndexType i_entity = 0; i_entity < number_of_entities; ++i_entity) {
        const auto value = rClusterMask.GetExpression().Evaluate(i_entity, i_entity, 0);
        if (value >= RequiredMinimumRedundancy) {
            temp.insert(temp.end(), *(rClusterMask.GetContainer().ptr_begin() + i_entity));
        }
    }

    if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        rModelPart.AddElements(temp.begin(), temp.end());
        std::vector<ModelPart::NodeType::Pointer> nodes;
        for (auto& r_element : temp) {
            for (auto itr_p_node_begin = r_element.GetGeometry().ptr_begin(); itr_p_node_begin != r_element.GetGeometry().ptr_end(); ++itr_p_node_begin) {
                nodes.push_back(*itr_p_node_begin);
            }
        }
        rModelPart.AddNodes(nodes.begin(), nodes.end());
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        rModelPart.AddConditions(temp.begin(), temp.end());
        std::vector<ModelPart::NodeType::Pointer> nodes;
        for (auto& r_condition : temp) {
            for (auto itr_p_node_begin = r_condition.GetGeometry().ptr_begin(); itr_p_node_begin != r_condition.GetGeometry().ptr_end(); ++itr_p_node_begin) {
                nodes.push_back(*itr_p_node_begin);
            }
        }
        rModelPart.AddNodes(nodes.begin(), nodes.end());
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        rModelPart.AddNodes(temp.begin(), temp.end());
    } else {
        static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported TContainerType.");
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
std::vector<typename ContainerExpression<TContainerType>::Pointer> MaskUtils::SmoothenMasks(
    const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rMaskList,
    const double Radius,
    const IndexType MaxLeafSize)
{
    KRATOS_TRY

    using position_adaptor_type = NanoFlannSingleContainerPositionAdapter<TContainerType>;

    using pointer_vector_type = typename position_adaptor_type::PointerVectorType;

    using distance_metric_type = typename nanoflann::metric_L2_Simple::traits<double, position_adaptor_type>::distance_t;

    using k_d_tree_index_type = nanoflann::KDTreeSingleIndexAdaptor<distance_metric_type, position_adaptor_type, 3>;

    using k_d_tree_thread_local_storage = NanoFlannKDTreeThreadLocalStorage<pointer_vector_type>;

    std::vector<typename ContainerExpression<TContainerType>::Pointer> smooth_mask_list;

    // do nothing if the input array is empty.
    if (rMaskList.empty()) {
        return smooth_mask_list;
    }

    const auto& r_container = rMaskList.front()->GetContainer();

    // first construct the Kd tree for all masks
    position_adaptor_type adaptor(&r_container);
    k_d_tree_index_type k_d_tree_index(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(MaxLeafSize, nanoflann::KDTreeSingleIndexAdaptorFlags::None, 0));
    k_d_tree_index.buildIndex();

    const double radius_square = Radius * Radius;

    SmoothClamper<TContainerType> clamper(0, 1);

    for (const auto& p_mask : rMaskList) {
        const auto& r_container = p_mask->GetContainer();
        const auto& r_mask_exp = p_mask->GetExpression();

        auto p_result = p_mask->Clone();
        auto p_result_exp = LiteralFlatExpression<double>::Create(r_container.size(), p_mask->GetItemShape());
        p_result->SetExpression(p_result_exp);

        IndexPartition<IndexType>(r_container.size()).for_each(k_d_tree_thread_local_storage(), [&k_d_tree_index, &r_container, &r_mask_exp, &p_result_exp, &clamper, radius_square](const IndexType Index, auto& rTLS){
            // search for entities within radius
            k_d_tree_index.radiusSearch(OptimizationUtils::GetEntityPosition(*(r_container.begin() + Index)).data().begin(), radius_square + std::numeric_limits<double>::epsilon(), rTLS.mNeighbourIndicesAndSquaredDistances, nanoflann::SearchParameters());

            auto& r_value = *(p_result_exp->begin() + Index);
            r_value = 1.0;
            for (const auto& r_result : rTLS.mNeighbourIndicesAndSquaredDistances) {
                if (r_mask_exp.Evaluate(r_result.first, r_result.first, 1) == 0) {
                    const auto distance_ratio = std::sqrt(r_result.second / radius_square);
                    r_value = clamper.ProjectForward(distance_ratio);
                    break;
                }
            }
        });

        smooth_mask_list.push_back(p_result);
    }

    return smooth_mask_list;

    KRATOS_CATCH("");
}

// template instantiations
#ifndef KRATOS_SI_APP_MASK_UTILS_INSTANTIATION
#define KRATOS_SI_APP_MASK_UTILS_INSTANTIATION(CONTAINER_TYPE)                                                                                                                             \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void MaskUtils::CheckCompatibility(                                                                                             \
        const ContainerExpression<CONTAINER_TYPE> &,                                                                                                                                       \
        const ContainerExpression<CONTAINER_TYPE> &);                                                                                                                                      \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) std::size_t MaskUtils::GetMaskSize(                                                                                             \
        const ContainerExpression<CONTAINER_TYPE> &, const IndexType);                                                                                                                     \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> MaskUtils::GetMask(                                                                         \
        const ContainerExpression<CONTAINER_TYPE> &);                                                                                                                                      \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) double MaskUtils::GetMaskThreshold(                                                                                             \
        const ContainerExpression<CONTAINER_TYPE> &);                                                                                                                                      \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> MaskUtils::GetMask(                                                                         \
        const ContainerExpression<CONTAINER_TYPE> &, const double);                                                                                                                        \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> MaskUtils::Union(                                                                           \
        const ContainerExpression<CONTAINER_TYPE> &,                                                                                                                                       \
        const ContainerExpression<CONTAINER_TYPE> &, const IndexType);                                                                                                                     \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> MaskUtils::Intersect(                                                                       \
        const ContainerExpression<CONTAINER_TYPE> &,                                                                                                                                       \
        const ContainerExpression<CONTAINER_TYPE> &, const IndexType);                                                                                                                     \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> MaskUtils::Subtract(                                                                        \
        const ContainerExpression<CONTAINER_TYPE> &,                                                                                                                                       \
        const ContainerExpression<CONTAINER_TYPE> &, const IndexType);                                                                                                                     \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> MaskUtils::Scale(                                                                           \
        const ContainerExpression<CONTAINER_TYPE> &,                                                                                                                                       \
        const ContainerExpression<CONTAINER_TYPE> &, const IndexType);                                                                                                                     \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) std::vector<std::tuple<std::vector<IndexType>, typename ContainerExpression<CONTAINER_TYPE>::Pointer>> MaskUtils::ClusterMasks( \
        const std::vector<ContainerExpression<CONTAINER_TYPE>> &, const IndexType);                                                                                                        \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) std::vector<IndexType> MaskUtils::GetMasksDividingReferenceMask(                                                                \
        const ContainerExpression<CONTAINER_TYPE> &,                                                                                                                                       \
        const std::vector<typename ContainerExpression<CONTAINER_TYPE>::Pointer> &,                                                                                                        \
        const IndexType);                                                                                                                                                                  \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) void MaskUtils::FillModelPartUsingClusterMask(                                                                                  \
        ModelPart&,                                                                                                                                                                        \
        const ContainerExpression<CONTAINER_TYPE>&,                                                                                                                                        \
        const IndexType);                                                                                                                                                                  \
    template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) std::vector<ContainerExpression<CONTAINER_TYPE>::Pointer> MaskUtils::SmoothenMasks<CONTAINER_TYPE>(                                             \
    const std::vector<ContainerExpression<CONTAINER_TYPE>::Pointer>&,                                                                                                                      \
    const double,                                                                                                                                                                          \
    const IndexType);                                                                                                                                                                      \

#endif

KRATOS_SI_APP_MASK_UTILS_INSTANTIATION(ModelPart::NodesContainerType)
KRATOS_SI_APP_MASK_UTILS_INSTANTIATION(ModelPart::ConditionsContainerType)
KRATOS_SI_APP_MASK_UTILS_INSTANTIATION(ModelPart::ElementsContainerType)

#undef KRATOS_SI_APP_MASK_UTILS_INSTANTIATION

} // namespace Kratos