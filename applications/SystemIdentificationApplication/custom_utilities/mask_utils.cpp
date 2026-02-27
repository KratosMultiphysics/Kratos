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
#include "utilities/data_type_traits.h"
#include "utilities/model_part_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/filtering/filter_utils.h"
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/smooth_clamper.h"

// Include base h
#include "mask_utils.h"

namespace Kratos {

void MaskUtils::CheckCompatibility(
    const TensorAdaptor<double>& rMask1,
    const TensorAdaptor<double>& rMask2)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rMask1.Shape().size() == 1)
        << "rMask1 should be a scalar tensor adaptor. [ shape of the given tensor adaptor = "
        << rMask1.Shape() << " ].\n";

    KRATOS_ERROR_IF_NOT(rMask2.Shape().size() == 1)
        << "rMask2 should be a scalar tensor adaptor. [ shape of the given tensor adaptor = "
        << rMask2.Shape() << " ].\n";

    KRATOS_ERROR_IF_NOT(rMask1.Size() == rMask2.Size())
        << "rMask1 and rMask2 entities size mismatch [ rMask1.size() = " << rMask1.Size()
        << ", rMask2.size() = " << rMask2.Size() << " ].\n";

    KRATOS_CATCH("");
}

std::size_t MaskUtils::GetMaskSize(
    const TensorAdaptor<double>& rMask,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rMask.Shape().size() == 1)
        << "Mask should be a scalar tensor adaptor. [ shape of the given tensor adaptor = "
        << rMask.Shape() << " ].\n";

    const auto input_data_view = rMask.ViewData();

    return IndexPartition<IndexType>(rMask.Size()).for_each<SumReduction<int>>([&input_data_view, RequiredMinimumRedundancy](const auto Index) {
        return input_data_view[Index] >= RequiredMinimumRedundancy;
    });

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer MaskUtils::GetMask(
    const TensorAdaptor<double>& rScalarTensorAdaptor)
{
    KRATOS_TRY

    const auto input_data_view = rScalarTensorAdaptor.ViewData();
    const auto number_of_entities = rScalarTensorAdaptor.Size();

    KRATOS_ERROR_IF_NOT(rScalarTensorAdaptor.Shape().size() == 1)
        << "rScalarTensorAdaptor should be a scalar tensor adaptor. [ shape of the given tensor adaptor = "
        << rScalarTensorAdaptor.Shape() << " ].\n";

    struct Data
    {
        IndexType mIndex;
        double mValue;
        bool operator<(const Data& rRight) const { return mValue < rRight.mValue; }
    };

    std::vector<Data> index_value_pairs_vector;
    index_value_pairs_vector.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&index_value_pairs_vector, &input_data_view](const auto Index) {
        index_value_pairs_vector[Index].mIndex = Index;
        index_value_pairs_vector[Index].mValue = input_data_view[Index];
    });

    // now sort tensor adaptor values
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

    auto p_result = rScalarTensorAdaptor.Clone();
    auto p_result_data_view = p_result->ViewData();
    std::ranges::fill(p_result_data_view, 0.0);

    IndexPartition<IndexType>(r_data.mIndex + 1).for_each([&p_result_data_view, &index_value_pairs_vector](const auto Index){
        p_result_data_view[index_value_pairs_vector[Index].mIndex] = 1;
    });

    return p_result;

    KRATOS_CATCH("");
}

double MaskUtils::GetMaskThreshold(const TensorAdaptor<double>& rScalarTensorAdaptor)
{
    KRATOS_TRY

    const auto input_data_view = rScalarTensorAdaptor.ViewData();
    const auto number_of_entities = rScalarTensorAdaptor.Size();

    KRATOS_ERROR_IF_NOT(rScalarTensorAdaptor.Shape().size() == 1)
        << "rScalarTensorAdaptor should be a scalar tensor adaptor. [ shape of the given tensor adaptor = "
        << rScalarTensorAdaptor.Shape() << " ].\n";

    struct Data
    {
        IndexType mIndex;
        double mValue;
        bool operator<(const Data& rRight) const { return mValue < rRight.mValue; }
    };

    std::vector<Data> index_value_pairs_vector;
    index_value_pairs_vector.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&index_value_pairs_vector, &input_data_view](const auto Index) {
        index_value_pairs_vector[Index].mIndex = Index;
        index_value_pairs_vector[Index].mValue = input_data_view[Index];
    });

    // now sort tensor adaptor values
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
        return 0.5 * (input_data_view[threshold_value_index_1] + input_data_view[threshold_value_index_2]);
    } else {
        const auto threshold_value_index = index_value_pairs_vector[r_data.mIndex].mIndex;
        return input_data_view[threshold_value_index];
    }

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer MaskUtils::GetMask(
    const TensorAdaptor<double>& rScalarTensorAdaptor,
    const double Threshold)
{
    KRATOS_TRY

    const auto input_data_view = rScalarTensorAdaptor.ViewData();
    const auto number_of_entities = input_data_view.size();

    KRATOS_ERROR_IF_NOT(rScalarTensorAdaptor.Shape().size() == 1)
        << "rScalarTensorAdaptor should be a scalar tensor adaptor. [ shape of the given tensor adaptor = "
        << rScalarTensorAdaptor.Shape() << " ].\n";

    auto p_result = rScalarTensorAdaptor.Clone();
    auto output_data_view = p_result->ViewData();

    IndexPartition<IndexType>(number_of_entities).for_each([&output_data_view, &input_data_view, Threshold](const auto Index) {
        output_data_view[Index] = input_data_view[Index] > Threshold;
    });

    return p_result;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer MaskUtils::Union(
    const TensorAdaptor<double>& rMask1,
    const TensorAdaptor<double>& rMask2,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rMask1, rMask2);

    const auto& r_mask_1_data_view = rMask1.ViewData();
    const auto& r_mask_2_data_view = rMask2.ViewData();
    const auto number_of_entities = r_mask_1_data_view.size();

    auto p_result = rMask1.Clone();
    auto output_data_view = p_result->ViewData();
    IndexPartition<IndexType>(number_of_entities).for_each([&output_data_view, &r_mask_1_data_view, &r_mask_2_data_view, RequiredMinimumRedundancy](const auto Index) {
        output_data_view[Index] = (r_mask_1_data_view[Index] >= RequiredMinimumRedundancy ||
                                   r_mask_2_data_view[Index] >= RequiredMinimumRedundancy)
                                                ? RequiredMinimumRedundancy
                                                : 0;
    });

    return p_result;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer  MaskUtils::Intersect(
    const TensorAdaptor<double>& rMask1,
    const TensorAdaptor<double>& rMask2,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rMask1, rMask2);

    const auto& r_mask_1_data_view = rMask1.ViewData();
    const auto& r_mask_2_data_view = rMask2.ViewData();
    const auto number_of_entities = r_mask_1_data_view.size();

    auto p_result = rMask1.Clone();
    auto output_data_view = p_result->ViewData();
    IndexPartition<IndexType>(number_of_entities).for_each([&output_data_view, &r_mask_1_data_view, &r_mask_2_data_view, RequiredMinimumRedundancy](const auto Index) {
            output_data_view[Index] = (r_mask_1_data_view[Index] >= RequiredMinimumRedundancy &&
                                       r_mask_2_data_view[Index] >= RequiredMinimumRedundancy)
                                                ? RequiredMinimumRedundancy
                                                : 0;
        });

    return p_result;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer MaskUtils::Subtract(
    const TensorAdaptor<double>& rMask1,
    const TensorAdaptor<double>& rMask2,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rMask1, rMask2);

    const auto& r_mask_1_data_view = rMask1.ViewData();
    const auto& r_mask_2_data_view = rMask2.ViewData();
    const auto number_of_entities = r_mask_1_data_view.size();

    auto p_result = rMask1.Clone();
    auto output_data_view = p_result->ViewData();
    IndexPartition<IndexType>(number_of_entities).for_each([&output_data_view, &r_mask_1_data_view, &r_mask_2_data_view, RequiredMinimumRedundancy](const auto Index) {
        output_data_view[Index] = (r_mask_1_data_view[Index] >= RequiredMinimumRedundancy &&
                                   r_mask_2_data_view[Index] < RequiredMinimumRedundancy)
                                                ? RequiredMinimumRedundancy
                                                : 0;
    });

    return p_result;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer MaskUtils::Scale(
    const TensorAdaptor<double>& rScalarTensorAdaptor,
    const TensorAdaptor<double>& rMask,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    CheckCompatibility(rScalarTensorAdaptor, rMask);

    const auto& input_data_view = rScalarTensorAdaptor.ViewData();
    const auto& mask_data_view = rMask.ViewData();
    const auto number_of_entities = input_data_view.size();

    auto p_result = rScalarTensorAdaptor.Clone();
    auto output_data_view = p_result->ViewData();
    IndexPartition<IndexType>(number_of_entities).for_each([&output_data_view, &input_data_view, &mask_data_view, RequiredMinimumRedundancy](const auto Index) {
        output_data_view[Index] = (mask_data_view[Index] >= RequiredMinimumRedundancy)
                                                ? input_data_view[Index]
                                                : 0;
    });

    return p_result;

    KRATOS_CATCH("");
}

std::vector<std::tuple<std::vector<IndexType>, TensorAdaptor<double>::Pointer>> MaskUtils::ClusterMasks(
    const std::vector<TensorAdaptor<double>::Pointer>& rMasksList,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    std::vector<std::tuple<std::vector<IndexType>, TensorAdaptor<double>::Pointer>> cluster_data;
    if (rMasksList.size() == 0) {
        return cluster_data;
    }

    const auto& r_front_ta = rMasksList.front();
    const IndexType number_of_entities = r_front_ta->Size();

    // Check if all the masks are compatible
    for (const auto& r_mask : rMasksList) {
        KRATOS_ERROR_IF_NOT(number_of_entities == r_mask->Size())
            << "Mismatch in mask size [required mask size = " << number_of_entities << ", "
            << " found one mask with size = " << r_mask->Size() << " ].\n";

        KRATOS_ERROR_IF_NOT(r_mask->Shape().size() == 1)
            << "Found a mask with a non scalar dimensionality which is not allowed [ required mask shape = {}, found shape = "
            << r_mask->Shape() << " ].\n";
    }

    std::vector<std::vector<IndexType>> domain_mask_indices(number_of_entities);
    IndexPartition<IndexType>(number_of_entities).for_each([&domain_mask_indices, &rMasksList, RequiredMinimumRedundancy](const auto Index) {
        auto& r_indices_list = domain_mask_indices[Index];
        for (IndexType i = 0; i < rMasksList.size(); ++i) {
            if (rMasksList[i]->ViewData()[Index] >= RequiredMinimumRedundancy) {
                r_indices_list.push_back(i);
            }
        }
    });

    // now find unique list of cluster indices
    std::vector<NDData<double>::Pointer> cluster_mask_nd_data;
    for (IndexType i = 0; i < number_of_entities; ++i) {
        auto& mask_indices = domain_mask_indices[i];
        std::sort(mask_indices.begin(), mask_indices.end());

        /// check against existing mask indices
        auto p_itr = std::find_if(cluster_data.begin(), cluster_data.end(), [&mask_indices](const auto& rData){ return std::get<0>(rData) == mask_indices; });
        if (p_itr == cluster_data.end()) {
            auto p_nd_data = Kratos::make_shared<NDData<double>>(Kratos::DenseVector<unsigned int>(1, number_of_entities));
            auto p_nd_data_view = p_nd_data->ViewData();
            std::fill(p_nd_data_view.begin(), p_nd_data_view.end(), 0.0);
            cluster_mask_nd_data.push_back(p_nd_data);

            auto p_cluster_ta = Kratos::make_shared<TensorAdaptor<double>>(r_front_ta->GetContainer(), p_nd_data, false);
            cluster_data.push_back(std::make_tuple(mask_indices, p_cluster_ta));
        }
    }

    // now fill in the cluster masks
    IndexPartition<IndexType>(number_of_entities).for_each([&cluster_data, &domain_mask_indices, &cluster_mask_nd_data](const auto Index) {
        const auto& mask_indices = domain_mask_indices[Index];
        auto p_itr = std::find_if(cluster_data.begin(), cluster_data.end(), [&mask_indices](const auto& rData){ return std::get<0>(rData) == mask_indices; });
        const auto cluster_index = std::distance(cluster_data.begin(), p_itr);
        cluster_mask_nd_data[cluster_index]->ViewData()[Index] = 1;
    });

    return cluster_data;

    KRATOS_CATCH("");
}

void MaskUtils::FillModelPartUsingClusterMask(
    ModelPart& rModelPart,
    const TensorAdaptor<double>& rClusterMask,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rModelPart.NumberOfNodes() == 0 || rModelPart.NumberOfElements() == 0 || rModelPart.NumberOfConditions())
        << "The " << rModelPart.FullName() << " is not empty.\n";

    KRATOS_ERROR_IF_NOT(rClusterMask.Shape().size() == 1)
        << "Cluster mask can only be with a scalar value for every entity. [ rClusterMask = " << rClusterMask << " ].\n";

    const auto number_of_entities = rClusterMask.Size();

    const auto cluster_data = rClusterMask.ViewData();
    std::visit([&cluster_data, &rModelPart, &rClusterMask, number_of_entities, RequiredMinimumRedundancy](auto pContainer){
        using container_type = BareType<decltype(*pContainer)>;

        container_type temp;
        for (IndexType i_entity = 0; i_entity < number_of_entities; ++i_entity) {
            const double value = cluster_data[i_entity];
            if (value >= RequiredMinimumRedundancy) {
                temp.insert(temp.end(), *(pContainer->ptr_begin() + i_entity));
            }
        }

        if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
            rModelPart.AddElements(temp.begin(), temp.end());
            std::vector<ModelPart::NodeType::Pointer> nodes;
            for (auto& r_element : temp) {
                for (auto itr_p_node_begin = r_element.GetGeometry().ptr_begin(); itr_p_node_begin != r_element.GetGeometry().ptr_end(); ++itr_p_node_begin) {
                    nodes.push_back(*itr_p_node_begin);
                }
            }
            rModelPart.AddNodes(nodes.begin(), nodes.end());
        } else if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
            rModelPart.AddConditions(temp.begin(), temp.end());
            std::vector<ModelPart::NodeType::Pointer> nodes;
            for (auto& r_condition : temp) {
                for (auto itr_p_node_begin = r_condition.GetGeometry().ptr_begin(); itr_p_node_begin != r_condition.GetGeometry().ptr_end(); ++itr_p_node_begin) {
                    nodes.push_back(*itr_p_node_begin);
                }
            }
            rModelPart.AddNodes(nodes.begin(), nodes.end());
        } else if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
            rModelPart.AddNodes(temp.begin(), temp.end());
        } else {
            KRATOS_ERROR << "Unsupported container type in the cluster mask tensor. "
                         << "Only supports tensor adaptors with conditions or elements [ tensor adaptor = "
                         << rClusterMask << " ].\n";
        }

    }, rClusterMask.GetContainer());

    KRATOS_CATCH("");
}

std::vector<IndexType> MaskUtils::GetMasksDividingReferenceMask(
    const TensorAdaptor<double>& rReferenceMask,
    const std::vector<TensorAdaptor<double>::Pointer>& rMasksList,
    const IndexType RequiredMinimumRedundancy)
{
    KRATOS_TRY

    const auto reference_mask_coverage = GetMaskSize(rReferenceMask, RequiredMinimumRedundancy);

    std::vector<IndexType> indices;
    for (IndexType i = 0; i < rMasksList.size(); ++i) {
        const auto& p_intersected_ta = Intersect(rReferenceMask, *rMasksList[i], RequiredMinimumRedundancy);
        const auto intersection_coverage = GetMaskSize(*p_intersected_ta, RequiredMinimumRedundancy);

        if (intersection_coverage > 0 && intersection_coverage < reference_mask_coverage) {
            indices.push_back(i);
        }
    }

    return indices;

    KRATOS_CATCH("");
}

std::vector<TensorAdaptor<double>::Pointer> MaskUtils::SmoothenMasks(
    const std::vector<TensorAdaptor<double>::Pointer>& rMaskList,
    const double Radius,
    const IndexType MaxLeafSize)
{
    KRATOS_TRY

    if (rMaskList.empty()) {
        return std::vector<TensorAdaptor<double>::Pointer>{};
    }

    return std::visit([&rMaskList, Radius, MaxLeafSize](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        std::vector<TensorAdaptor<double>::Pointer> smooth_mask_list;

        if constexpr(IsInList<container_type, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {

            using position_adaptor_type = NanoFlannSingleContainerPositionAdapter<container_type>;

            using pointer_vector_type = typename position_adaptor_type::PointerVectorType;

            using distance_metric_type = typename nanoflann::metric_L2_Simple::traits<double, position_adaptor_type>::distance_t;

            using k_d_tree_index_type = nanoflann::KDTreeSingleIndexAdaptor<distance_metric_type, position_adaptor_type, 3>;

            using k_d_tree_thread_local_storage = NanoFlannKDTreeThreadLocalStorage<pointer_vector_type>;

            smooth_mask_list.reserve(rMaskList.size());

            // first construct the Kd tree for all masks
            position_adaptor_type adaptor(&*pContainer);

            k_d_tree_index_type k_d_tree_index(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(MaxLeafSize, nanoflann::KDTreeSingleIndexAdaptorFlags::None, 0));
            k_d_tree_index.buildIndex();

            const double radius_square = Radius * Radius;

            SmoothClamper clamper(0, 1);

            for (const auto& p_mask : rMaskList) {
                const auto& r_container = *std::get<typename container_type::Pointer>(p_mask->GetContainer());

                KRATOS_ERROR_IF_NOT(&r_container == &*pContainer)
                    << "Container mismatch in tensor adaptors [ reference tensor adaptor =  "
                    << *rMaskList.front() << ", other tensor adaptor = " << *p_mask << " ].\n";

                KRATOS_ERROR_IF_NOT(p_mask->Shape().size() == 1)
                    << "Mask tensor adaptors should be scalar tensor adaptors. [ tensor adaptor = " << *p_mask << " ].\n";

                const auto mask_data_view = p_mask->ViewData();

                auto p_result = p_mask->Clone();
                auto result_data_view = p_result->ViewData();

                IndexPartition<IndexType>(r_container.size()).for_each(k_d_tree_thread_local_storage(), [&k_d_tree_index, &r_container, &mask_data_view, &result_data_view, &clamper, radius_square](const IndexType Index, auto& rTLS){
                    // search for entities within radius
                    k_d_tree_index.radiusSearch(OptimizationUtils::GetEntityPosition(*(r_container.begin() + Index)).data().begin(), radius_square + std::numeric_limits<double>::epsilon(), rTLS.mNeighbourIndicesAndSquaredDistances, nanoflann::SearchParameters());

                    auto& r_value = result_data_view[Index];
                    r_value = 1.0;
                    for (const auto& r_result : rTLS.mNeighbourIndicesAndSquaredDistances) {
                        if (mask_data_view[r_result.first] == 0) {
                            const auto distance_ratio = std::sqrt(r_result.second / radius_square);
                            r_value = clamper.ProjectForward(distance_ratio);
                            break;
                        }
                    }
                });

                smooth_mask_list.push_back(p_result);
            }
        } else {
            KRATOS_ERROR << "Unsupported container type found in mask tensor adaptor = " << *rMaskList.front() << " .\n";
        }

        return std::vector<TensorAdaptor<double>::Pointer>{};
    }, rMaskList.front()->GetContainer());

    KRATOS_CATCH("");
}

} // namespace Kratos