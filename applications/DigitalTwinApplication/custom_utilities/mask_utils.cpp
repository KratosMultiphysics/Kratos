//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

// System includes
#include <numeric>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "expression/literal_flat_expression.h"

// Application includes

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
ContainerExpression<TContainerType> MaskUtils::Substract(
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

// template instantiations
#ifndef KRATOS_DT_APP_MASK_UTILS_INSTANTIATION
#define KRATOS_DT_APP_MASK_UTILS_INSTANTIATION(CONTAINER_TYPE)         \
    template void MaskUtils::CheckCompatibility(                       \
        const ContainerExpression<CONTAINER_TYPE>&,                    \
        const ContainerExpression<CONTAINER_TYPE>&);                   \
    template std::size_t MaskUtils::GetMaskSize(                       \
        const ContainerExpression<CONTAINER_TYPE>&, const IndexType);  \
    template ContainerExpression<CONTAINER_TYPE> MaskUtils::GetMask(   \
        const ContainerExpression<CONTAINER_TYPE>&);                   \
    template ContainerExpression<CONTAINER_TYPE> MaskUtils::GetMask(   \
        const ContainerExpression<CONTAINER_TYPE>&, const double);     \
    template ContainerExpression<CONTAINER_TYPE> MaskUtils::Union(     \
        const ContainerExpression<CONTAINER_TYPE>&,                    \
        const ContainerExpression<CONTAINER_TYPE>&, const IndexType);  \
    template ContainerExpression<CONTAINER_TYPE> MaskUtils::Intersect( \
        const ContainerExpression<CONTAINER_TYPE>&,                    \
        const ContainerExpression<CONTAINER_TYPE>&, const IndexType);  \
    template ContainerExpression<CONTAINER_TYPE> MaskUtils::Substract( \
        const ContainerExpression<CONTAINER_TYPE>&,                    \
        const ContainerExpression<CONTAINER_TYPE>&, const IndexType);
#endif

KRATOS_DT_APP_MASK_UTILS_INSTANTIATION(ModelPart::NodesContainerType)
KRATOS_DT_APP_MASK_UTILS_INSTANTIATION(ModelPart::ConditionsContainerType)
KRATOS_DT_APP_MASK_UTILS_INSTANTIATION(ModelPart::ElementsContainerType)

#undef KRATOS_DT_APP_MASK_UTILS_INSTANTIATION

} // namespace Kratos