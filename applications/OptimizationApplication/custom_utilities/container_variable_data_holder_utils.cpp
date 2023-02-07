//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <variant>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Include base h
#include "container_variable_data_holder_utils.h"

namespace Kratos
{

template<class TContainerType>
double ContainerVariableDataHolderUtils::EntityMaxNormL2(const ContainerVariableDataHolderBase<TContainerType>& rContainer)
{
    if (rContainer.GetDataDimension() == 0) {
        return 0.0;
    }

    const auto& r_data = rContainer.GetData();
    const IndexType data_dimension = rContainer.GetDataDimension();

    const IndexType number_of_entities = r_data.size() / data_dimension;

    return std::sqrt(rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(number_of_entities).for_each<MaxReduction<double>>([&](const IndexType Index) {
        double value = 0.0;
        const IndexType local_start = Index * data_dimension;
        for (IndexType i = 0; i < data_dimension; ++i) {
            value += std::pow(r_data[local_start + i], 2);
        }
        return value;
    })));
}

template<class TContainerType>
double ContainerVariableDataHolderUtils::NormInf(const ContainerVariableDataHolderBase<TContainerType>& rContainer)
{
    const auto& r_data = rContainer.GetData();
    return rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(r_data.size()).for_each<MaxReduction<double>>([&](const IndexType Index) {
        return r_data[Index];
    }));
}

double ContainerVariableDataHolderUtils::NormInf(
    const CollectiveVariableDataHolder& rContainer)
{
    double max_norm = 0.0;
    for (const auto& p_variable_data_container : rContainer.GetVariableDataHolders()) {
        std::visit([&](auto&& v) { max_norm = std::max(max_norm, NormInf(*v));}, p_variable_data_container);
    }

    return max_norm;
}

template<class TContainerType>
double ContainerVariableDataHolderUtils::NormL2(
    const ContainerVariableDataHolderBase<TContainerType>& rContainer)
{
    const double local_l2_norm_square = IndexPartition<IndexType>(rContainer.GetData().size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return rContainer.GetData()[Index];
    });

    return std::sqrt(rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(local_l2_norm_square));
}

double ContainerVariableDataHolderUtils::NormL2(
    const CollectiveVariableDataHolder& rContainer)
{
    double l2_norm_square = 0.0;
    for (const auto& p_variable_data_container : rContainer.GetVariableDataHolders()) {
        std::visit([&](auto&& v) { l2_norm_square += std::pow(NormL2(*v), 2);}, p_variable_data_container);
    }

    return std::sqrt(l2_norm_square);
}

template<class TContainerType>
double ContainerVariableDataHolderUtils::InnerProduct(
    const ContainerVariableDataHolderBase<TContainerType>& rContainer1,
    const ContainerVariableDataHolderBase<TContainerType>& rContainer2)
{
    const auto& r_data_1 = rContainer1.GetData();
    const auto& r_data_2 = rContainer2.GetData();

    KRATOS_ERROR_IF(r_data_1.size() != r_data_2.size())
        << "Data size mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    KRATOS_ERROR_IF(rContainer1.GetModelPart() != rContainer2.GetModelPart())
        << "Model part mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    return rContainer1.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(IndexPartition<IndexType>(r_data_1.size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return r_data_1[Index] * r_data_2[Index];
    }));
}

double ContainerVariableDataHolderUtils::InnerProduct(
    const CollectiveVariableDataHolder& rContainer1,
    const CollectiveVariableDataHolder& rContainer2)
{
    KRATOS_ERROR_IF_NOT(rContainer1.IsCompatibleWith(rContainer2))
        << "Unsupported collective variable data holders provided for \"+\" operation."
        << "\nLeft operand : " << rContainer1 << "\nRight operand: " << rContainer2 << std::endl;

    double inner_product_value = 0.0;
    for (IndexType i = 0; i < rContainer1.GetVariableDataHolders().size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            inner_product_value += InnerProduct(*v, *std::get<v_type>(rContainer2.GetVariableDataHolders()[i]));
        }, rContainer1.GetVariableDataHolders()[i]);
    }
    return inner_product_value;
}

// template instantiations
#define INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ContainerType)                                                                                                                \
    template double ContainerVariableDataHolderUtils::EntityMaxNormL2(const ContainerVariableDataHolderBase<ContainerType>&);                                                       \
    template double ContainerVariableDataHolderUtils::NormInf(const ContainerVariableDataHolderBase<ContainerType>&);                                                               \
    template double ContainerVariableDataHolderUtils::NormL2(const ContainerVariableDataHolderBase<ContainerType>&);                                                                \
    template double ContainerVariableDataHolderUtils::InnerProduct(const ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ContainerType>&);

INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType)

#undef INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE

}