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

// template instantiations
#define INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ContainerType)                                                                                                                \
    template double ContainerVariableDataHolderUtils::EntityMaxNormL2(const ContainerVariableDataHolderBase<ContainerType>&);                                                       \
    template double ContainerVariableDataHolderUtils::NormInf(const ContainerVariableDataHolderBase<ContainerType>&);                                                               \
    template double ContainerVariableDataHolderUtils::InnerProduct(const ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ContainerType>&);

INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType)

#undef INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE

}