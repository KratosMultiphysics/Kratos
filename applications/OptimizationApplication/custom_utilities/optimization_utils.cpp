//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes

// Include base h
#include "optimization_utils.h"

namespace Kratos
{
template<>
IndexType OptimizationUtils::GetLocalSize<double>(const IndexType DomainSize)
{
    return 1;
}

template<>
IndexType OptimizationUtils::GetLocalSize<array_1d<double, 3>>(const IndexType DomainSize)
{
    return DomainSize;
}

template<>
void OptimizationUtils::AssignValue(
    const double& rValue,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    Vector& rOutput)
{
    rOutput[VectoStartingIndex] = rValue;
}

template<>
void OptimizationUtils::AssignValue(
    const array_1d<double, 3>& rValue,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    Vector& rOutput)
{
    rOutput[VectoStartingIndex + ValueComponentIndex] = rValue[ValueComponentIndex];
}

template<class TContainerType>
void OptimizationUtils::GetContainerIds(
    const TContainerType& rContainer,
    std::vector<IndexType>& rOutput)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    rOutput.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        rOutput[Index] = (rContainer.begin() + Index)->Id();
    });

    KRATOS_CATCH("")
}

template<class TContainerType, class TDataType>
void OptimizationUtils::GetContainerVariableToVector(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const IndexType DomainSize,
    Vector& rOutput)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();
    const IndexType local_size = GetLocalSize<TDataType>(DomainSize);

    rOutput.resize(number_of_entities * local_size);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        const auto& values = (rContainer.begin() + Index)->GetValue(rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            AssignValue(values, i, Index * local_size, rOutput);
        }
    });

    KRATOS_CATCH("")
}

template<class TContainerType>
void OptimizationUtils::GetContainerPropertiesVariableToVector(
    const TContainerType& rContainer,
    const Variable<double> & rVariable,
    Vector& rOutput)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    rOutput.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        rOutput[Index] = (rContainer.begin() + Index)->GetProperties()[rVariable];
    });

    KRATOS_CATCH("")
}

template<class TContainerType>
GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    if (rContainer.size() > 0) {
        const auto first_geometry_type = rContainer.begin()->GetGeometry().GetGeometryType();
        const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity) -> IndexType {
            return first_geometry_type == rEntity.GetGeometry().GetGeometryType();
        });

        if (rDataCommunicator.AndReduceAll(local_value)) {
            return first_geometry_type;
        } else {
            return GeometryData::KratosGeometryType::Kratos_generic_type;
        }
    } else {
        return GeometryData::KratosGeometryType::Kratos_generic_type;
    }

    KRATOS_CATCH("")
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAllContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.AndReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MaxReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.OrReduceAll(local_value);

    KRATOS_CATCH("");
}

// template instantiations
template void OptimizationUtils::GetContainerIds(const ModelPart::NodesContainerType&, std::vector<IndexType>&);
template void OptimizationUtils::GetContainerIds(const ModelPart::ConditionsContainerType&, std::vector<IndexType>&);
template void OptimizationUtils::GetContainerIds(const ModelPart::ElementsContainerType&, std::vector<IndexType>&);

template void OptimizationUtils::GetContainerVariableToVector(const ModelPart::NodesContainerType&, const Variable<double>&, const IndexType, Vector&);
template void OptimizationUtils::GetContainerVariableToVector(const ModelPart::ConditionsContainerType&, const Variable<double>&, const IndexType, Vector&);
template void OptimizationUtils::GetContainerVariableToVector(const ModelPart::ElementsContainerType&, const Variable<double>&, const IndexType, Vector&);

template void OptimizationUtils::GetContainerVariableToVector(const ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, Vector&);
template void OptimizationUtils::GetContainerVariableToVector(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, Vector&);
template void OptimizationUtils::GetContainerVariableToVector(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, Vector&);

template GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template void OptimizationUtils::GetContainerPropertiesVariableToVector(const ModelPart::ConditionsContainerType&, const Variable<double>&, Vector&);
template void OptimizationUtils::GetContainerPropertiesVariableToVector(const ModelPart::ElementsContainerType&, const Variable<double>&, Vector&);

template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

}