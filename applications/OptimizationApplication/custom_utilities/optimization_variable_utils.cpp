//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Reza Najian Asl, https://github.com/RezaNajian
//                   Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes

// Include base h
#include "optimization_variable_utils.h"

namespace Kratos
{

template<class TContainerType>
void OptimizationVariableUtils::GetContainerIds(
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

template<class TContainerType>
void OptimizationVariableUtils::GetContainerVariableToMatrix(
    const TContainerType& rContainer,
    const Variable<array_1d<double, 3>> & rVariable,
    Matrix& rOutput)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    rOutput.resize(number_of_entities, 3);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        row(rOutput, Index) = (rContainer.begin() + Index)->GetValue(rVariable);
    });

    KRATOS_CATCH("")
}

template<class TContainerType>
void OptimizationVariableUtils::GetContainerVariableToVector(
    const TContainerType& rContainer,
    const Variable<double> & rVariable,
    Vector& rOutput)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    rOutput.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        rOutput[Index] = (rContainer.begin() + Index)->GetValue(rVariable);
    });

    KRATOS_CATCH("")
}

template<class TContainerType>
void OptimizationVariableUtils::GetContainerPropertiesVariableToVector(
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
GeometryData::KratosGeometryType OptimizationVariableUtils::GetContainerEntityGeometryType(
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
bool OptimizationVariableUtils::IsVariableExistsInAllContainerProperties(
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
bool OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties(
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
template void OptimizationVariableUtils::GetContainerIds(const ModelPart::NodesContainerType&, std::vector<IndexType>&);
template void OptimizationVariableUtils::GetContainerIds(const ModelPart::ConditionsContainerType&, std::vector<IndexType>&);
template void OptimizationVariableUtils::GetContainerIds(const ModelPart::ElementsContainerType&, std::vector<IndexType>&);

template void OptimizationVariableUtils::GetContainerVariableToMatrix(const ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, Matrix&);
template void OptimizationVariableUtils::GetContainerVariableToMatrix(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, Matrix&);
template void OptimizationVariableUtils::GetContainerVariableToMatrix(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, Matrix&);

template void OptimizationVariableUtils::GetContainerVariableToVector(const ModelPart::NodesContainerType&, const Variable<double>&, Vector&);
template void OptimizationVariableUtils::GetContainerVariableToVector(const ModelPart::ConditionsContainerType&, const Variable<double>&, Vector&);
template void OptimizationVariableUtils::GetContainerVariableToVector(const ModelPart::ElementsContainerType&, const Variable<double>&, Vector&);

template GeometryData::KratosGeometryType OptimizationVariableUtils::GetContainerEntityGeometryType(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template GeometryData::KratosGeometryType OptimizationVariableUtils::GetContainerEntityGeometryType(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template void OptimizationVariableUtils::GetContainerPropertiesVariableToVector(const ModelPart::ConditionsContainerType&, const Variable<double>&, Vector&);
template void OptimizationVariableUtils::GetContainerPropertiesVariableToVector(const ModelPart::ElementsContainerType&, const Variable<double>&, Vector&);

template bool OptimizationVariableUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationVariableUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationVariableUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationVariableUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template bool OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationVariableUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

}