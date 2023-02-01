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
#include <cmath>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

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
void OptimizationUtils::AssignValueToVector(
    const double& rValue,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    Vector& rOutput)
{
    rOutput[VectoStartingIndex] = rValue;
}

template<>
void OptimizationUtils::AssignValueToVector(
    const array_1d<double, 3>& rValue,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    Vector& rOutput)
{
    rOutput[VectoStartingIndex + ValueComponentIndex] = rValue[ValueComponentIndex];
}

template<>
void OptimizationUtils::AssignValueFromVector(
    double& rOutput,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    const Vector& rInput)
{
    rOutput = rInput[VectoStartingIndex];
}

template<>
void OptimizationUtils::AssignValueFromVector(
    array_1d<double, 3>& rOutput,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    const Vector& rInput)
{
    rOutput[ValueComponentIndex] = rInput[VectoStartingIndex + ValueComponentIndex];
}

// void BroadcastVector(
//     Vector& rSourceValue,
//     const DataCommunicator& rCommunicator,
//     const IndexType SourceRank)
// {
//     if (!rCommunicator.IsDistributed()) {
//         return;
//     } else {
//         rCommunicator.Broadcast(rSourceValue, SourceRank);
//     }
// }

// void BroadcastMatrix(
//     Matrix& rSourceValue,
//     const DataCommunicator& rCommunicator,
//     const IndexType SourceRank)
// {
//     if (!rCommunicator.IsDistributed()) {
//         return;
//     } else {
//         rCommunicator.Broadcast(rSourceValue, SourceRank);
//     }
// }

template<class TContainerType>
void OptimizationUtils::GetContainerIds(
    std::vector<IndexType>& rOutput,
    const TContainerType& rContainer)
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
    Vector& rOutput,
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const IndexType DomainSize)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();
    const IndexType local_size = GetLocalSize<TDataType>(DomainSize);

    rOutput.resize(number_of_entities * local_size);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        const auto& values = (rContainer.begin() + Index)->GetValue(rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            AssignValueToVector(values, i, Index * local_size, rOutput);
        }
    });

    KRATOS_CATCH("")
}

template<class TDataType>
void OptimizationUtils::GetHistoricalContainerVariableToVector(
    Vector& rOutput,
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const IndexType DomainSize)
{
    KRATOS_TRY

    const IndexType number_of_entities = rModelPart.NumberOfNodes();
    const IndexType local_size = GetLocalSize<TDataType>(DomainSize);

    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable))
        << rVariable.Name() << " is not found in the nodal solution step variables container of "
        << rModelPart.FullName() << ".\n";

    rOutput.resize(number_of_entities * local_size);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        const auto& values = (rModelPart.NodesBegin() + Index)->FastGetSolutionStepValue(rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            AssignValueToVector(values, i, Index * local_size, rOutput);
        }
    });

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
void OptimizationUtils::GetContainerPropertiesVariableToVector(
    Vector& rOutput,
    const TContainerType& rContainer,
    const Variable<TDataType> & rVariable,
    const IndexType DomainSize)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();
    const IndexType local_size = GetLocalSize<TDataType>(DomainSize);

    rOutput.resize(number_of_entities * local_size);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        const auto& values = (rContainer.begin() + Index)->GetProperties().GetValue(rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            AssignValueToVector(values, i, Index * local_size, rOutput);
        }
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

double OptimizationUtils::CalculateVectorL2Norm(const Vector& rInput)
{
    KRATOS_TRY

    const auto square_sum = IndexPartition<IndexType>(rInput.size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return rInput[0] * rInput[0];
    });

    return std::sqrt(square_sum);

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
void OptimizationUtils::AssignVectorToContainerPropertiesVariable(
    TContainerType& rContainer,
    const Variable<TDataType>& rPropertiesVariable,
    const IndexType DomainSize,
    const Vector& rValues)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();
    const IndexType local_size = GetLocalSize<TDataType>(DomainSize);

    KRATOS_ERROR_IF_NOT(number_of_entities * local_size == rValues.size())
        << "Container size and values size mismatch. [ Values size = "
        << rValues.size() << ", container size = " << number_of_entities
        << ", local size = " << local_size << " ].\n";

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        TDataType value;
        for (IndexType i = 0; i < local_size; ++i) {
            AssignValueFromVector(value, i, Index * local_size, rValues);
        }
        (rContainer.begin() + Index)->GetProperties().SetValue(rPropertiesVariable, value);
    });

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
void OptimizationUtils::AssignVectorToContainerVariable(
    TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const IndexType DomainSize,
    const Vector& rValues)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();
    const IndexType local_size = GetLocalSize<TDataType>(DomainSize);

    KRATOS_ERROR_IF_NOT(number_of_entities * local_size == rValues.size())
        << "Container size and values size mismatch. [ Values size = "
        << rValues.size() << ", container size = " << number_of_entities
        << ", local size = " << local_size << " ].\n";

    VariableUtils().SetNonHistoricalVariableToZero(rVariable, rContainer);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        auto& value = (rContainer.begin() + Index)->GetValue(rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            AssignValueFromVector(value, i, Index * local_size, rValues);
        }
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(
    ModelPart& rModelPart,
    TContainerType& rContainer)
{
    KRATOS_TRY

    auto properties_id = block_for_each<MaxReduction<IndexType>>(rModelPart.PropertiesArray(), [](const auto pProperties) {
        return pProperties->Id();
    });

    // creation of properties is done in serial
    for (auto& r_entity : rContainer) {
        auto p_properties = rModelPart.CreateNewProperties(++properties_id);
        const auto& element_properties = r_entity.GetProperties();
        *p_properties = element_properties;
        r_entity.SetProperties(p_properties);
    }

    KRATOS_CATCH("");
}

template<class TDataType>
void OptimizationUtils::AssignVectorToContainerHistoricalVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const IndexType DomainSize,
    const Vector& rValues)
{
    KRATOS_TRY

    const IndexType number_of_nodes = rModelPart.NumberOfNodes();
    const IndexType local_size = GetLocalSize<TDataType>(DomainSize);

    KRATOS_ERROR_IF_NOT(number_of_nodes * local_size == rValues.size())
        << "Number of nodes and values size mismatch. [ Values size = "
        << rValues.size() << ", number of nodes = " << number_of_nodes
        << ", local size = " << local_size << " ].\n";

    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable))
        << rVariable.Name() << " is not found in the nodal solution step variables container of "
        << rModelPart.FullName() << ".\n";

    IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType Index){
        auto& value = (rModelPart.NodesBegin() + Index)->FastGetSolutionStepValue(rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            AssignValueFromVector(value, i, Index * local_size, rValues);
        }
    });

    KRATOS_CATCH("");
}

void OptimizationUtils::AddVectors(
    Vector& rOutput,
    const Vector& rA,
    const Vector& rB)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rA.size() == rB.size())
        << "Vector size mismatch. [ rA.size() = " << rA.size()
        << ", rB.size() = " << rB.size() << " ].\n";

    if (rOutput.size() != rA.size()) {
        rOutput.resize(rA.size(), false);
    }

    IndexPartition<IndexType>(rOutput.size()).for_each([&](const IndexType Index) {
        rOutput[Index] = rA[Index] + rB[Index];
    });

    KRATOS_CATCH("");
}

void OptimizationUtils::SubstractVectors(
    Vector& rOutput,
    const Vector& rA,
    const Vector& rB)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rA.size() == rB.size())
        << "Vector size mismatch. [ rA.size() = " << rA.size()
        << ", rB.size() = " << rB.size() << " ].\n";

    if (rOutput.size() != rA.size()) {
        rOutput.resize(rA.size(), false);
    }

    IndexPartition<IndexType>(rOutput.size()).for_each([&](const IndexType Index) {
        rOutput[Index] = rA[Index] - rB[Index];
    });

    KRATOS_CATCH("");
}

void OptimizationUtils::MultiplyVector(
    Vector& rOutput,
    const Vector& rA,
    const double Multiplier)
{
    KRATOS_TRY

    if (rOutput.size() != rA.size()) {
        rOutput.resize(rA.size(), false);
    }

    IndexPartition<IndexType>(rOutput.size()).for_each([&](const IndexType Index) {
        rOutput[Index] = rA[Index] * Multiplier;
    });

    KRATOS_CATCH("");
}

void OptimizationUtils::DivideVector(
    Vector& rOutput,
    const Vector& rA,
    const double Divisor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(Divisor == 0.0) << "Division by zero.\n";

    if (rOutput.size() != rA.size()) {
        rOutput.resize(rA.size(), false);
    }

    IndexPartition<IndexType>(rOutput.size()).for_each([&](const IndexType Index) {
        rOutput[Index] = rA[Index] / Divisor;
    });

    KRATOS_CATCH("");
}

double OptimizationUtils::NormInf(
    const Vector& rOutput)
{
    KRATOS_TRY

    return IndexPartition<IndexType>(rOutput.size()).for_each<MaxReduction<double>>([&](const IndexType Index) {
        return rOutput[Index];
    });

    KRATOS_CATCH("")
}

// template instantiations
template void OptimizationUtils::GetContainerIds(std::vector<IndexType>&, const ModelPart::NodesContainerType&);
template void OptimizationUtils::GetContainerIds(std::vector<IndexType>&, const ModelPart::ConditionsContainerType&);
template void OptimizationUtils::GetContainerIds(std::vector<IndexType>&, const ModelPart::ElementsContainerType&);

template void OptimizationUtils::GetHistoricalContainerVariableToVector(Vector&, const ModelPart&, const Variable<double>&, const IndexType);
template void OptimizationUtils::GetHistoricalContainerVariableToVector(Vector&, const ModelPart&, const Variable<array_1d<double, 3>>&, const IndexType);

template void OptimizationUtils::GetContainerVariableToVector(Vector&, const ModelPart::NodesContainerType&, const Variable<double>&, const IndexType);
template void OptimizationUtils::GetContainerVariableToVector(Vector&, const ModelPart::ConditionsContainerType&, const Variable<double>&, const IndexType);
template void OptimizationUtils::GetContainerVariableToVector(Vector&, const ModelPart::ElementsContainerType&, const Variable<double>&, const IndexType);

template void OptimizationUtils::GetContainerVariableToVector(Vector&, const ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const IndexType);
template void OptimizationUtils::GetContainerVariableToVector(Vector&, const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType);
template void OptimizationUtils::GetContainerVariableToVector(Vector&, const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType);

template void OptimizationUtils::GetContainerPropertiesVariableToVector(Vector&, const ModelPart::ConditionsContainerType&, const Variable<double>&, const IndexType);
template void OptimizationUtils::GetContainerPropertiesVariableToVector(Vector&, const ModelPart::ElementsContainerType&, const Variable<double>&, const IndexType);

template void OptimizationUtils::GetContainerPropertiesVariableToVector(Vector&, const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType);
template void OptimizationUtils::GetContainerPropertiesVariableToVector(Vector&, const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType);

template GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ConditionsContainerType&);
template void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ElementsContainerType&);

template void OptimizationUtils::AssignVectorToContainerVariable(ModelPart::NodesContainerType&, const Variable<double>&, const IndexType, const Vector&);
template void OptimizationUtils::AssignVectorToContainerVariable(ModelPart::ConditionsContainerType&, const Variable<double>&, const IndexType, const Vector&);
template void OptimizationUtils::AssignVectorToContainerVariable(ModelPart::ElementsContainerType&, const Variable<double>&, const IndexType, const Vector&);

template void OptimizationUtils::AssignVectorToContainerVariable(ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, const Vector&);
template void OptimizationUtils::AssignVectorToContainerVariable(ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, const Vector&);
template void OptimizationUtils::AssignVectorToContainerVariable(ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, const Vector&);

template void OptimizationUtils::AssignVectorToContainerPropertiesVariable(ModelPart::ConditionsContainerType&, const Variable<double>&, const IndexType, const Vector&);
template void OptimizationUtils::AssignVectorToContainerPropertiesVariable(ModelPart::ElementsContainerType&, const Variable<double>&, const IndexType, const Vector&);

template void OptimizationUtils::AssignVectorToContainerPropertiesVariable(ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, const Vector&);
template void OptimizationUtils::AssignVectorToContainerPropertiesVariable(ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const IndexType, const Vector&);

template void OptimizationUtils::AssignVectorToContainerHistoricalVariable(ModelPart&, const Variable<double>&, const IndexType, const Vector&);
template void OptimizationUtils::AssignVectorToContainerHistoricalVariable(ModelPart&, const Variable<array_1d<double, 3>>&, const IndexType, const Vector&);

template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

}