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
#include <sstream>
#include <limits>

// Project includes
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "container_data.h"

namespace Kratos {

namespace ContainerDataHelperUtilities
{

void AssignValueToVector(
    Vector& rOutput,
    const IndexType StartIndex,
    const IndexType ComponentIndex,
    const double& rValue)
{
    rOutput[StartIndex] = rValue;
}

void AssignValueToVector(
    Vector& rOutput,
    const IndexType StartIndex,
    const IndexType ComponentIndex,
    const array_1d<double, 3>& rValue)
{
    rOutput[StartIndex + ComponentIndex] = rValue[ComponentIndex];
}
}

ContainerData::ContainerData(
    ModelPart& rModelPart,
    const ContainerDataType& rContainerDataType)
    : mrModelPart(rModelPart),
      mrContainerDataType(rContainerDataType)
{
}

ContainerData::ContainerData(const ContainerData& rOther)
    : mrModelPart(rOther.mrModelPart),
      mrContainerDataType(rOther.mrContainerDataType)
{
    this->mData.resize(rOther.mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index){
        this->mData[Index] = rOther.mData[Index];
    });
}

template<class TDataType>
void ContainerData::AssignDataToContainerVariable(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << mrModelPart.FullName() << ".\n";

    switch (mrContainerDataType) {
        case ContainerDataType::NodalHistorical:
            OptimizationUtils::AssignVectorToContainerHistoricalVariable(mrModelPart, rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE], mData);
            break;
        case ContainerDataType::NodalNonHistorical:
            OptimizationUtils::AssignVectorToContainerVariable(mrModelPart.Nodes(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE], mData);
            break;
        case ContainerDataType::ConditionNonHistorical:
            OptimizationUtils::AssignVectorToContainerVariable(mrModelPart.Conditions(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE], mData);
            break;
        case ContainerDataType::ConditionProperties:
            OptimizationUtils::AssignVectorToContainerPropertiesVariable(mrModelPart.Conditions(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE], mData);
            break;
        case ContainerDataType::ElementNonHistorical:
            OptimizationUtils::AssignVectorToContainerVariable(mrModelPart.Elements(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE], mData);
            break;
        case ContainerDataType::ElementProperties:
            OptimizationUtils::AssignVectorToContainerPropertiesVariable(mrModelPart.Elements(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE], mData);
            break;
    }

    KRATOS_CATCH("");
}

template<class TDataType>
void ContainerData::ReadDataFromContainerVariable(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << mrModelPart.FullName() << ".\n";

    switch (mrContainerDataType) {
        case ContainerDataType::NodalHistorical:
            OptimizationUtils::GetHistoricalContainerVariableToVector(mData, mrModelPart, rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);
            break;
        case ContainerDataType::NodalNonHistorical:
            OptimizationUtils::GetContainerVariableToVector(mData, mrModelPart.Nodes(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);
            break;
        case ContainerDataType::ConditionNonHistorical:
            OptimizationUtils::GetContainerVariableToVector(mData, mrModelPart.Conditions(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);
            break;
        case ContainerDataType::ConditionProperties:
            OptimizationUtils::GetContainerPropertiesVariableToVector(mData, mrModelPart.Conditions(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);
            break;
        case ContainerDataType::ElementNonHistorical:
            OptimizationUtils::GetContainerVariableToVector(mData, mrModelPart.Elements(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);
            break;
        case ContainerDataType::ElementProperties:
            OptimizationUtils::GetContainerPropertiesVariableToVector(mData, mrModelPart.Elements(), rVariable, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);
            break;
    }

    KRATOS_CATCH("");
}

template<class TDataType>
void ContainerData::SetDataForContainerVariable(const Variable<TDataType>& rVariable, const TDataType& rValue)
{

    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << mrModelPart.FullName() << ".\n";

    IndexType number_of_entities = 0;
    const IndexType local_size = OptimizationUtils::GetLocalSize<TDataType>(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);

    switch (mrContainerDataType) {
        case ContainerDataType::NodalHistorical:
        case ContainerDataType::NodalNonHistorical:
            number_of_entities = mrModelPart.NumberOfNodes();
            break;
        case ContainerDataType::ConditionNonHistorical:
        case ContainerDataType::ConditionProperties:
            number_of_entities = mrModelPart.NumberOfConditions();
            break;
        case ContainerDataType::ElementNonHistorical:
        case ContainerDataType::ElementProperties:
            number_of_entities = mrModelPart.NumberOfElements();
            break;
    }

    if (this->mData.size() != local_size * number_of_entities) {
        this->mData.resize(local_size * number_of_entities, false);
    }

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index) {
        for (IndexType i = 0; i < local_size; ++i) {
            ContainerDataHelperUtilities::AssignValueToVector(this->mData, Index * local_size, i, rValue);
        }
    });
}

ContainerData ContainerData::Clone()
{
    return ContainerData(*this);
}

bool ContainerData::IsSameContainer(const ContainerData& rOther) const
{
    return this->mrModelPart == rOther.mrModelPart && this->mrContainerDataType == rOther.mrContainerDataType;
}

double ContainerData::NormInf() const
{
    return mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(this->mData.size()).for_each<MaxReduction<double>>([&](const IndexType Index) {
        return this->mData[Index];
    }));
}

double ContainerData::InnerProduct(const ContainerData& rOther) const
{
    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->IsSameContainer(rOther))
        << "Data containers are of not same type.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    // TODO: Make it MPI compatible with SumAll, but need to take care of ghost nodes in the case nodal historical/non historical
    //       containers are used
    return IndexPartition<IndexType>(this->mData.size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return this->mData[Index] * rOther.mData[Index];
    });
}

Vector& ContainerData::GetData()
{
    return mData;
}

const Vector& ContainerData::GetData() const
{
    return mData;
}

const ContainerData::ContainerDataType& ContainerData::GetContainerDataType() const
{
    return mrContainerDataType;
}

const ModelPart& ContainerData::GetModelPart() const
{
    return mrModelPart;
}

ModelPart& ContainerData::GetModelPart()
{
    return mrModelPart;
}

std::string ContainerData::Info() const
{
    std::stringstream msg;
    msg << "ContainerData [ ModelPartName: " << mrModelPart.FullName() << ", ContainerType: ";
    switch (mrContainerDataType) {
        case ContainerDataType::NodalHistorical:
            msg << "NodalHistorical";
            break;
        case ContainerDataType::NodalNonHistorical:
            msg << "NodalNonHistorical";
            break;
        case ContainerDataType::ConditionNonHistorical:
            msg << "ConditionNonHistorical";
            break;
        case ContainerDataType::ConditionProperties:
            msg << "ConditionProperties";
            break;
        case ContainerDataType::ElementNonHistorical:
            msg << "ElementNonHistorical";
            break;
        case ContainerDataType::ElementProperties:
            msg << "ElementProperties";
            break;
    }
    msg << ", DataSize = " << mData.size() << " ]";
    return msg.str();
}

void ContainerData::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

std::variant<
    std::reference_wrapper<ModelPart::NodesContainerType>,
    std::reference_wrapper<ModelPart::ConditionsContainerType>,
    std::reference_wrapper<ModelPart::ElementsContainerType>> ContainerData::GetContainer()
{
    switch (mrContainerDataType) {
        case ContainerDataType::NodalHistorical:
        case ContainerDataType::NodalNonHistorical:
            return mrModelPart.Nodes();
        case ContainerDataType::ConditionNonHistorical:
        case ContainerDataType::ConditionProperties:
            return mrModelPart.Conditions();
        case ContainerDataType::ElementNonHistorical:
        case ContainerDataType::ElementProperties:
            return mrModelPart.Elements();
    }

    KRATOS_ERROR << "Unsupported container data type.";
    return mrModelPart.Nodes();
}

std::variant<
    std::reference_wrapper<const ModelPart::NodesContainerType>,
    std::reference_wrapper<const ModelPart::ConditionsContainerType>,
    std::reference_wrapper<const ModelPart::ElementsContainerType>> ContainerData::GetContainer() const
{
    switch (mrContainerDataType) {
        case ContainerDataType::NodalHistorical:
        case ContainerDataType::NodalNonHistorical:
            return mrModelPart.Nodes();
        case ContainerDataType::ConditionNonHistorical:
        case ContainerDataType::ConditionProperties:
            return mrModelPart.Conditions();
        case ContainerDataType::ElementNonHistorical:
        case ContainerDataType::ElementProperties:
            return mrModelPart.Elements();
    }

    KRATOS_ERROR << "Unsupported container data type.";
    return mrModelPart.Nodes();
}

ContainerData ContainerData::operator+(const ContainerData& rOther) const
{
    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->IsSameContainer(rOther))
        << "Data containers are of not same type.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerData result(this->mrModelPart, this->mrContainerDataType);
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] + rOther.mData[Index];
    });


    return result;
}

ContainerData& ContainerData::operator+=(const ContainerData& rOther)
{
    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +=.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->IsSameContainer(rOther))
        << "Data containers are of not same type.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] += rOther.mData[Index];
    });

    return *this;
}

ContainerData ContainerData::operator-(const ContainerData& rOther) const
{
    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for -.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->IsSameContainer(rOther))
        << "Data containers are of not same type.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerData result(this->mrModelPart, this->mrContainerDataType);
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] - rOther.mData[Index];
    });

    return result;
}

ContainerData& ContainerData::operator-=(const ContainerData& rOther)
{
    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for -=.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->IsSameContainer(rOther))
        << "Data containers are of not same type.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] -= rOther.mData[Index];
    });

    return *this;
}

ContainerData ContainerData::operator*(const double Value) const
{
    ContainerData result(this->mrModelPart, this->mrContainerDataType);
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] * Value;
    });

    return result;
}

ContainerData& ContainerData::operator*=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] *= Value;
    });

    return *this;
}

ContainerData ContainerData::operator/(const double Value) const
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      divisor           : " << Value << "\n";

    ContainerData result(this->mrModelPart, this->mrContainerDataType);
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] / Value;
    });

    return result;
}

ContainerData& ContainerData::operator/=(const double Value)
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      divisor           : " << Value << "\n";

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] /= Value;
    });

    return *this;
}

ContainerData& ContainerData::operator=(const ContainerData& rOther)
{
    KRATOS_ERROR_IF_NOT(this->IsSameContainer(rOther))
        << "Data containers are of not same type.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    if (this->mData.size() != rOther.mData.size()) {
        this->mData.resize(rOther.mData.size(), false);
    }

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] = rOther.mData[Index];
    });

    return *this;
}

// template instantiations

template void ContainerData::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerData::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);

template void ContainerData::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerData::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);

template void ContainerData::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerData::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

} // namespace Kratos