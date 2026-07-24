#include "containers/multi_level_data_value_container.h"

namespace Kratos
{

MultiLevelDataValueAccesor::UniquePointer MultiLevelDataValueAccesor::Clone() const
{
    KRATOS_ERROR << "Clone method must be implemented in derived class." << std::endl;
}

std::size_t MultiLevelDataValueAccesor::GetIndex(
    const std::size_t Index1) const
{
    KRATOS_ERROR << "GetIndex method must be implemented in derived class." << std::endl;
}

std::size_t MultiLevelDataValueAccesor::GetIndex(
    const std::size_t Index1,
    const std::size_t Index2) const
{
    KRATOS_ERROR << "GetIndex method must be implemented in derived class." << std::endl;
}

std::size_t MultiLevelDataValueAccesor::GetIndex(
    const std::size_t Index1,
    const std::size_t Index2,
    const std::size_t Index3) const
{
    KRATOS_ERROR << "GetIndex method must be implemented in derived class." << std::endl;
}

std::size_t MultiLevelDataValueAccesor::size() const
{
    KRATOS_ERROR << "size method must be implemented in derived class." << std::endl;
}

MultiLevelDataValueContainer::MultiLevelDataValueContainer(
    const MultiLevelDataValueContainer& rOther)
    : mData(rOther.mData)
{}

MultiLevelDataValueContainer::MultiLevelDataValueContainer(
    MultiLevelDataValueContainer&& rOther) noexcept
    : mData(std::move(rOther.mData))
{}

MultiLevelDataValueContainer& MultiLevelDataValueContainer::operator=(
    const MultiLevelDataValueContainer& rOther)
{
    if (this != &rOther) {
        mData = rOther.mData;
    }
    return *this;
}

MultiLevelDataValueContainer& MultiLevelDataValueContainer::operator=(
    MultiLevelDataValueContainer&& rOther) noexcept
{
    if (this != &rOther) {
        mData = std::move(rOther.mData);
    }
    return *this;
}

const MultiLevelDataValueAccesor::UniquePointer& 
  MultiLevelDataValueContainer::GetIndexAccesor(
    const VariableData& rThisVariable) const
{
    auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.Key()));
    if (it_i != mData.end()) {
        return (it_i->GetIndexAccesor());
    } else {
        KRATOS_ERROR << "The variable was not found in database, please do AddVariable first." << std::endl;
    }
}

} // namespace Kratos