#include "utilities/layered_thickness_data_container.h"

namespace Kratos
{
LayeredThicknessDataContainer::LayeredThicknessDataContainer(
    const std::vector<DataValueContainer>& rContainerVector)
    : mContainerVector(rContainerVector)
{
}

DataValueContainer& LayeredThicknessDataContainer ::Get(
    const std::size_t Index)
{
    return mContainerVector.at(Index);
}

const DataValueContainer& LayeredThicknessDataContainer ::Get(
    const std::size_t Index) const
{
    return mContainerVector.at(Index);
}

// Set function to modify the DataValueContainer at a specific Index
void LayeredThicknessDataContainer ::Set(
    const std::size_t Index,
    const DataValueContainer& rValue)
{
    if (Index < mContainerVector.size()) {
        mContainerVector[Index] = rValue;
    }
    else {
        KRATOS_ERROR << "Index out of range" << std::endl;
    }
}

// Function to add a new DataValueContainer to the vector
void LayeredThicknessDataContainer ::Add(
    const DataValueContainer& rValue)
{
    mContainerVector.push_back(rValue);
}

// Function to get the size of the container vector
std::size_t LayeredThicknessDataContainer ::Size() const
{
    return mContainerVector.size();
}

}