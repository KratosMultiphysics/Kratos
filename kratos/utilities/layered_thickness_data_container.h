#pragma once

#include "containers/data_value_container.h"

namespace Kratos
{

/**
 * @brief This class is used to store data in a layered object,
 * where each layer is represented by a DataValueContainer object.
 */
class KRATOS_API(KRATOS_CORE) LayeredThicknessDataContainer 
{
public:

    /**
     * @brief Pointer definition of LayeredThicknessDataContainer .
     */
    KRATOS_CLASS_POINTER_DEFINITION(DataValueContainer);

    /**
     * @brief Constructor with std::vector<DataValueContainer>.
     * @param rContainerVector Reference to a vector of DataValueContainer objects.
     */
    LayeredThicknessDataContainer(
        const std::vector<DataValueContainer>& rContainerVector);

    /**
     * @brief Default constructor.
     */
    LayeredThicknessDataContainer() = default;

    /**
     * @brief Copy constructor.
     * @param rOther Reference to another LayeredThicknessDataContainer  object.
     */
    LayeredThicknessDataContainer(
        const LayeredThicknessDataContainer & rOther) = default;

    /**
     * @brief Move constructor.
     * @param rOther Rvalue reference to another LayeredThicknessDataContainer  object.
     */
    LayeredThicknessDataContainer(
        LayeredThicknessDataContainer && rOther) noexcept = default;

    /**
     * @brief Copy assignment operator.
     * @param rOther Reference to another LayeredThicknessDataContainer  object.
     * @return Reference to the current object.
     */
    LayeredThicknessDataContainer & operator=(
        const LayeredThicknessDataContainer & rOther) = default;

    /**
     * @brief Move assignment operator.
     * @param rOther Rvalue reference to another LayeredThicknessDataContainer  object.
     * @return Reference to the current object.
     */
    LayeredThicknessDataContainer & operator=(
        LayeredThicknessDataContainer && rOther) noexcept = default;

    /**
     * @brief Destructor.
     */
    virtual ~LayeredThicknessDataContainer () = default;

    /**
     * @brief Get function to access the DataValueContainer at a specific index.
     * @param Index The index of the DataValueContainer to access.
     * @return Reference to the DataValueContainer at the specified index.
     */
    DataValueContainer& Get(
        const std::size_t Index);

    /**
     * @brief Const version of the Get function.
     * @param Index The index of the DataValueContainer to access.
     * @return Const reference to the DataValueContainer at the specified index.
     */
    const DataValueContainer& Get(
        const std::size_t Index) const;

    /**
     * @brief Set function to modify the DataValueContainer at a specific index.
     * @param Index The index of the DataValueContainer to modify.
     * @param rValue Reference to the new DataValueContainer value.
     */
    void Set(
        const std::size_t Index,
        const DataValueContainer& rValue);

    /**
     * @brief Function to add a new DataValueContainer to the vector.
     * @param rValue Reference to the DataValueContainer to add.
     */
    void Add(
        const DataValueContainer& rValue);

    /**
     * @brief Function to get the size of the container vector.
     * @return The size of the container vector.
     */
    std::size_t Size() const;

    // The following methods must exist to compile but are not used
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_ERROR << "Not implemented" << std::endl;
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_ERROR << "Not implemented" << std::endl;
    }

protected:
    /**
     * @brief The vector of DataValueContainer objects.
     */
    std::vector<DataValueContainer> mContainerVector;

};

inline std::istream& operator >> (
    std::istream& rIStream,
    LayeredThicknessDataContainer & rThis)
{
    return rIStream;
}

inline std::ostream& operator << (
    std::ostream& rOStream,
    const LayeredThicknessDataContainer & rThis)
{
    rOStream << "LayeredThicknessDataContainer " << std::endl;

    return rOStream;
}

} // namespace Kratos