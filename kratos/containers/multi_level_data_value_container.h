#pragma once

// System includes
#include <vector>

// Kratos includes
#include "includes/kratos_components.h"
#include "containers/variable.h"

namespace Kratos 
{

class KRATOS_API(KRATOS_CORE) MultiLevelDataValueAccesor
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiLevelDataValueAccesor);

    virtual std::size_t GetIndex(
        const std::size_t Index1) const;

    virtual std::size_t GetIndex(
        const std::size_t Index1,
        const std::size_t Index2) const;

    virtual std::size_t GetIndex(
        const std::size_t Index1,
        const std::size_t Index2,
        const std::size_t Index3) const;

    virtual std::size_t size() const;

    virtual MultiLevelDataValueAccesor::UniquePointer Clone() const;

    virtual ~MultiLevelDataValueAccesor() = default;

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {}

    virtual void load(Serializer& rSerializer)
    {}
};

/**
 * @class MultiLevelDataValueContainer
 * @brief A container class for managing multi-level data values.
 * 
 * This class provides functionalities to store, retrieve, and manage data values
 * associated with different variables at multiple levels. It supports adding variables,
 * checking for their existence, and getting/setting their values using index accessors.
 * 
 * @details The container uses a nested DataBlock class to manage individual blocks of data
 * associated with specific variables. Each DataBlock is responsible for allocating memory,
 * storing values, and providing access to these values based on indices.
 * 
 * The container supports copy and move semantics, allowing for efficient copying and moving
 * of data between instances.
 */
class KRATOS_API(KRATOS_CORE) MultiLevelDataValueContainer
{
public:


    KRATOS_CLASS_POINTER_DEFINITION(MultiLevelDataValueContainer);

    MultiLevelDataValueContainer() = default;

    MultiLevelDataValueContainer(const MultiLevelDataValueContainer& rOther);

    MultiLevelDataValueContainer(MultiLevelDataValueContainer&& rOther) noexcept;

    MultiLevelDataValueContainer& operator=(const MultiLevelDataValueContainer& rOther);

    MultiLevelDataValueContainer& operator=(MultiLevelDataValueContainer&& rOther) noexcept;

    virtual ~MultiLevelDataValueContainer() = default;

    /**
     * @brief Adds a variable to the multi-level data container.
     *
     * This function adds a variable along with its associated index accessor to the multi-level data container.
     *
     * @tparam TDataType The data type of the variable.
     * @param rThisVariable The variable to be added.
     * @param pIndexAccesor A unique pointer to the index accessor associated with the variable.
     */
    template<class TDataType>
    void AddVariable(
        const Variable<TDataType>& rThisVariable,
        MultiLevelDataValueAccesor::UniquePointer pIndexAccesor)
    {
        if (Has(rThisVariable)) {
            KRATOS_ERROR << "The variable is already in the database. Remove it first if you want to use a different accesor" << std::endl;
        }
        auto p_source_variable = &rThisVariable.GetSourceVariable();
        mData.emplace_back(p_source_variable, std::move(pIndexAccesor));
    }

    /**
     * @brief Checks if the container has the specified variable.
     * 
     * This function searches through the container to determine if it contains
     * the variable specified by @rThisVariable. It uses the variable's source
     * key to perform the search.
     * 
     * @tparam TDataType The type of the data stored in the variable.
     * @param rThisVariable The variable to check for in the container.
     * @return true if the container has the specified variable, false otherwise.
     */
    template<class TDataType>
    bool Has(
        const Variable<TDataType>& rThisVariable) const
    {
        auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        return it_i != mData.end();
    }

    /**
     * @brief Removes a variable from the multi-level data container.
     * 
     * This function removes the specified variable from the multi-level data container.
     * If the variable is not found in the container, no action is taken.
     * 
     * @tparam TDataType The type of the data stored in the variable.
     * @param rThisVariable The variable to be removed from the container.
     */
    template<class TDataType>
    void RemoveVariable(
        const Variable<TDataType>& rThisVariable)
    {
        auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        if (it_i != mData.end()) {
            mData.erase(it_i);
        }
    }

    /**
     * @brief Retrieves the index accessor for a given variable.
     * 
     * This function searches for the specified variable in the data container and returns
     * its associated index accessor if found. If the variable is not found, an error is thrown.
     * 
     * @param rThisVariable The variable for which the index accessor is to be retrieved.
     * @return Const reference to the index accessor for the specified variable.
     * @throws std::runtime_error If the variable is not found in the data container.
     */
    const MultiLevelDataValueAccesor::UniquePointer& GetIndexAccesor(
        const VariableData& rThisVariable) const;

    /**
     * @brief Retrieves the value of the specified data type.
     * 
     * This function template is used to get the value of a specific data type (TDataType)
     * from the multi-level data container. The data type can be any type that is stored
     * within the container.
     * 
     * @tparam TDataType The type of data to be retrieved.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param Indices A variadic list of indices specifying the location in the multi-level data container.
     * @return TDataType The value of the specified data type.
     */
    template<class TDataType, typename... Args> 
    TDataType& GetValue(
        const Variable<TDataType>& rThisVariable,
        Args&&... Indices)
    {
        auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        if (it_i != mData.end()) {
            const std::size_t index = it_i->GetIndex(std::forward<Args>(Indices)...);
            return *(it_i->pValue<TDataType>(index) + rThisVariable.GetComponentIndex());
        } else {
            KRATOS_ERROR << "The variable was not found in database, please do AddVariable first." << std::endl;
        }
    }

    /**
     * @brief Retrieves the value of the specified data type. (const version)
     * 
     * This function template is used to get the value of a specific data type (TDataType)
     * from the multi-level data container. The data type can be any type that is stored
     * within the container.
     * 
     * @tparam TDataType The type of data to be retrieved.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param Indices A variadic list of indices specifying the location in the multi-level data container.
     * @return TDataType The value of the specified data type.
     */
    template<class TDataType, typename... Args> 
    const TDataType& GetValue(
        const Variable<TDataType>& rThisVariable,
        Args&&... Indices) const
    {
        auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        if (it_i != mData.end()) {
            const std::size_t index = it_i->GetIndex(std::forward<Args>(Indices)...);
            return *(it_i->pValue<TDataType>(index) + rThisVariable.GetComponentIndex());
        }

        return rThisVariable.Zero();
    }

    /**
     * @brief Sets the value of a variable in the multi-level data container.
     * 
     * This function sets the value of a specified variable in the multi-level data container
     * at the given indices. If the variable is not found in the container, an error is thrown.
     * 
     * @tparam TDataType The type of the data to be set.
     * @param rThisVariable The variable whose value is to be set.
     * @param rValue The value to be set for the specified variable.
     * @param Indices A variadic list of indices specifying the location in the multi-level data container.
     * 
     * @throws std::runtime_error If the variable is not found in the data container.
     */
    template<class TDataType, typename... Args> 
    void SetValue(
        const Variable<TDataType>& rThisVariable,
        const TDataType& rValue,
        Args&&... Indices)
    {
        auto it_data = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        if (it_data != mData.end()) {
            const std::size_t index = it_data->GetIndex(std::forward<Args>(Indices)...);
            KRATOS_DEBUG_ERROR_IF(index >= it_data->size()) << "Index out of range" << std::endl;
            *(it_data->pValue<TDataType>(index) + rThisVariable.GetComponentIndex()) = rValue;
        } else {
            KRATOS_ERROR << "The variable was not found in database, please do AddVariable first." << std::endl;
        }
    }


protected:

    class DataBlock
    {
        public:

        using BlockType=double;

        DataBlock() = default;

        DataBlock(
            const VariableData* pVariable,
            MultiLevelDataValueAccesor::UniquePointer pIndexAccesor)
            : mpVariable(pVariable),
              mpIndexAccesor(std::move(pIndexAccesor)),
              mSize(mpIndexAccesor->size())
        {
            mpValues = Allocate(pVariable, mSize);
            const std::size_t allocation_size = mpVariable->Size() / sizeof(BlockType);
            for(std::size_t i = 0; i < mSize; i++) {
                mpVariable->AssignZero(mpValues + i * allocation_size);
            }
        }

        DataBlock(const DataBlock& rOther)
            : mpVariable(rOther.mpVariable),
              mpIndexAccesor(rOther.mpIndexAccesor ? 
                rOther.mpIndexAccesor->Clone() :
                nullptr),
              mSize(rOther.mSize)
        {
            const std::size_t allocation_size = mpVariable->Size() / sizeof(BlockType);
            mpValues = new BlockType[mSize * allocation_size];
            for(std::size_t i = 0; i < mSize; i++) {
                mpVariable->Copy(rOther.mpValues + i * allocation_size, mpValues + i * allocation_size);
            }
        }

        DataBlock(DataBlock&& rOther) noexcept
            : mpVariable(std::exchange(rOther.mpVariable, nullptr)),
              mpIndexAccesor(std::exchange(rOther.mpIndexAccesor, nullptr)),
              mSize(std::exchange(rOther.mSize, 0)),
              mpValues(std::exchange(rOther.mpValues, nullptr))
        {
        }

        DataBlock& operator=(const DataBlock& rOther)
        {
            if (this != &rOther) {
                mpVariable = rOther.mpVariable;
                mpIndexAccesor = rOther.mpIndexAccesor ? 
                    std::make_unique<MultiLevelDataValueAccesor>(*(rOther.mpIndexAccesor->Clone())) :
                    nullptr;
                mSize = rOther.mSize;

                const std::size_t allocation_size = mpVariable->Size() / sizeof(BlockType);
                delete[] mpValues;
                mpValues = new BlockType[mSize * allocation_size];
                for (std::size_t i = 0; i < mSize; i++) {
                    mpVariable->Copy(mpValues + i * allocation_size, rOther.mpValues + i * allocation_size);
                }
            }
            return *this;
        }

        DataBlock& operator=(DataBlock&& rOther) noexcept
        {
            if (this != &rOther) {
                mpVariable = std::exchange(rOther.mpVariable, nullptr);
                mpIndexAccesor = std::exchange(rOther.mpIndexAccesor, nullptr);
                mSize = std::exchange(rOther.mSize, 0);
                delete[] mpValues;
                mpValues = std::exchange(rOther.mpValues, nullptr);
            }
            return *this;
        }

        std::size_t size() const
        {
            return mSize;
        }

        const VariableData& GetVariable() const
        {
            return *mpVariable;
        }

        template<typename... Args>
        std::size_t GetIndex(
            Args&&... Indices) const
        {
            return mpIndexAccesor->GetIndex(std::forward<Args>(Indices)...);
        }

        template<class TDataType>
        TDataType* pValue(std::size_t Index)
        {
            std::size_t offset = Index * mpVariable->Size()/sizeof(BlockType);
            TDataType* p = reinterpret_cast<TDataType*>(mpValues + offset);
            return p;
        }

        template<class TDataType>
        const TDataType* pValue(std::size_t Index) const
        {
            std::size_t offset = Index * mpVariable->Size()/sizeof(BlockType);
            const TDataType* p = reinterpret_cast<const TDataType*>(mpValues + offset);
            return p;
        }

        const BlockType* pPointer(std::size_t Index) const
        {
            std::size_t offset = Index * mpVariable->Size()/sizeof(BlockType);
            return mpValues + offset;
        }

        virtual ~DataBlock()
        {
            if (mpValues == nullptr) {
                return;
            }
            const std::size_t allocation_size = mpVariable->Size()/sizeof(BlockType);
            for(std::size_t i = 0; i < mSize; i++) {
                mpVariable->Destruct(mpValues + i * allocation_size);
            }

            delete[] mpValues;
        }

        const MultiLevelDataValueAccesor::UniquePointer& GetIndexAccesor() const
        {
            return mpIndexAccesor;
        }

    private:

        BlockType* Allocate(const VariableData* pVariable , std::size_t Size)
        {
            std::size_t allocation_size = pVariable->Size()/sizeof(BlockType);
            return new BlockType[Size*allocation_size];
        }

        const VariableData* mpVariable;
        MultiLevelDataValueAccesor::UniquePointer mpIndexAccesor;
        std::size_t mSize;
        BlockType* mpValues;

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const
        {
            rSerializer.save("VariableName", mpVariable->Name());
            rSerializer.save("IndexAccesor", mpIndexAccesor);
            rSerializer.save("Size", mSize);
            for (std::size_t i = 0; i < mSize; i++) {
                std::size_t offset = i * mpVariable->Size()/sizeof(BlockType);
                mpVariable->Save(rSerializer, mpValues + offset);
            }
        }

        virtual void load(Serializer& rSerializer)
        {
            std::string name;
            rSerializer.load("VariableName", name);
            mpVariable = KratosComponents<VariableData>::pGet(name);
            rSerializer.load("IndexAccesor", mpIndexAccesor);
            rSerializer.load("Size", mSize);
            mpValues = Allocate(mpVariable, mSize);
            for (std::size_t i = 0; i < mSize; i++) {
                std::size_t offset = i * mpVariable->Size()/sizeof(BlockType);
                mpVariable->Load(rSerializer, mpValues + offset);
            }
        }
    };

    std::vector<DataBlock> mData;

    /**
    * @brief Functor class used to check if a `ValueType` has a specific index key.
    * @details The `IndexCheck` class is designed to be used with algorithms like `std::find_if` 
    * to search for a `ValueType` with a specific source key.
    */
    class IndexCheck
    {
        std::size_t mI; /// The source key index to be checked against.

    public:

        /**
        * @brief Constructor that initializes the functor with a specific source key index.
        * @param I The source key index to be checked against.
        */
        explicit IndexCheck(std::size_t I) : mI(I) {}

        /**
        * @brief Overloaded function call operator to compare the `ValueType`'s source key
        *  with the stored index.
        * @param I The `ValueType` whose source key is to be compared.
        * @return True if the `ValueType`'s source key matches the stored index, otherwise false.
        */
        bool operator()(const DataBlock& I)
        {
            return I.GetVariable().SourceKey() == mI;
        }
    };

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("Data", mData);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("Data", mData);
    }

};


inline std::istream& operator >> (
    std::istream& rIStream,
    MultiLevelDataValueContainer& rThis)
{
    return rIStream;
}

inline std::ostream& operator << (
    std::ostream& rOStream,
    const MultiLevelDataValueContainer& rThis)
{
    rOStream << "[MultiLevelDataValueContainer]" << std::endl;
    return rOStream;
}


} // namespace kratos 