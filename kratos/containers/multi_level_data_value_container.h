#pragma once

#include <vector>

#include "containers/data_value_container.h"
#include "containers/variable.h"

namespace Kratos 
{


class MultiLevelDataValueContainer
{
public:

    /// class to hold the values
    class DataBlock
    {
        public:

        using BlockType=double;

        DataBlock(const VariableData* pVariable, std::size_t Size)
            : mpVariable(pVariable)
            , mSize(Size)
            , mpValues(Allocate(pVariable, Size))
        {
        }

        const VariableData& GetVariable() const
        {
            return *mpVariable;
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

        // destructor
        ~DataBlock()
        {
            for(std::size_t i = 0; i < mSize; i++)
                mpVariable->Delete(mpValues + i);

            delete[] mpValues;
        }

        // private:

        BlockType* Allocate(const VariableData* pVariable , std::size_t Size)
        {
            std::size_t allocation_size = pVariable->Size()/sizeof(BlockType);
            KRATOS_WATCH(allocation_size)
            KRATOS_WATCH(Size)
            return new BlockType[Size*allocation_size];
        }

        const VariableData* mpVariable;
        std::size_t mSize;
        BlockType* mpValues;
    };


    /// Type of the container used for variables
    using ContainerType = std::vector<DataBlock>;

    MultiLevelDataValueContainer()
    {
    }

    template<class TDataType, class TAccessIndex> 
    TDataType& GetValue(
        const Variable<TDataType>& rThisVariable,
        const TAccessIndex& rIndexAccesor,
        const std::vector<std::size_t>& rIndeces)
    {
        auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        if (it_i != mData.end()) {
            return *(it_i->pValue<TDataType>(rIndexAccesor.GetIndex(rIndeces)) + rThisVariable.GetComponentIndex());
        }

#ifdef KRATOS_DEBUG
        if(OpenMPUtils::IsInParallel() != 0)
            KRATOS_ERROR << "attempting to do a GetValue for: " << rThisVariable << " unfortunately the variable is not in the database and the operations is not threadsafe (this function is being called from within a parallel region)" << std::endl;
#endif

        auto p_source_variable = &rThisVariable.GetSourceVariable();
        mData.emplace_back(p_source_variable,rIndexAccesor.Size());

        return *(mData.back().pValue<TDataType>(rIndexAccesor.GetIndex(rIndeces))); // + rThisVariable.GetComponentIndex());
    }

    template<class TDataType, class TAccessIndex> 
    const TDataType& GetValue(
        const Variable<TDataType>& rThisVariable,
        const TAccessIndex& rIndexAccesor,
        const std::vector<std::size_t>& rIndeces) const
    {
        auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        if (it_i != mData.end()) {
            return *(it_i->pValue<TDataType>(rIndexAccesor.GetIndex(rIndeces)) + rThisVariable.GetComponentIndex());
        }

        return rThisVariable.Zero();
    }

    template<class TDataType, class TAccessIndex> 
    void SetValue(
        const Variable<TDataType>& rThisVariable,
        const TAccessIndex& rIndexAccesor,
        const std::vector<std::size_t>& rIndeces,
        const TDataType& rValue)
    {
        auto it_i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()));
        if (it_i != mData.end()) {
            *(it_i->pValue<TDataType>(rIndexAccesor.GetIndex(rIndeces)) + rThisVariable.GetComponentIndex()) = rValue;
        }

        auto p_source_variable = &rThisVariable.GetSourceVariable();
        mData.emplace_back(p_source_variable, rIndexAccesor.Size());
        const std::size_t component_index = rThisVariable.GetComponentIndex();
        const std::size_t accesor_index = rIndexAccesor.GetIndex(rIndeces);
        TDataType* p_value = mData.back().pValue<TDataType>(accesor_index); //missing rThisVariable.GetComponentIndex()
        *p_value = rValue;
    }

protected:
    ContainerType mData;

private:
    /**
    * @brief Functor class used to check if a `ValueType` has a specific index key.
    * @details The `IndexCheck` class is designed to be used with algorithms like `std::find_if` to search for a `ValueType` with a specific source key.
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
        * @brief Overloaded function call operator to compare the `ValueType`'s source key with the stored index.
        * @param I The `ValueType` whose source key is to be compared.
        * @return True if the `ValueType`'s source key matches the stored index, otherwise false.
        */
        bool operator()(const DataBlock& I)
        {
            return I.GetVariable().SourceKey() == mI;
        }
    };

};


} // namespace kratos 