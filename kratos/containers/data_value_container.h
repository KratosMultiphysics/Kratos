//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <cstddef>

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/kratos_components.h"
#include "includes/exception.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class DataValueContainer
 * @ingroup KratosCore
 * @brief Container for storing data values associated with variables.
 * @details This class provides a container for storing data values associated with variables.
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) DataValueContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Define local flag
    KRATOS_DEFINE_LOCAL_FLAG(OVERWRITE_OLD_VALUES);

    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(DataValueContainer);

    /// Type of the container used for variables
    using ValueType = std::pair<const VariableData*, void*>;

    /// Type of the container used for variables
    using ContainerType = std::vector<ValueType>;

    /// Type of the container used for variables
    using iterator = ContainerType::iterator;

    /// Type of the container used for variables
    using const_iterator = ContainerType::const_iterator;

    /// Type of the container used for variables
    using SizeType = ContainerType::size_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataValueContainer() noexcept = default;

    /// Move constructor.
    DataValueContainer(DataValueContainer&&) noexcept = default;

    /// Copy constructor.
    DataValueContainer(DataValueContainer const& rOther)
    {
        for(const_iterator i = rOther.mData.begin() ; i != rOther.mData.end() ; ++i)
            mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
    }

    /// Destructor.
    virtual ~DataValueContainer()
    {
        for(iterator i = mData.begin() ; i != mData.end() ; ++i)
            i->first->Delete(i->second);
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Accessor operator for retrieving a data value by a VariableData.
     * @details This operator allows you to retrieve a data value by providing a VariableData object.
     * @tparam TDataType The data type of the value to retrieve.
     * @param rThisVariable The VariableData object specifying the variable.
     * @return A reference to the data value.
     */
    template<class TDataType>
    const TDataType& operator()(const VariableData& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    /**
     * @brief Accessor operator for retrieving a data value by a Variable.
     * @details This operator allows you to retrieve a data value by providing a Variable object.
     * @tparam TDataType The data type of the value to retrieve.
     * @param rThisVariable The Variable object specifying the variable.
     * @return A reference to the data value.
     */
    template<class TDataType>
    TDataType& operator()(const Variable<TDataType>& rThisVariable)
    {
        return GetValue<TDataType>(rThisVariable);
    }

    /**
     * @brief Accessor operator for retrieving a data value by a Variable (const version).
     * @details This operator allows you to retrieve a data value by providing a Variable object.
     * @tparam TDataType The data type of the value to retrieve.
     * @param rThisVariable The Variable object specifying the variable.
     * @return A reference to the data value.
     */
    template<class TDataType>
    const TDataType& operator()(const Variable<TDataType>& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    /**
     * @brief Index operator for retrieving a data value by a VariableData.
     * @details This operator allows you to retrieve a data value by providing a VariableData object.
     * @tparam TDataType The data type of the value to retrieve.
     * @param rThisVariable The VariableData object specifying the variable.
     * @return A reference to the data value.
     */
    template<class TDataType>
    TDataType& operator[](const VariableData& rThisVariable)
    {
        return GetValue<TDataType>(rThisVariable);
    }

    /**
     * @brief Index operator for retrieving a data value by a VariableData (const version).
     * @details This operator allows you to retrieve a data value by providing a VariableData object.
     * @tparam TDataType The data type of the value to retrieve.
     * @param rThisVariable The VariableData object specifying the variable.
     * @return A reference to the data value.
     */
    template<class TDataType>
    const TDataType& operator[](const VariableData& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    /**
     * @brief Index operator for retrieving a data value by a Variable.
     * @details This operator allows you to retrieve a data value by providing a Variable object.
     * @tparam TDataType The data type of the value to retrieve.
     * @param rThisVariable The Variable object specifying the variable.
     * @return A reference to the data value.
     */
    template<class TDataType>
    TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue<TDataType>(rThisVariable);
    }

    /**
     * @brief Index operator for retrieving a data value by a Variable (const version).
     * @details This operator allows you to retrieve a data value by providing a Variable object.
     * @tparam TDataType The data type of the value to retrieve.
     * @param rThisVariable The Variable object specifying the variable.
     * @return A reference to the data value.
     */
    template<class TDataType>
    const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue<TDataType>(rThisVariable);
    }

    /**
     * @brief Iterator pointing to the beginning of the container.
     * @return An iterator pointing to the beginning of the container.
     */
    iterator begin()
    {
        return mData.begin();
    }

    /**
     * @brief Const iterator pointing to the beginning of the container.
     * @return A const iterator pointing to the beginning of the container.
     */
    const_iterator begin() const
    {
        return mData.begin();
    }

    /**
     * @brief Iterator pointing to the end of the container.
     * @return An iterator pointing to the end of the container.
     */
    iterator end()
    {
        return mData.end();
    }

    /**
     * @brief Const iterator pointing to the end of the container.
     * @return A const iterator pointing to the end of the container.
     */
    const_iterator end() const
    {
        return mData.end();
    }

    /**
     * @brief Move assignment operator
     * @param[in] rOther The temporary or explicitly moved DataValueContainer object to transfer resources from.
     */
    DataValueContainer& operator=(DataValueContainer&& rOther) noexcept = default;

    /**
     * @brief Assignment operator for copying data from another DataValueContainer.
     * @details This operator allows you to assign the contents of another DataValueContainer to this container.
     * @param rOther The DataValueContainer to copy data from.
     * @return A reference to the modified DataValueContainer.
     */
    DataValueContainer& operator=(const DataValueContainer& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Gets the value associated with a given variable.
     * @tparam TDataType The data type of the variable.
     * @param rThisVariable The variable for which the value is to be retrieved.
     * @return Reference to the value of the provided variable.
     */
    template<class TDataType>
    TDataType& GetValue(const Variable<TDataType>& rThisVariable)
    {
        typename ContainerType::iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()))) != mData.end()) {
            return *(static_cast<TDataType*>(i->second) + rThisVariable.GetComponentIndex());
        }

    #ifdef KRATOS_SMP_OPENMP
        KRATOS_DEBUG_ERROR_IF(static_cast<bool>(omp_in_parallel())) << "Attempting to do a GetValue for: " << rThisVariable << " unfortunately the variable is not in the database and the operations is not threadsafe (this function is being called from within a parallel region)" << std::endl;
    #endif

        auto p_source_variable = &rThisVariable.GetSourceVariable();
        mData.push_back(ValueType(p_source_variable,p_source_variable->Clone(p_source_variable->pZero())));

        return *(static_cast<TDataType*>(mData.back().second) + rThisVariable.GetComponentIndex());
    }

    /**
     * @brief Gets the value associated with a given variable (const version).
     * @tparam TDataType The data type of the variable.
     * @param rThisVariable The variable for which the value is to be retrieved.
     * @return Const reference to the value of the provided variable.
     * @todo Make the variable of the constant version consistent with the one of the "classical" one
     */
    template<class TDataType>
    const TDataType& GetValue(const Variable<TDataType>& rThisVariable) const
    {
        typename ContainerType::const_iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()))) != mData.end()) {
            return *(static_cast<const TDataType*>(i->second) + rThisVariable.GetComponentIndex());
        }

        return rThisVariable.Zero();
    }

    /**
     * @brief Gets the size of the data container.
     * @return Size of the data container.
     */
    SizeType Size()
    {
        return mData.size();
    }

    /**
     * @brief Sets the value for a given variable.
     * @tparam TDataType The data type of the variable.
     * @param rThisVariable The variable for which the value is to be set.
     * @param rValue The value to be set for the variable.
     */
    template<class TDataType>
    void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
    {
        typename ContainerType::iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey())))  != mData.end()) {
            *(static_cast<TDataType*>(i->second) + rThisVariable.GetComponentIndex()) = rValue;
        } else {
            auto p_source_variable = &rThisVariable.GetSourceVariable();
            mData.push_back(ValueType(p_source_variable,p_source_variable->Clone(p_source_variable->pZero())));
            *(static_cast<TDataType*>(mData.back().second) + rThisVariable.GetComponentIndex()) = rValue;
        }
    }

    /**
     * @brief Erases the value associated with a given variable.
     * @tparam TDataType The data type of the variable.
     * @param rThisVariable The variable whose associated value is to be erased.
     */
    template<class TDataType>
    void Erase(const Variable<TDataType>& rThisVariable)
    {
        typename ContainerType::iterator i;

        if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey()))) != mData.end()) {
            i->first->Delete(i->second);
            mData.erase(i);
        }
    }

    /**
     * @brief Clears the entire data container.
     */
    void Clear()
    {
        for(ContainerType::iterator i = mData.begin() ; i != mData.end() ; i++)
            i->first->Delete(i->second);

        mData.clear();
    }

    /**
     * @brief Merges this data container with another data container.
     * @param rOther The other data container to be merged.
     * @param Options Flags for the merging options.
     */
    void Merge(const DataValueContainer& rOther, const Flags Options);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if the data container has a value associated with a given variable.
     * @tparam TDataType The data type of the variable.
     * @param rThisVariable The variable for which the check is to be made.
     * @return True if the variable has an associated value in the container, otherwise false.
     */
    template<class TDataType>
    bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return (std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.SourceKey())) != mData.end());
    }

    /**
     * @brief Checks if the data container is empty.
     * @return True if the container is empty, otherwise false.
     */
    bool IsEmpty() const
    {
        return mData.empty();
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Retrieves a string representation of the data value container.
     * @return A string that describes the data value container.
     */
    virtual std::string Info() const
    {
        return std::string("data value container");
    }

    /**
     * @brief Outputs a brief description of the data value container to a given stream.
     * @param rOStream The output stream to which the information should be printed.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "data value container";
    }

    /**
     * @brief Outputs the detailed data contents of the data value container to a given stream.
     * @details For each data entry, it prints the variable type and its associated value.
     * @param rOStream The output stream to which the data should be printed.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
        for(const_iterator i = mData.begin() ; i != mData.end() ; ++i) {
            rOStream <<"    ";
            i->first->Print(i->second, rOStream);
            rOStream << std::endl;
        }
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@{

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
        bool operator()(const ValueType& I)
        {
            return I.first->SourceKey() == mI;
        }
    };

    ///@}
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ContainerType mData; /// The data container considered

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    /**
     * @class Serializer
     * @brief A fried class responsible for handling the serialization process.
     */
    friend class Serializer;

    /**
     * @brief Extract the object's state and uses the Serializer to store it.
     * @param rSerializer Serializer instance to be used for saving.
     */
    virtual void save(Serializer& rSerializer) const;

    /**
     * @brief Set the object's state based on data provided by the Serializer.
     * @param rSerializer Serializer instance to be used for loading.
     */
    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class DataValueContainer

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DataValueContainer& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DataValueContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.