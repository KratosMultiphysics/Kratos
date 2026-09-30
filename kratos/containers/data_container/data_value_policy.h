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
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

#pragma once

// System includes
#include <type_traits>

// External includes

// Project includes
#include "containers/variable.h"
#include "containers/data_container/data_accessor.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class DataValuePolicyBase
 * @ingroup KratosCore
 * @brief Type-erased interface defining how the values of a DataChunk are stored and managed.
 * @details A value policy encapsulates all type-dependent operations on the values kept in a DataChunk (allocation, copy, assignment, zero-initialization, destruction and printing) behind a type-erased interface working on void pointers, together with the comparison semantics used by DataContainer to match chunks: IsSameType (same dynamic policy type, ignoring configuration details such as the zero value or the number of layers), IsSame (same dynamic type and configuration), IsCompatible (whether the policy can be used with a given Variable, i.e. the variable value type matches the policy value type) and IsSparse (whether the policy stores data sparsely; sparse chunks start empty and grow via DataContainer::UpdateSparseStorage / DataContainer::AddToSparseStorage).
 *
 * Policies are usually not used directly: they are passed to DataContainer::Add, which clones them into the created DataChunk. The class provides no thread-safety guarantees; policy objects are immutable after construction and safe to share for reading.
 * @see DataValuePolicy, LayeredDataValuePolicy, SparseDataValuePolicy, DataChunk, DataContainer
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) DataValuePolicyBase
{
public:
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param SizeOfData Size of one value in bytes. Must be greater than zero.
     */
    DataValuePolicyBase(std::size_t SizeOfData);

    /// Destructor.
    virtual ~DataValuePolicyBase() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Allocate a new value copy-constructed from the given source.
     * @param pSource Pointer to the source value.
     * @return Pointer to the newly allocated value (owned by the caller, release with Delete).
     */
    virtual void* CloneData(const void* pSource) const = 0;

    /**
     * @brief Copy-construct a value at an already allocated destination (placement new).
     * @param pSource Pointer to the source value.
     * @param pDestination Pointer to raw memory of at least SizeOfData() bytes.
     * @return Pointer to the constructed value.
     */
    virtual void* Copy(const void* pSource, void* pDestination) const = 0;

    /**
     * @brief Assign the source value to an already constructed destination value.
     * @param pSource Pointer to the source value.
     * @param pDestination Pointer to the destination value.
     */
    virtual void Assign(const void* pSource, void* pDestination) const = 0;

    /**
     * @brief Construct the policy zero value at the given destination (placement new).
     * @param pDestination Pointer to raw memory of at least SizeOfData() bytes.
     */
    virtual void AssignZero(void* pDestination) const = 0;

    /**
     * @brief Destroy and deallocate a value previously created with CloneData or Allocate.
     * @param pSource Pointer to the value to delete.
     */
    virtual void Delete(void* pSource) const = 0;

    /**
     * @brief Call the destructor of the value without releasing its memory.
     * @param pSource Pointer to the value to destruct.
     */
    virtual void Destruct(void* pSource) const = 0;

    /**
     * @brief Print the value to the given stream.
     * @param pSource Pointer to the value to print.
     * @param rOStream The output stream.
     */
    virtual void Print(const void* pSource, std::ostream& rOStream) const = 0;

    /**
     * @brief Print the value with a leading "Data: " label to the given stream.
     * @param pSource Pointer to the value to print.
     * @param rOStream The output stream.
     */
    virtual void PrintData(const void* pSource, std::ostream& rOStream) const = 0;

    /**
     * @brief Allocate a new default-constructed value.
     * @param pData Output: receives the pointer to the newly allocated value (owned by the caller, release with Delete).
     */
    virtual void Allocate(void** pData) const = 0;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Get the size of one value in bytes.
     * @return The size of the data to be allocated per value.
     */
    virtual std::size_t SizeOfData() const;

    /**
     * @brief Check if the given policy has the same dynamic type as this one.
     * @details Configuration details (e.g. the zero value) are ignored; see IsSame for a full comparison.
     * @param rOther The policy to compare against.
     * @return true if both policies are of the same dynamic type.
     */
    virtual bool IsSameType(const DataValuePolicyBase& rOther) const = 0;

    /**
     * @brief Check if the given policy is the same as this one (type and configuration).
     * @param rOther The policy to compare against.
     * @return true if both policies are of the same dynamic type with identical configuration.
     */
    virtual bool IsSame(const DataValuePolicyBase& rOther) const = 0;

    /**
     * @brief Check if this policy can be used with the given variable.
     * @details The base overload always returns false; DataValuePolicy<TValueType> shadows it, returning true when the variable value type matches the policy value type. Note this is a non-virtual template shadowing pattern: the check resolves on the static type of the policy, which is the concrete type at the DataContainer::Add call site.
     * @tparam TVariableValueType The value type of the variable.
     * @param rVariable The variable to check against.
     * @return false (in the base class).
     */
    template<typename TVariableValueType>
    bool IsCompatible(const Variable<TVariableValueType>& rVariable) const
    {
        return false;
    }

    /**
     * @brief Check if the policy is expected to work with sparse data.
     * @return false unless overridden (see SparseDataValuePolicy).
     */
    virtual bool IsSparse() const;

    /**
     * @brief Get the accessor of the sparse index variable driving the sparse storage.
     * @details Only meaningful for sparse policies; the base implementation raises an error.
     * @return The accessor of the sparse index variable.
     */
    virtual DataAccessor<int> GetSparseIndexAccessor() const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::size_t mSizeOfData; /// Size of the data to be allocated in bytes

    ///@}

}; // Class DataValuePolicyBase

/**
 * @class DataValuePolicy
 * @ingroup KratosCore
 * @brief Value policy implementing dense, contiguous storage of values of type TValueType.
 * @details Implements the type-erased operations of DataValuePolicyBase for a concrete value type, and holds the zero value used to initialize newly allocated entries of a DataChunk (note this zero belongs to the policy, not to the Variable; it can be customized per policy instance, e.g. DataValuePolicy<int>(-1) for sparse index variables). Two DataValuePolicy instances IsSame only when their zero values match.
 * @see DataValuePolicyBase, LayeredDataValuePolicy, SparseDataValuePolicy, DataChunk
 * @tparam TValueType The type of the data value to be stored.
 * @author Pooyan Dadvand
 */
template<typename TValueType>
class DataValuePolicy : public DataValuePolicyBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Type alias for the data type
    using ValueType = TValueType;

    /// Span type for the data. Used to access step data in a DataChunk
    using StepSpanType = Kratos::span<TValueType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @details The default value of rZero relies on TValueType()'s default constructor producing a meaningful zero, which holds for scalar types, std::string and std::vector, but not for lightweight ublas-derived types such as array_1d, Vector or Matrix, whose default constructors deliberately leave their storage uninitialized for performance; callers using such value types must pass an explicit rZero.
     * @param rZero The value used to zero-initialize entries governed by this policy.
     */
    DataValuePolicy(TValueType const& rZero = TValueType())
        : DataValuePolicyBase(sizeof(TValueType)), mZero(rZero)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone the policy.
     * @return A newly allocated copy of this policy.
     */
    virtual std::unique_ptr<DataValuePolicy<TValueType>> Clone() const
    {
        return std::make_unique<DataValuePolicy<TValueType>>(*this);
    }

    void* CloneData(const void* pSource) const override
    {
        return new TValueType(*static_cast<const TValueType*>(pSource));
    }

    void* Copy(const void* pSource, void* pDestination) const override
    {
        return new(pDestination) TValueType(*static_cast<const TValueType*>(pSource));
    }

    void Assign(const void* pSource, void* pDestination) const override
    {
        (*static_cast<TValueType*>(pDestination)) = (*static_cast<const TValueType*>(pSource));
    }

    void AssignZero(void* pDestination) const override
    {
        new(pDestination) TValueType(mZero);
    }

    void Delete(void* pSource) const override
    {
        delete static_cast<TValueType*>(pSource);
    }

    void Destruct(void* pSource) const override
    {
        static_cast<TValueType*>(pSource)->~TValueType();
    }

    void Print(const void* pSource, std::ostream& rOStream) const override
    {
        rOStream << *static_cast<const TValueType*>(pSource);
    }

    void PrintData(const void* pSource, std::ostream& rOStream) const override
    {
        rOStream << "Data: " << *static_cast<const TValueType*>(pSource);
    }

    void Allocate(void** pData) const override
    {
        *pData = new TValueType();
    }

    ///@}
    ///@name Inquiry
    ///@{

    bool IsSameType(const DataValuePolicyBase& rOther) const override
    {
        return dynamic_cast<const DataValuePolicy<TValueType>*>(&rOther) != nullptr;
    }

    bool IsSame(const DataValuePolicyBase& rOther) const override
    {
        // Check if the other policy is of the same type
        const DataValuePolicy<TValueType>* p_other_policy = dynamic_cast<const DataValuePolicy<TValueType>*>(&rOther);
        if (!p_other_policy) {
            return false; // different type
        }
        // check if all data members are the same
        return (mZero == p_other_policy->mZero &&
                SizeOfData() == p_other_policy->SizeOfData());
    }

    /**
     * @brief Check if this policy can be used with the given variable.
     * @details Shadows DataValuePolicyBase::IsCompatible; returns true only when the variable value type matches the policy value type.
     * @tparam TVariableValueType The value type of the variable.
     * @param rVariable The variable to check against.
     * @return true if the variable value type matches TValueType.
     */
    template<typename TVariableValueType>
    bool IsCompatible(const Variable<TVariableValueType>& rVariable) const
    {
        if constexpr (std::is_same_v<TValueType, TVariableValueType>) {
            return true;
        }

        return false;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the zero value of this policy.
     * @return The value used to zero-initialize entries governed by this policy.
     */
    const TValueType& Zero() const
    {
        return mZero;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    TValueType mZero; /// Zero value of the data type

    ///@}

}; // Class DataValuePolicy

/**
 * @class LayeredDataValuePolicy
 * @ingroup KratosCore
 * @brief Value policy storing a fixed number of layers (a std::vector<TValueType>) per entity.
 * @details Each entry governed by this policy is a std::vector<TValueType> of NumberOfLayers() elements, all zero-initialized with the given per-layer zero value. Two layered policies IsSame only when their number of layers matches (IsSameType ignores it).
 * @see DataValuePolicy, DataChunk
 * @tparam TValueType The type of the per-layer value.
 * @author Pooyan Dadvand
 */
template<typename TValueType>
class LayeredDataValuePolicy : public DataValuePolicy<std::vector<TValueType>>
{
public:
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param NumberOfLayers The number of layers to store per entity. Must be greater than zero.
     * @param rZero The per-layer zero value.
     */
    LayeredDataValuePolicy(std::size_t NumberOfLayers, TValueType const& rZero = TValueType())
        : DataValuePolicy<std::vector<TValueType>>(std::vector<TValueType>(NumberOfLayers, rZero)),
          mNumberOfLayers(NumberOfLayers)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(NumberOfLayers > 0) << "Number of layers must be greater than zero." << std::endl;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone the policy.
     * @return A newly allocated copy of this policy.
     */
    std::unique_ptr<DataValuePolicy<std::vector<TValueType>>> Clone() const override
    {
        return std::make_unique<LayeredDataValuePolicy<TValueType>>(*this);
    }

    ///@}
    ///@name Inquiry
    ///@{

    bool IsSameType(const DataValuePolicyBase& rOther) const override
    {
        return dynamic_cast<const LayeredDataValuePolicy<TValueType>*>(&rOther) != nullptr;
    }

    bool IsSame(const DataValuePolicyBase& rOther) const override
    {
        // Check if the other policy is of the same type
        const LayeredDataValuePolicy<TValueType>* p_other_policy = dynamic_cast<const LayeredDataValuePolicy<TValueType>*>(&rOther);
        if (!p_other_policy) {
            return false; // different type
        }
        // check if all data members are the same
        return (mNumberOfLayers == p_other_policy->mNumberOfLayers);
    }

    /**
     * @brief Check if this policy can be used with the given variable.
     * @details A layered policy is compatible with variables of the per-layer value type (not with variables of the stored std::vector type).
     * @tparam TVariableValueType The value type of the variable.
     * @param rVariable The variable to check against.
     * @return true if the variable value type matches TValueType.
     */
    template<typename TVariableValueType>
    bool IsCompatible(const Variable<TVariableValueType>& rVariable) const
    {
        if constexpr (std::is_same_v<TValueType, TVariableValueType>) {
            return true;
        }

        return false;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the number of layers stored per entity.
     * @return The number of layers.
     */
    std::size_t NumberOfLayers() const
    {
        return mNumberOfLayers;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::size_t mNumberOfLayers; /// Number of layers to be stored

    ///@}

}; // Class LayeredDataValuePolicy

///@}

///@} addtogroup block

} // namespace Kratos
