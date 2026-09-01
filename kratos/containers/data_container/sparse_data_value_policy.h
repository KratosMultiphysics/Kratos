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

// External includes

// Project includes
#include "containers/data_container/data_accessor.h"
#include "containers/data_container/data_value_policy.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class SparseDataValuePolicy
 * @ingroup KratosCore
 * @brief Value policy for chunks that store values only for a sparse subset of entities.
 * @details A sparse chunk does not allocate one value per entity; instead it is driven by a dense index variable (a Variable<int> stored in the same DataContainer, referenced through the index accessor passed at construction): entities with a non-negative index own the value at that position of the sparse chunk, entities with a negative index own no value. Chunks created with a sparse policy start empty (zero entities, see DataContainer::Add) and grow via DataContainer::UpdateSparseStorage (rebuild to the number of active indices) or DataContainer::AddToSparseStorage (append entries for new entities, preserving existing values).
 * @see DataValuePolicy, DataAccessor, DataContainer
 * @tparam TValueType The type of the data value to be stored.
 * @author Pooyan Dadvand
 */
template<typename TValueType>
class SparseDataValuePolicy : public DataValuePolicy<TValueType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Type of the accessor of the sparse index variable
    using IndexAccessorType = DataAccessor<int>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param IndexAccessor Accessor of the Variable<int> holding the sparse index for each entity (non-negative for active entries, negative otherwise). It must belong to the same DataContainer in which this policy is used.
     * @param rZero The value used to zero-initialize entries governed by this policy.
     */
    SparseDataValuePolicy(IndexAccessorType IndexAccessor, const TValueType& rZero = TValueType())
        : DataValuePolicy<TValueType>(rZero),
          mIndexAccessor(IndexAccessor)
    {
    }

    /// Copy constructor.
    SparseDataValuePolicy(const SparseDataValuePolicy<TValueType>& rOther)
        : DataValuePolicy<TValueType>(rOther),
          mIndexAccessor(rOther.mIndexAccessor)
    {
    }

    /// Destructor.
    ~SparseDataValuePolicy() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone the policy.
     * @return A newly allocated copy of this policy (preserving the index accessor).
     */
    std::unique_ptr<DataValuePolicy<TValueType>> Clone() const override
    {
        return std::make_unique<SparseDataValuePolicy<TValueType>>(*this);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if the policy is expected to work with sparse data.
     * @return true.
     */
    bool IsSparse() const override
    {
        return true;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the accessor of the sparse index variable driving the sparse storage.
     * @return The accessor of the sparse index variable.
     */
    DataAccessor<int> GetSparseIndexAccessor() const override
    {
        return mIndexAccessor;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    IndexAccessorType mIndexAccessor; /// Accessor of the sparse index variable

    ///@}

}; // Class SparseDataValuePolicy

///@}

///@} addtogroup block

} // namespace Kratos
