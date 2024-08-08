//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#pragma once

// System includes
#include <cmath>

// Project includes
#include "includes/data_communicator.h"
#include "utilities/parallel_utilities.h"

// Application incldues
#include "custom_utilities/data_type_traits.h"

namespace Kratos
{

///@name Classes
///@{

class GenericReductionUtilities
{
private:
    ///@name Private class definitions
    ///@{

    template<class TOperationType>
    class GenericReducer
    {
    public:
        using return_type = TOperationType;

        return_type mValue;

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(const return_type& rValue){
            mValue.Execute(rValue);
        }

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(const GenericReducer<TOperationType>& rOther)
        {
            KRATOS_CRITICAL_SECTION
            mValue.Execute(rOther.mValue);
        }
    };

    ///@}

public:
    ///@name Static operations
    ///@{

    template<class TDataContainerType, class TNormType, template <class T1> class TOperationType, bool IdRequired, int TPower = 1>
    static TOperationType<typename TNormType::template ResultantValueType<typename TDataContainerType::DataType>> GenericReduction(
        const DataCommunicator& rDataCommunicator,
        const TDataContainerType& rDataContainer,
        const TNormType& rNorm)
    {
        KRATOS_TRY

        using entity_data_type = typename TDataContainerType::DataType;

        using norm_resultant_type = typename TNormType::template ResultantValueType<typename TDataContainerType::DataType>;

        using operation_type = TOperationType<norm_resultant_type>;

        using reducer_type = GenericReducer<operation_type>;

        auto local_value = IndexPartition<IndexType>(rDataContainer.Size()).for_each<reducer_type>(entity_data_type{}, [&rDataContainer, &rNorm](const IndexType Index, entity_data_type& rTLS) {
            rDataContainer.GetValue(rTLS, Index);
            auto norm_value = rNorm.Evaluate(rTLS);
            for (IndexType i = 0; i < DataTypeTraits<norm_resultant_type>::Size(norm_value); ++i) {
                DataTypeTraits<norm_resultant_type>::GetComponent(norm_value, i) = std::pow(DataTypeTraits<norm_resultant_type>::GetComponent(norm_value, i), TPower);
            }

            if constexpr(!IdRequired) {
                return operation_type(norm_value);
            } else {
                return operation_type(norm_value, rDataContainer.GetId(Index));
            }

        });

        local_value.Synchronize(rDataCommunicator);

        return local_value;

        KRATOS_CATCH("");
    }

    ///@}
};

} // namespace Kratos