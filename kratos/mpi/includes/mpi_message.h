//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#pragma once

#include <string>
#include <vector>
#include <type_traits>
#include "mpi.h"

#include "containers/array_1d.h"
#include "containers/flags.h"
#include "utilities/data_type_traits.h"

namespace Kratos
{

namespace Internals {

template<class TDataType> struct MPIDataType;

template<> struct MPIDataType<int>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_INT;
    }
};

template<> struct MPIDataType<unsigned int>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_UNSIGNED;
    }
};

template<> struct MPIDataType<long unsigned int>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_UNSIGNED_LONG;
    }
};

template<> struct MPIDataType<double>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_DOUBLE;
    }
};

template<> struct MPIDataType<bool>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_C_BOOL;
    }
};

template<> struct MPIDataType<char>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_CHAR;
    }
};

template<> struct MPIDataType<int64_t>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_INT64_T;
    }
};

}

template<class TDataType> class MPIMessage
{
public:
    ///@name Type definitions
    ///@{

    using MessageDataType = TDataType;

    using MPIMessageDataTypeTraits = DataTypeTraits<TDataType>;

    using SubDataType = typename MPIMessageDataTypeTraits::ValueType;

    using PrimitiveDataType = typename MPIMessageDataTypeTraits::PrimitiveType;

    static constexpr bool HasContiguousPrimitiveData = MPIMessageDataTypeTraits::HasContiguousPrimitiveData;

    static constexpr bool HasDynamicMemoryAllocation = MPIMessageDataTypeTraits::HasDynamicMemoryAllocation;

    ///@}
    ///@name Public operations
    ///@{

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValue)
    {
        return mData.GetData(rValue);
    }

    inline const void* Buffer(const MessageDataType& rValue)
    {
        return  mData.GetData(rValue);
    }

    inline int Size(const MessageDataType& rValue)
    {
        return MPIMessageDataTypeTraits::Size(rValue);
    }

    inline int SubDataTypeSize(const MessageDataType& rValue)
    {
        if constexpr(!std::is_same_v<SubDataType, TDataType>) {
            return rValue.size() == 0 ? 0 : MPIMessage<SubDataType>().Size(rValue.front());
        } else {
            return 0;
        }
    }

    inline std::vector<unsigned int> Shape(const MessageDataType& rValue)
    {
        return MPIMessageDataTypeTraits::Shape(rValue);
    }

    inline bool Resize(MessageDataType& rValue, const std::vector<unsigned int>& rShape)
    {
        return MPIMessageDataTypeTraits::Reshape(rValue, rShape);
    }

    inline void Update(MessageDataType& rValue)
    {
        mData.UpdateValues(rValue);
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    DataBuffer<TDataType> mData;

    ///@}
};

} // namespace Kratos
