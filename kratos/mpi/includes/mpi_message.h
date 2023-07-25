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
#include "mpi.h"

#include "containers/array_1d.h"
#include "containers/flags.h"

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

template<class TDataType> class MPIMessage;

template<class TDataType> class MPIMessage
{
public:
    using MessageDataType = TDataType;

    using PrimitiveDataType = TDataType;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValue)
    {
        return &rValue;
    }

    static inline const void* Buffer(const MessageDataType& rValue)
    {
        return &rValue;
    }

    static inline int Size(const MessageDataType&)
    {
        return 1;
    }
};

template<class TDataType, std::size_t Dimension> class MPIMessage<array_1d<TDataType, Dimension>>
{
public:
    using MessageDataType = array_1d<TDataType, Dimension>;

    using PrimitiveDataType = TDataType;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValues)
    {
        return rValues.data().data();
    }

    static inline const void* Buffer(const MessageDataType& rValues)
    {
        return rValues.data().data();
    }

    static inline int Size(const MessageDataType& rValues)
    {
        return Dimension;
    }
};

template<> class MPIMessage<std::string>
{
public:
    using MessageDataType = std::string;

    using PrimitiveDataType = char;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValue)
    {
        return rValue.data();
    }

    static inline const void* Buffer(const MessageDataType& rValue)
    {
        return rValue.data();
    }

    static inline int Size(const MessageDataType& rValues)
    {
        return rValues.size();
    }
};

template<class TDataType> class MPIMessage<std::vector<TDataType>>
{
public:
    using MessageDataType = std::vector<TDataType>;

    using PrimitiveDataType = typename MPIMessage<TDataType>::PrimitiveDataType;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValue)
    {
        return rValue.data();
    }

    static inline const void* Buffer(const MessageDataType& rValue)
    {
        return rValue.data();
    }

    static inline int Size(const MessageDataType& rValue)
    {
        return rValue.size();
    }
};

} // namespace Kratos
