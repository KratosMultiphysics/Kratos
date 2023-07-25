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

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValue)
    {
        return &rValue;
    }

    inline const void* Buffer(const MessageDataType& rValue)
    {
        return &rValue;
    }

    inline int Size(const MessageDataType&)
    {
        return 1;
    }

    inline std::vector<unsigned int> Shape(const MessageDataType&)
    {
        return {};
    }

    inline bool Resize(MessageDataType&, const std::vector<unsigned int>&)
    {
        return false;
    }
};

template<class TDataType, std::size_t Dimension> class MPIMessage<array_1d<TDataType, Dimension>>
{
public:
    using MessageDataType = array_1d<TDataType, Dimension>;

    using PrimitiveDataType = TDataType;

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValues)
    {
        return rValues.data().data();
    }

    inline const void* Buffer(const MessageDataType& rValues)
    {
        return rValues.data().data();
    }

    inline int Size(const MessageDataType& rValues)
    {
        return Dimension;
    }

    inline std::vector<unsigned int> Shape(const MessageDataType&)
    {
        return {static_cast<unsigned int>(Dimension)};
    }

    inline bool Resize(MessageDataType&, const std::vector<unsigned int>&)
    {
        return false;
    }
};

template<> class MPIMessage<std::string>
{
public:
    using MessageDataType = std::string;

    using PrimitiveDataType = char;

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValue)
    {
        return rValue.data();
    }

    inline const void* Buffer(const MessageDataType& rValue)
    {
        return rValue.data();
    }

    inline int Size(const MessageDataType& rValues)
    {
        return rValues.size();
    }

    inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size())};
    }

    inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() != 1) << "Invalid shape provided for std::string.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0]);
            return true;
        }

        return false;
    }
};

template<> class MPIMessage<Vector>
{
public:
    using MessageDataType = Vector;

    using PrimitiveDataType = double;

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValue)
    {
        return rValue.data().begin();
    }

    inline const void* Buffer(const MessageDataType& rValue)
    {
        return rValue.data().begin();
    }

    inline int Size(const MessageDataType& rValues)
    {
        return rValues.size();
    }

    inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size())};
    }

    inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() != 1) << "Invalid shape provided for Vector.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0], false);
            return true;
        }

        return false;
    }
};

template<> class MPIMessage<Matrix>
{
public:
    using MessageDataType = Matrix;

    using PrimitiveDataType = double;

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValue)
    {
        return rValue.data().begin();
    }

    inline const void* Buffer(const MessageDataType& rValue)
    {
        return rValue.data().begin();
    }

    inline int Size(const MessageDataType& rValues)
    {
        return rValues.size1() * rValues.size1();
    }

    inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size1()), static_cast<unsigned int>(rValues.size2())};
    }

    inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() != 2) << "Invalid shape provided for Matrix.";

        if (rValues.size1() != rShape[0] || rValues.size2() != rShape[1]) {
            rValues.resize(rShape[0], rShape[0], false);
            return true;
        }

        return false;
    }
};

template<class TDataType> class MPIMessage<std::vector<TDataType>>
{
public:
    using MessageDataType = std::vector<TDataType>;

    using PrimitiveDataType = typename MPIMessage<TDataType>::PrimitiveDataType;

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValue)
    {
        return rValue.data();
    }

    inline const void* Buffer(const MessageDataType& rValue)
    {
        return rValue.data();
    }

    inline int Size(const MessageDataType& rValue)
    {
        return rValue.size();
    }

    inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size())};
    }

    inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() != 1) << "Invalid shape provided for std::vector.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0]);
            return true;
        }

        return false;
    }
};

} // namespace Kratos
