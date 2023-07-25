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

}

template<class TDataType> class MPIMessage
{
public:
    using MessageDataType = TDataType;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<TDataType>::DataType();
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

    static inline std::vector<unsigned int> Shape(const MessageDataType&)
    {
        return {};
    }

    static inline bool Resize(MessageDataType&, const std::vector<unsigned int>&)
    {
        return false;
    }

    static constexpr bool IsContiguous()
    {
        return true;
    }
};

template<class TDataType, std::size_t Dimension> class MPIMessage<array_1d<TDataType, Dimension>>
{
public:
    using MessageDataType = array_1d<TDataType, Dimension>;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<TDataType>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValues)
    {
        return rValues.data().data();
    }

    static inline const void* Buffer(const MessageDataType& rValues)
    {
        return rValues.data().data();
    }

    static inline int Size(const MessageDataType&)
    {
        return Dimension;
    }

    static inline std::vector<unsigned int> Shape(const MessageDataType&)
    {
        return {static_cast<unsigned int>(Dimension)};
    }

    static inline bool Resize(MessageDataType&, const std::vector<unsigned int>&)
    {
        return false;
    }

    static constexpr bool IsContiguous()
    {
        return true;
    }
};

template<> class MPIMessage<Vector>
{
public:
    using MessageDataType = Vector;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<double>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValues)
    {
        return rValues.data().begin();
    }

    static inline const void* Buffer(const MessageDataType& rValues)
    {
        return rValues.data().begin();
    }

    static inline int Size(const MessageDataType& rValues)
    {
        return rValues.size();
    }

    static inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size())};
    }

    static inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rShape.size() == 1) << "Invalid shape for vector.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0], false);
            return true;
        }

        return false;
    }

    static constexpr bool IsContiguous()
    {
        return true;
    }
};

template<> class MPIMessage<Matrix>
{
public:
    using MessageDataType = Matrix;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<double>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValues)
    {
        return rValues.data().begin();
    }

    static inline const void* Buffer(const MessageDataType& rValues)
    {
        return rValues.data().begin();
    }

    static inline int Size(const MessageDataType& rValues)
    {
        return rValues.size1() * rValues.size2();
    }

    static inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size1()), static_cast<unsigned int>(rValues.size2())};
    }

    static inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rShape.size() == 2) << "Invalid shape for matrix.";

        if (rValues.size1() != rShape[0] || rValues.size2() != rShape[1]) {
            rValues.resize(rShape[0], rShape[1], false);
            return true;
        }

        return false;
    }

    static constexpr bool IsContiguous()
    {
        return true;
    }
};

template<class TDataType> class MPIMessage<std::vector<TDataType>>
{
public:
    using MessageDataType = std::vector<TDataType>;

    static inline MPI_Datatype DataType()
    {
        return MPIMessage<TDataType>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValues)
    {
        return rValues.data();
    }

    static inline const void* Buffer(const MessageDataType& rValues)
    {
        return rValues.data();
    }

    static inline int Size(const MessageDataType& rValues)
    {
        return (rValues.size() == 0 ? 0 : rValues.size() * MPIMessage<TDataType>::Size(rValues.front()));
    }

    static inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size())};
    }

    static inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rShape.size() == 1) << "Invalid shape for std::vector.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0]);
            return true;
        }

        return false;
    }

    static constexpr bool IsContiguous()
    {
        return std::is_arithmetic_v<TDataType>;
    }
};

template<> class MPIMessage<std::string>
{
public:
    using MessageDataType = std::string;

    static inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<char>::DataType();
    }

    static inline void* Buffer(MessageDataType& rValues)
    {
        return rValues.data();
    }

    static inline const void* Buffer(const MessageDataType& rValues)
    {
        return rValues.data();
    }

    static inline int Size(const MessageDataType& rValues)
    {
        return rValues.size();
    }

    static inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size())};
    }

    static inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rShape.size() == 1) << "Invalid shape for std::string.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0]);
            return true;
        }

        return false;
    }

    static constexpr bool IsContiguous()
    {
        return true;
    }
};

// template<> class MPIMessage<int>: public Internals::ValueMessage<int>, public Internals::MPIDataType<int> {};
// template<> class MPIMessage<unsigned int>: public Internals::ValueMessage<unsigned int>, public Internals::MPIDataType<unsigned int> {};
// template<> class MPIMessage<long unsigned int>: public Internals::ValueMessage<long unsigned int>, public Internals::MPIDataType<long unsigned int> {};
// template<> class MPIMessage<double>: public Internals::ValueMessage<double>, public Internals::MPIDataType<double> {};
// template<> class MPIMessage<bool>: public Internals::ValueMessage<bool>, public Internals::MPIDataType<bool> {};
// template<> class MPIMessage<char>: public Internals::ValueMessage<char>, public Internals::MPIDataType<char> {};
// template<> class MPIMessage<Flags::BlockType>: public Internals::ValueMessage<Flags::BlockType>, public Internals::MPIDataType<bool> {};

// template<> class MPIMessage<std::string>: public Internals::StringMessage, public Internals::MPIDataType<char> {};

// template<class ValueType> class MPIMessage< std::vector<ValueType> >: public Internals::VectorMessage<ValueType>, public Internals::MPIDataType<ValueType> {};
// template<class ValueType, std::size_t Dimension> class MPIMessage<array_1d<ValueType,Dimension>>: public Internals::ArrayMessage<array_1d<ValueType,Dimension>>, public Internals::MPIDataType<ValueType> {};
// template<> class MPIMessage<Vector>: public Internals::ArrayMessage<Vector>, public Internals::MPIDataType<double> {};
// template<> class MPIMessage<Matrix>: public Internals::MatrixMessage, public Internals::MPIDataType<double> {};

} // namespace Kratos
