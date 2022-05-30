//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Jordi Cotela
//                   Philipp Bucher (https://github.com/philbucher)
//
// Ported from "kratos/mpi/includes/mpi_message.h"
//

#ifndef CO_SIM_IO_MPI_MESSAGE_INCLUDED
#define CO_SIM_IO_MPI_MESSAGE_INCLUDED

// System includes
#include <string>
#include <array>
#include <vector>

// External includes
#include "mpi.h"

namespace CoSimIO {
namespace Internals {

template<class TDataType> struct MPIDataType;

template<> struct MPIDataType<int>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_INT;
    }
    static constexpr int LengthPerObject = 1;
};

template<> struct MPIDataType<unsigned int>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_UNSIGNED;
    }
    static constexpr int LengthPerObject = 1;
};

template<> struct MPIDataType<long unsigned int>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_UNSIGNED_LONG;
    }
    static constexpr int LengthPerObject = 1;
};

template<> struct MPIDataType<double>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_DOUBLE;
    }
    static constexpr int LengthPerObject = 1;
};

template<> struct MPIDataType<std::string>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_CHAR;
    }
    static constexpr int LengthPerObject = 1;
};

template<> struct MPIDataType<bool>
{
    static inline MPI_Datatype DataType()
    {
        return MPI_C_BOOL;
    }
    static constexpr int LengthPerObject = 1;
};

template<class TDataType> class ValueMessage
{
public:
    static inline void* Buffer(TDataType& rValue)
    {
        return &rValue;
    }

    static inline const void* Buffer(const TDataType& rValue)
    {
        return &rValue;
    }

    static inline int Size(const TDataType&)
    {
        return 1;
    }
};

template<class TDataType> class VectorMessage
{
public:
    static inline void* Buffer(std::vector<TDataType>& rValue)
    {
        return rValue.data();
    }

    static inline const void* Buffer(const std::vector<TDataType>& rValue)
    {
        return rValue.data();
    }

    static inline int Size(const std::vector<TDataType>& rValue)
    {
        return rValue.size() * MPIDataType<TDataType>::LengthPerObject;
    }
};

template<class TDataType, std::size_t Dimension> class ArrayMessage
{
public:
    static inline void* Buffer(std::array<TDataType,Dimension>& rValues)
    {
        return rValues.data();
    }

    static inline const void* Buffer(const std::array<TDataType,Dimension>& rValues)
    {
        return rValues.data();
    }

    static inline int Size(const std::array<TDataType,Dimension>& rValues)
    {
        return Dimension * MPIDataType<TDataType>::LengthPerObject;
    }
};

class StringMessage
{
public:
    static inline void* Buffer(std::string& rValue)
    {
        // Note: this uses the fact that the C++11 standard defines std::strings to
        // be contiguous in memory to perform MPI communication (based on char*) in place.
        // In older C++, this cannot be expected to be always the case, so a copy of the
        // string would be required.
        return const_cast<char *>(rValue.data());
    }

    static inline const void* Buffer(const std::string& rValue)
    {
        return rValue.data();
    }

    static inline int Size(const std::string& rValues)
    {
        return rValues.size();
    }
};

} // namespace Internal

template<class TDataType> class MPIMessage;

template<> class MPIMessage<int>: public Internals::ValueMessage<int>, public Internals::MPIDataType<int> {};
template<> class MPIMessage<unsigned int>: public Internals::ValueMessage<unsigned int>, public Internals::MPIDataType<unsigned int> {};
template<> class MPIMessage<long unsigned int>: public Internals::ValueMessage<long unsigned int>, public Internals::MPIDataType<long unsigned int> {};
template<> class MPIMessage<double>: public Internals::ValueMessage<double>, public Internals::MPIDataType<double> {};
template<> class MPIMessage<bool>: public Internals::ValueMessage<bool>, public Internals::MPIDataType<bool> {};

template<> class MPIMessage<std::string>: public Internals::StringMessage, public Internals::MPIDataType<std::string> {};

template<class ValueType> class MPIMessage< std::vector<ValueType> >: public Internals::VectorMessage<ValueType>, public Internals::MPIDataType<ValueType> {};
template<class ValueType, std::size_t Dimension> class MPIMessage<std::array<ValueType,Dimension>>: public Internals::ArrayMessage<ValueType,Dimension>, public Internals::MPIDataType<ValueType> {};

} // namespace CoSimIO

#endif // CO_SIM_IO_MPI_MESSAGE_INCLUDED