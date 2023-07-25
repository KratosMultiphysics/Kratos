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

    static constexpr bool HasDynamicMemoryAllocation = false;

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

    inline void Update(MessageDataType&)
    {
    }
};

template<class TDataType, std::size_t Dimension> class MPIMessage<array_1d<TDataType, Dimension>>
{
public:
    using MessageDataType = array_1d<TDataType, Dimension>;

    using PrimitiveDataType = TDataType;

    static constexpr bool HasDynamicMemoryAllocation = false;

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

    inline void Update(MessageDataType&)
    {
    }
};

template<> class MPIMessage<std::string>
{
public:
    using MessageDataType = std::string;

    using PrimitiveDataType = char;

    static constexpr bool HasDynamicMemoryAllocation = true;

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
        KRATOS_ERROR_IF_NOT(rShape.size() == 1) << "Invalid shape provided for std::string.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0]);
            return true;
        }

        return false;
    }

    inline void Update(MessageDataType&)
    {
    }
};

template<> class MPIMessage<Vector>
{
public:
    using MessageDataType = Vector;

    using PrimitiveDataType = double;

    static constexpr bool HasDynamicMemoryAllocation = true;

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
        KRATOS_ERROR_IF_NOT(rShape.size() == 1) << "Invalid shape provided for Vector.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0], false);
            return true;
        }

        return false;
    }

    inline void Update(MessageDataType&)
    {
    }
};

template<> class MPIMessage<Matrix>
{
public:
    using MessageDataType = Matrix;

    using PrimitiveDataType = double;

    static constexpr bool HasDynamicMemoryAllocation = true;

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
        return rValues.size1() * rValues.size2();
    }

    inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size1()), static_cast<unsigned int>(rValues.size2())};
    }

    inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() == 2) << "Invalid shape provided for Matrix.";

        if (rValues.size1() != rShape[0] || rValues.size2() != rShape[1]) {
            rValues.resize(rShape[0], rShape[1], false);
            return true;
        }

        return false;
    }

    inline void Update(MessageDataType&)
    {
    }
};

template<class TDataType> class MPIMessage<std::vector<TDataType>>
{
public:
    using MessageDataType = std::vector<TDataType>;

    using PrimitiveDataType = typename MPIMessage<TDataType>::PrimitiveDataType;

    static constexpr bool IsContiguous = std::is_same_v<TDataType, PrimitiveDataType>;

    static constexpr bool HasDynamicMemoryAllocation = true;

    inline MPI_Datatype DataType()
    {
        return Internals::MPIDataType<PrimitiveDataType>::DataType();
    }

    inline void* Buffer(MessageDataType& rValue)
    {
        if constexpr(IsContiguous) {
            return rValue.data();
        } else {
            if (rValue.size() > 0) {
                MPIMessage<TDataType> sub_mpi_message;
                const auto sub_data_type_size = sub_mpi_message.Size(rValue.front());

                const auto primitive_data_size = Size(rValue);
                if (mTempValues.size() != static_cast<unsigned int>(primitive_data_size)) {
                    mTempValues.resize(primitive_data_size);
                }

                unsigned int local_index = 0;
                for (unsigned int i = 0; i < rValue.size(); ++i) {
                    PrimitiveDataType* begin = static_cast<PrimitiveDataType*>(sub_mpi_message.Buffer(rValue[i]));
                    for (int j = 0; j < sub_data_type_size; ++j) {
                        mTempValues[local_index++] = *(begin++);
                    }
                }
            } else {
                mTempValues.resize(0);
            }

            return mTempValues.data();
        }
    }

    inline const void* Buffer(const MessageDataType& rValue)
    {
        if constexpr(IsContiguous) {
            return rValue.data();
        } else {
            if (rValue.size() > 0) {
                MPIMessage<TDataType> sub_mpi_message;
                const auto sub_data_type_size = sub_mpi_message.Size(rValue.front());

                const auto primitive_data_size = Size(rValue);
                if (mTempValues.size() != static_cast<unsigned int>(primitive_data_size)) {
                    mTempValues.resize(primitive_data_size);
                }

                unsigned int local_index = 0;
                for (unsigned int i = 0; i < rValue.size(); ++i) {
                    PrimitiveDataType const* begin = static_cast<PrimitiveDataType const*>(sub_mpi_message.Buffer(rValue[i]));
                    for (int j = 0; j < sub_data_type_size; ++j) {
                        mTempValues[local_index++] = *(begin++);
                    }
                }
            } else {
                mTempValues.resize(0);
            }

            return mTempValues.data();
        }
    }

    inline int Size(const MessageDataType& rValue)
    {
        return (rValue.size() == 0 ? 0 : rValue.size() * MPIMessage<TDataType>().Size(rValue.front()));
    }

    inline std::vector<unsigned int> Shape(const MessageDataType& rValues)
    {
        return {static_cast<unsigned int>(rValues.size())};
    }

    inline bool Resize(MessageDataType& rValues, const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() == 1) << "Invalid shape provided for std::vector.";

        if (rValues.size() != rShape[0]) {
            rValues.resize(rShape[0]);
            return true;
        }

        return false;
    }

    inline void Update(MessageDataType& rValue)
    {
        if constexpr(!IsContiguous) {
            KRATOS_ERROR_IF(rValue.size() == 0 && mTempValues.size() != 0)
                << "Size mismatch. [ rValue.size() = " << rValue.size()
                << ", mTempValues.size() = " << mTempValues.size() << " ].\n";

            KRATOS_ERROR_IF(rValue.size() != 0 && mTempValues.size() != static_cast<unsigned int>(Size(rValue)))
                << "Size mismatch. [ rValue.size() = " << rValue.size()
                << ", mTempValues.size() = " << mTempValues.size() << " ].\n";

            if (rValue.size() > 0) {
                MPIMessage<TDataType> sub_mpi_message;
                const auto sub_data_type_size = sub_mpi_message.Size(rValue.front());

                unsigned int local_index = 0;
                for (unsigned int i = 0; i < rValue.size(); ++i) {
                    PrimitiveDataType* begin = static_cast<PrimitiveDataType*>(sub_mpi_message.Buffer(rValue[i]));
                    for (int j = 0; j < sub_data_type_size; ++j) {
                        *(begin++) = mTempValues[local_index++];
                    }
                }
            }
        }
    }
private:
    std::vector<PrimitiveDataType> mTempValues;
};

} // namespace Kratos
