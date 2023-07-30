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
#include <algorithm>
#include <type_traits>

// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

template<class TDataType>
struct DataTypeTraits
{
    using DataType = TDataType;

    using RawDataType = TDataType;

    static constexpr bool HasDynamicMemoryAllocation = false;

    static void Initialize(
        TDataType& rValue,
        const TDataType InitializationValue)
    {
        static_assert(std::is_arithmetic_v<TDataType>, "Unsupported data type.");
        rValue = InitializationValue;
    }

    static RawDataType& GetComponent(
        TDataType& rValue,
        const IndexType Index)
    {
        return rValue;
    }

    static RawDataType GetComponent(
        const TDataType& rValue,
        const IndexType Index)
    {
        return rValue;
    }

    static IndexType Size(const TDataType& rValue)
    {
        return 1;
    }
};

template<class TDataType, std::size_t TSize>
struct DataTypeTraits<array_1d<TDataType, TSize>>
{
    using DataType = array_1d<TDataType, TSize>;

    using RawDataType = TDataType;

    static constexpr bool HasDynamicMemoryAllocation = false;

    static void Initialize(
        DataType& rValue,
        const TDataType InitializationValue)
    {
        std::fill(rValue.begin(), rValue.end(), InitializationValue);
    }

    static RawDataType& GetComponent(
        DataType& rValue,
        const IndexType Index)
    {
        return rValue[Index];
    }

    static RawDataType GetComponent(
        const DataType& rValue,
        const IndexType Index)
    {
        return rValue[Index];
    }

    static IndexType Size(const DataType& rValue)
    {
        return TSize;
    }
};

template<>
struct DataTypeTraits<Vector>
{
    using DataType = Vector;

    using RawDataType = double;

    static constexpr bool HasDynamicMemoryAllocation = true;

    static void Initialize(
        DataType& rValue,
        const double InitializationValue)
    {
        std::fill(rValue.begin(), rValue.end(), InitializationValue);
    }

    static RawDataType& GetComponent(
        DataType& rValue,
        const IndexType Index)
    {
        return rValue[Index];
    }

    static RawDataType GetComponent(
        const DataType& rValue,
        const IndexType Index)
    {
        return rValue[Index];
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size();
    }
};

template<>
struct DataTypeTraits<Matrix>
{
    using DataType = Matrix;

    using RawDataType = double;

    static constexpr bool HasDynamicMemoryAllocation = true;

    static void Initialize(
        DataType& rValue,
        const double InitializationValue)
    {
        std::fill(rValue.data().begin(), rValue.data().end(), InitializationValue);
    }

    static RawDataType& GetComponent(
        DataType& rValue,
        const IndexType Index)
    {
        return rValue.data()[Index];
    }

    static RawDataType GetComponent(
        const DataType& rValue,
        const IndexType Index)
    {
        return rValue.data()[Index];
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size1() * rValue.size2();
    }
};

template<class TDataType>
struct DataTypeTraits<std::vector<TDataType>>
{
    using DataType = std::vector<TDataType>;

    using RawDataType = TDataType;

    static constexpr bool HasDynamicMemoryAllocation = true;

    static void Initialize(
        DataType& rValue,
        const TDataType InitializationValue)
    {
        std::fill(rValue.begin(), rValue.end(), InitializationValue);
    }

    static RawDataType& GetComponent(
        DataType& rValue,
        const IndexType Index)
    {
        return rValue[Index];
    }

    static RawDataType GetComponent(
        const DataType& rValue,
        const IndexType Index)
    {
        return rValue[Index];
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size();
    }
};

} // namespace Kratos