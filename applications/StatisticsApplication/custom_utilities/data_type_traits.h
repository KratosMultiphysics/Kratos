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
#include "includes/data_communicator.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

template<class TDataType>
struct DataTypeTraits
{
    using DataType = TDataType;

    using RawDataType = TDataType;

    using VectorType = std::vector<RawDataType>;

    static void Initialize(
        TDataType& rValue,
        const TDataType InitializationValue)
    {
        static_assert(std::is_arithmetic_v<TDataType>, "Unsupported data type.");
        rValue = InitializationValue;
    }

    static bool Resize(
        TDataType& rValue,
        const TDataType rTarget)
    {
        return false;
    }

    static IndexType Size(const TDataType& rValue)
    {
        return 1;
    }

    static bool SynchronizeSize(
        TDataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        return false;
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != 1) {
            rOutput.resize(1);
        }
        rOutput[0] = rValue;
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        rValue = rOriginValues[0];
    }
};

template<class TDataType, std::size_t TSize>
struct DataTypeTraits<array_1d<TDataType, TSize>>
{
    using DataType = array_1d<TDataType, TSize>;

    using RawDataType = TDataType;

    using VectorType = std::vector<RawDataType>;

    static void Initialize(
        DataType& rValue,
        const TDataType InitializationValue)
    {
        std::fill(rValue.begin(), rValue.end(), InitializationValue);
    }

    static bool Resize(
        DataType& rValue,
        const DataType& rTarget)
    {
        return false;
    }

    static IndexType Size(const DataType& rValue)
    {
        return TSize;
    }

    static bool SynchronizeSize(
        DataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        return false;
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != TSize) {
            rOutput.resize(TSize);
        }
        std::copy(rValue.begin(), rValue.end(), rOutput.begin());
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        std::copy(rOriginValues.begin(), rOriginValues.end(), rValue.begin());
    }
};

template<>
struct DataTypeTraits<Vector>
{
    using DataType = Vector;

    using RawDataType = double;

    using VectorType = std::vector<RawDataType>;

    static void Initialize(
        DataType& rValue,
        const double InitializationValue)
    {
        std::fill(rValue.begin(), rValue.end(), InitializationValue);
    }

    static bool Resize(
        DataType& rValue,
        const DataType& rTarget)
    {
        if (rValue.size() != rTarget.size()) {
            rValue.resize(rTarget.size(), false);
            return true;
        } else {
            return false;
        }
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size();
    }

    static bool SynchronizeSize(
        DataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        // first communicate all the sizes
        const int world_size = rDataCommunicator.Size();
        std::vector<int> local_size(1UL, rValue.size());
        std::vector<int> sizes(world_size);
        rDataCommunicator.AllGather(local_size, sizes);

        if (rValue.size() == 0) {
            // now find non zero size to set this ranks empty vector
            for (const auto size : sizes) {
                if (size > 0) {
                    rValue.resize(size, false);
                    return true;
                }
            }
        }

        return false;
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != Size(rValue)) {
            rOutput.resize(Size(rValue));
        }
        std::copy(rValue.begin(), rValue.end(), rOutput.begin());
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        std::copy(rOriginValues.begin(), rOriginValues.end(), rValue.begin());
    }
};

template<>
struct DataTypeTraits<Matrix>
{
    using DataType = Matrix;

    using RawDataType = double;

    using VectorType = std::vector<RawDataType>;

    static void Initialize(
        DataType& rValue,
        const double InitializationValue)
    {
        std::fill(rValue.data().begin(), rValue.data().end(), InitializationValue);
    }

    static bool Resize(
        DataType& rValue,
        const DataType& rTarget)
    {
        if (rValue.size1() != rTarget.size1() || rValue.size2() != rTarget.size2()) {
            rValue.resize(rTarget.size1(), rTarget.size2(), false);
            return true;
        } else {
            return false;
        }
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size1() * rValue.size2();
    }

    static bool SynchronizeSize(
        DataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        // first communicate all the sizes
        const int world_size = rDataCommunicator.Size();
        std::vector<std::vector<int>> sizes(world_size);
        sizes = rDataCommunicator.AllGatherv(std::vector<int>{static_cast<int>(rValue.size1()), static_cast<int>(rValue.size2())});

        if (rValue.size1() == 0 || rValue.size2() == 0) {
            // now find non zero size to set this ranks empty vector
            for (const auto& size_info : sizes) {
                if (size_info[0] > 0 && size_info[1] > 0) {
                    rValue.resize(size_info[0], size_info[1], false);
                    return true;
                }
            }
        }

        return false;
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != Size(rValue)) {
            rOutput.resize(Size(rValue));
        }
        std::copy(rValue.data().begin(), rValue.data().end(), rOutput.begin());
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        std::copy(rOriginValues.begin(), rOriginValues.end(), rValue.data().begin());
    }
};

template<class TDataType>
struct DataTypeTraits<std::vector<TDataType>>
{
    using DataType = std::vector<TDataType>;

    using RawDataType = TDataType;

    using VectorType = std::vector<RawDataType>;

    static void Initialize(
        DataType& rValue,
        const TDataType InitializationValue)
    {
        std::fill(rValue.begin(), rValue.end(), InitializationValue);
    }

    static bool Resize(
        DataType& rValue,
        const DataType& rTarget)
    {
        if (rValue.size() != rTarget.size()) {
            rValue.resize(rTarget.size(), false);
            return true;
        } else {
            return false;
        }
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size();
    }

    static bool SynchronizeSize(
        DataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        // first communicate all the sizes
        const int world_size = rDataCommunicator.Size();
        std::vector<int> local_size(1UL, rValue.size());
        std::vector<int> sizes(world_size);
        rDataCommunicator.AllGather(local_size, sizes);

        if (rValue.size() == 0) {
            // now find non zero size to set this ranks empty vector
            for (const auto size : sizes) {
                if (size > 0) {
                    rValue.resize(size, false);
                    return true;
                }
            }
        }

        return false;
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != Size(rValue)) {
            rOutput.resize(Size(rValue));
        }
        std::copy(rValue.begin(), rValue.end(), rOutput.begin());
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        std::copy(rOriginValues.begin(), rOriginValues.end(), rValue.begin());
    }
};

} // namespace Kratos