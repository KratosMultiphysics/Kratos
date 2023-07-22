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
#include <variant>

// Project includes
#include "containers/array_1d.h"
#include "includes/data_communicator.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

using SupportedDataType = std::variant<
                        int,
                        double,
                        array_1d<double, 3>,
                        array_1d<double, 4>,
                        array_1d<double, 6>,
                        array_1d<double, 9>,
                        Vector,
                        Matrix>;

using IndicesType = std::vector<IndexType>;

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
        const IndicesType& rShape)
    {
        return false;
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

    static IndicesType Shape(const TDataType& rValue)
    {
        return {};
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
        FillToVector(rOutput.begin(), rValue);
    }

    template<class TIteratorType>
    static void FillToVector(
        TIteratorType Begin,
        const DataType& rValue)
    {
        *Begin = rValue;
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        FillFromVector(rValue, rOriginValues.begin(), rOriginValues.end());
    }

    template<class TIteratorType>
    static void FillFromVector(
        DataType& rValue,
        TIteratorType Begin,
        TIteratorType End)
    {
        rValue = *Begin;
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
        const IndicesType& rShape)
    {
        return false;
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

    static IndicesType Shape(const array_1d<TDataType, TSize>& rValue)
    {
        return {static_cast<IndexType>(TSize)};
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
        FillToVector(rOutput.begin(), rValue);
    }

    template<class TIteratorType>
    static void FillToVector(
        TIteratorType Begin,
        const DataType& rValue)
    {
        std::copy(rValue.begin(), rValue.end(), Begin);
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        FillFromVector(rValue, rOriginValues.begin(), rOriginValues.end());
    }

    template<class TIteratorType>
    static void FillFromVector(
        DataType& rValue,
        TIteratorType Begin,
        TIteratorType End)
    {
        std::copy(Begin, End, rValue.begin());
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
        const IndicesType& rShape)
    {
        KRATOS_ERROR_IF(rShape.size() != 1) << "Vectors should be reszed with shape having one dimension.";

        if (rValue.size() != rShape[0]) {
            rValue.resize(rShape[0], false);
            return true;
        } else {
            return false;
        }
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

    static IndicesType Shape(const DataType& rValue)
    {
        return {static_cast<IndexType>(rValue.size())};
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size();
    }

    static bool SynchronizeSize(
        DataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        KRATOS_TRY

        // first communicate all the sizes
        const int world_size = rDataCommunicator.Size();
        std::vector<unsigned int> local_size(1UL, rValue.size());
        std::vector<unsigned int> sizes(world_size);
        rDataCommunicator.AllGather(local_size, sizes);

        bool is_resized = false;
        if (rValue.size() == 0) {
            // now find non zero size to set this ranks empty vector
            for (const auto size : sizes) {
                if (size > 0) {
                    rValue.resize(size, false);
                    is_resized = true;
                    break;
                }
            }
        }

        // check all sizes in all ranks
        KRATOS_ERROR_IF_NOT(std::all_of(sizes.begin(), sizes.end(), [&rValue](const auto Size) {
                                    return Size == 0 || (Size > 0 && Size == rValue.size());
                            })) << "All the ranks should have the same vector size.\n";

        return is_resized;

        KRATOS_CATCH("");
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != Size(rValue)) {
            rOutput.resize(Size(rValue));
        }
        FillToVector(rOutput.begin(), rValue);
    }

    template<class TIteratorType>
    static void FillToVector(
        TIteratorType Begin,
        const DataType& rValue)
    {
        std::copy(rValue.begin(), rValue.end(), Begin);
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        FillFromVector(rValue, rOriginValues.begin(), rOriginValues.end());
    }

    template<class TIteratorType>
    static void FillFromVector(
        DataType& rValue,
        TIteratorType Begin,
        TIteratorType End)
    {
        std::copy(Begin, End, rValue.begin());
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
        const IndicesType& rShape)
    {
        KRATOS_ERROR_IF(rShape.size() != 2) << "Matrix should be reszed with shape having two dimensions.";

        if (rValue.size1() != rShape[0] || rValue.size2() != rShape[1]) {
            rValue.resize(rShape[0], rShape[1], false);
            return true;
        } else {
            return false;
        }
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

    static IndicesType Shape(const DataType& rValue)
    {
        return {static_cast<IndexType>(rValue.size1()), static_cast<IndexType>(rValue.size2())};
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size1() * rValue.size2();
    }

    static bool SynchronizeSize(
        DataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        KRATOS_TRY

        // first communicate all the sizes
        const unsigned int world_size = rDataCommunicator.Size();
        std::vector<std::vector<unsigned int>> sizes(world_size);
        sizes = rDataCommunicator.AllGatherv(std::vector<unsigned int>{static_cast<unsigned int>(rValue.size1()), static_cast<unsigned int>(rValue.size2())});

        bool is_resized = false;
        if (rValue.size1() == 0 || rValue.size2() == 0) {
            // now find non zero size to set this ranks empty vector
            for (const auto& size_info : sizes) {
                if (size_info[0] > 0 && size_info[1] > 0) {
                    rValue.resize(size_info[0], size_info[1], false);
                    is_resized = true;
                    break;
                }
            }
        }

        // check all sizes in all ranks
        KRATOS_ERROR_IF_NOT(std::all_of(sizes.begin(), sizes.end(), [&rValue](const auto SizeInfo) {
                                return SizeInfo[0] == 0 || SizeInfo[1] == 0 || (SizeInfo[0] > 0 && SizeInfo[0] == rValue.size1() && SizeInfo[1] == rValue.size2());
                            })) << "All the ranks should have the same matrix size.\n";

        return is_resized;

        KRATOS_CATCH("");
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != Size(rValue)) {
            rOutput.resize(Size(rValue));
        }
        FillToVector(rOutput.begin(), rValue);
    }

    template<class TIteratorType>
    static void FillToVector(
        TIteratorType Begin,
        const DataType& rValue)
    {
        std::copy(rValue.data().begin(), rValue.data().end(), Begin);
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        FillFromVector(rValue, rOriginValues.begin(), rOriginValues.end());
    }

    template<class TIteratorType>
    static void FillFromVector(
        DataType& rValue,
        TIteratorType Begin,
        TIteratorType End)
    {
        std::copy(Begin, End, rValue.data().begin());
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
        const IndicesType& rShape)
    {
        KRATOS_ERROR_IF(rShape.size() != 1) << "Vectors should be reszed with shape having one dimension.";

        if (rValue.size() != rShape[0]) {
            rValue.resize(rShape[0], false);
            return true;
        } else {
            return false;
        }
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

    static IndicesType Shape(const DataType& rValue)
    {
        return {static_cast<IndexType>(rValue.size())};
    }

    static IndexType Size(const DataType& rValue)
    {
        return rValue.size();
    }

    static bool SynchronizeSize(
        DataType& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        KRATOS_TRY

        // first communicate all the sizes
        const int world_size = rDataCommunicator.Size();
        std::vector<unsigned int> local_size(1UL, rValue.size());
        std::vector<unsigned int> sizes(world_size);
        rDataCommunicator.AllGather(local_size, sizes);

        bool is_resized = false;
        if (rValue.size() == 0) {
            // now find non zero size to set this ranks empty vector
            for (const auto size : sizes) {
                if (size > 0) {
                    rValue.resize(size, false);
                    is_resized = true;
                    break;
                }
            }
        }

        // check all sizes in all ranks
        KRATOS_ERROR_IF(std::all_of(sizes.begin(), sizes.end(), [&rValue](const auto Size) {
            return Size == 0 || (Size > 0 && Size == rValue.size());
        })) << "All the ranks should have the same vector size.\n";

        return is_resized;

        KRATOS_CATCH("");
    }

    static void FillToVector(
        VectorType& rOutput,
        const DataType& rValue)
    {
        if (rOutput.size() != Size(rValue)) {
            rOutput.resize(Size(rValue));
        }
        FillToVector(rOutput.begin(), rValue);
    }

    template<class TIteratorType>
    static void FillToVector(
        TIteratorType Begin,
        const DataType& rValue)
    {
        std::copy(rValue.begin(), rValue.end(), Begin);
    }

    static void FillFromVector(
        DataType& rValue,
        const VectorType& rOriginValues)
    {
        FillFromVector(rValue, rOriginValues.begin(), rOriginValues.end());
    }

    template<class TIteratorType>
    static void FillFromVector(
        DataType& rValue,
        TIteratorType Begin,
        TIteratorType End)
    {
        std::copy(Begin, End, rValue.begin());
    }
};

} // namespace Kratos