//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <algorithm>
#include <string>
#include <type_traits>
#include <vector>

// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos {

template<class TDataType> class DataTypeTraits
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = TDataType;

    using ValueType = TDataType;

    using PrimitiveType = TDataType;

    static constexpr bool HasContiguousPrimitiveData = true;

    static constexpr bool HasDynamicMemoryAllocation = false;

    ///@}
    ///@name Public static operations
    ///@{

    static inline unsigned int Size(const ContainerType&)
    {
        return 1;
    }

    static inline std::vector<unsigned int> Shape(const ContainerType&)
    {
        return {};
    }

    static inline bool Reshape(
        ContainerType&,
        const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF(rShape.size() != 0)
            << "Invalid shape given for primitive data type [ Expected shape = [], provided shape = "
            << rShape << " ].\n";
        return false;
    }

    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        return &rValue;
    }

    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        return &rValue;
    }

    inline static void FillToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rValue)
    {
        *pContiguousDataBegin = rValue;
    }

    inline static void FillFromContiguousData(
        ContainerType& rValue,
        PrimitiveType const * pContiguousDataBegin)
    {
        rValue = *pContiguousDataBegin;
    }

    ///@}
};

template<class TDataType, std::size_t Dimension> class DataTypeTraits<array_1d<TDataType, Dimension>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = array_1d<TDataType, Dimension>;

    using ValueType = TDataType;

    using PrimitiveType = typename DataTypeTraits<ValueType>::PrimitiveType;

    static constexpr bool HasContiguousPrimitiveData = std::is_same_v<PrimitiveType, ValueType>;

    static constexpr bool HasDynamicMemoryAllocation = DataTypeTraits<ValueType>::HasDynamicMemoryAllocation;

    ///@}
    ///@name Public static operations
    ///@{

    static inline unsigned int Size(const ContainerType& rContainer)
    {
        return Dimension * DataTypeTraits<ValueType>::Size(rContainer[0]);
    }

    static inline std::vector<unsigned int> Shape(const ContainerType& rContainer)
    {
        auto shape = DataTypeTraits<ValueType>::Shape(rContainer[0]);
        shape.insert(shape.begin(), Dimension);
        return shape;
    }

    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() >= 1 && rShape[0] == static_cast<unsigned int>(Dimension))
            << "Invalid shape given for array_1d data type [ Expected shape = [" << Dimension << "], provided shape = "
            << rShape << " ].\n";

        bool is_reshaped = false;

        if constexpr(DataTypeTraits<ValueType>::HasDynamicMemoryAllocation) {
            std::vector<unsigned int> sub_data_type_shape(rShape.size() - 1);
            std::copy(rShape.begin() + 1, rShape.end(), sub_data_type_shape.begin());
            std::for_each(rContainer.begin(), rContainer.end(), [&is_reshaped, &sub_data_type_shape](auto& rValue) { is_reshaped = DataTypeTraits<ValueType>::Reshape(rValue, sub_data_type_shape) || is_reshaped; });
        }

        return is_reshaped;
    }

    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data().begin();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data().begin();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static void FillToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        const auto stride = DataTypeTraits<ValueType>::Size(rContainer[0]);
        for (unsigned int i = 0; i < Dimension; ++i) {
            DataTypeTraits<ValueType>::FillToContiguousData(pContiguousDataBegin + i * stride, rContainer[i]);
        }
    }

    inline static void FillFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        const auto stride = DataTypeTraits<ValueType>::Size(rContainer[0]);
        for (unsigned int i = 0; i < Dimension; ++i) {
            DataTypeTraits<ValueType>::FillFromContiguousData(rContainer[i], pContiguousDataBegin + i * stride);
        }
    }

    ///@}
};

template<class TDataType> class DataTypeTraits<DenseVector<TDataType>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = DenseVector<TDataType>;

    using ValueType = TDataType;

    using PrimitiveType = typename DataTypeTraits<ValueType>::PrimitiveType;

    static constexpr bool HasContiguousPrimitiveData = std::is_same_v<PrimitiveType, ValueType>;

    static constexpr bool HasDynamicMemoryAllocation = true;

    ///@}
    ///@name Public static operations
    ///@{

    static inline unsigned int Size(const ContainerType& rValue)
    {
        return (rValue.size() == 0 ? 0 : rValue.size() * DataTypeTraits<ValueType>::Size(rValue[0]));
    }

    static inline std::vector<unsigned int> Shape(const ContainerType& rValue)
    {
        std::vector<unsigned int> shape;
        if (rValue.size() == 0) {
            shape = DataTypeTraits<ValueType>::Shape(ValueType{});
        } else {
            shape = DataTypeTraits<ValueType>::Shape(rValue[0]);
        }
        shape.insert(shape.begin(), rValue.size());
        return shape;
    }

    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() >= 1) << "Invalid shape given for DenseVector data type.";

        bool is_reshaped = false;

        if (rContainer.size() != rShape[0]) {
            rContainer.resize(rShape[0], false);
            is_reshaped = true;
        }

        if constexpr(DataTypeTraits<ValueType>::HasDynamicMemoryAllocation) {
            std::vector<unsigned int> sub_data_type_shape(rShape.size() - 1);
            std::copy(rShape.begin() + 1, rShape.end(), sub_data_type_shape.begin());
            std::for_each(rContainer.begin(), rContainer.end(), [&is_reshaped, &sub_data_type_shape](auto& rValue) { is_reshaped = DataTypeTraits<ValueType>::Reshape(rValue, sub_data_type_shape) || is_reshaped; });
        }

        return is_reshaped;
    }

    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data().begin();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data().begin();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static void FillToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        if (rContainer.size() != 0) {
            const auto stride = DataTypeTraits<ValueType>::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                DataTypeTraits<ValueType>::FillToContiguousData(pContiguousDataBegin + i * stride, rContainer[i]);
            }
        }
    }

    inline static void FillFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        if (rContainer.size() != 0) {
            const auto stride = DataTypeTraits<ValueType>::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                DataTypeTraits<ValueType>::FillFromContiguousData(rContainer[i], pContiguousDataBegin + i * stride);
            }
        }
    }

    ///@}
};

template<class TDataType> class DataTypeTraits<DenseMatrix<TDataType>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = DenseMatrix<TDataType>;

    using ValueType = TDataType;

    using PrimitiveType = typename DataTypeTraits<ValueType>::PrimitiveType;

    static constexpr bool HasContiguousPrimitiveData = std::is_same_v<PrimitiveType, ValueType>;

    static constexpr bool HasDynamicMemoryAllocation = true;

    ///@}
    ///@name Public static operations
    ///@{

    static inline unsigned int Size(const ContainerType& rValue)
    {
        return (rValue.size1() == 0 || rValue.size2() == 0 ? 0 : rValue.size1() * rValue.size2() * DataTypeTraits<ValueType>::Size(rValue.data()[0]));
    }

    static inline std::vector<unsigned int> Shape(const ContainerType& rValue)
    {
        std::vector<unsigned int> shape;
        if (rValue.size1() == 0 || rValue.size2() == 0) {
            shape = DataTypeTraits<ValueType>::Shape(ValueType{});
        } else {
            shape = DataTypeTraits<ValueType>::Shape(rValue.data()[0]);
        }
        shape.insert(shape.begin(), rValue.size2());
        shape.insert(shape.begin(), rValue.size1());
        return shape;
    }

    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() >= 2) << "Invalid shape given for DenseMatrix data type.";

        bool is_reshaped = false;

        if (rContainer.size1() != rShape[0] || rContainer.size2() != rShape[1]) {
            rContainer.resize(rShape[0], rShape[1], false);
            is_reshaped = true;
        }

        if constexpr(DataTypeTraits<ValueType>::HasDynamicMemoryAllocation) {
            std::vector<unsigned int> sub_data_type_shape(rShape.size() - 2);
            std::copy(rShape.begin() + 2, rShape.end(), sub_data_type_shape.begin());
            std::for_each(rContainer.data().begin(), rContainer.data().end(), [&is_reshaped, &sub_data_type_shape](auto& rValue) { is_reshaped = DataTypeTraits<ValueType>::Reshape(rValue, sub_data_type_shape) || is_reshaped; });
        }

        return is_reshaped;
    }

    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data().begin();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data().begin();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static void FillToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        if (rContainer.size1() != 0 && rContainer.size2() != 0) {
            const auto stride = DataTypeTraits<ValueType>::Size(rContainer(0, 0));
            for (unsigned int i = 0; i < rContainer.size1() * rContainer.size2(); ++i) {
                DataTypeTraits<ValueType>::FillToContiguousData(pContiguousDataBegin + i * stride, rContainer.data()[i]);
            }
        }
    }

    inline static void FillFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        if (rContainer.size1() != 0 && rContainer.size2() != 0) {
            const auto stride = DataTypeTraits<ValueType>::Size(rContainer(0, 0));
            for (unsigned int i = 0; i < rContainer.size1() * rContainer.size2(); ++i) {
                DataTypeTraits<ValueType>::FillFromContiguousData(rContainer.data()[i], pContiguousDataBegin + i * stride);
            }
        }
    }

    ///@}
};

template<> class DataTypeTraits<std::string>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = std::string;

    using ValueType = char;

    using PrimitiveType = char;

    static constexpr bool HasContiguousPrimitiveData = true;

    static constexpr bool HasDynamicMemoryAllocation = true;

    ///@}
    ///@name Public static operations
    ///@{

    static inline unsigned int Size(const ContainerType& rValue)
    {
        return rValue.size();
    }

    static inline std::vector<unsigned int> Shape(const ContainerType& rValue)
    {
        return {static_cast<unsigned int>(rValue.size())};
    }

    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() == 1) << "Invalid shape given for std::string data type.";

        bool is_reshaped = false;

        if (rContainer.size() != rShape[0]) {
            rContainer.resize(rShape[0], false);
            is_reshaped = true;
        }

        return is_reshaped;
    }

    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        return rValue.data();
    }

    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        return rValue.data();
    }

    inline static void FillToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        for (unsigned int i = 0; i < rContainer.size(); ++i) {
            *(pContiguousDataBegin++) = rContainer[i];
        }
    }

    inline static void FillFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        for (unsigned int i = 0; i < rContainer.size(); ++i) {
            rContainer[i] = *(pContiguousDataBegin++);
        }
    }

    ///@}
};

template<class TDataType> class DataTypeTraits<std::vector<TDataType>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = std::vector<TDataType>;

    using ValueType = TDataType;

    using PrimitiveType = typename DataTypeTraits<ValueType>::PrimitiveType;

    static constexpr bool HasContiguousPrimitiveData = std::is_same_v<PrimitiveType, ValueType>;

    static constexpr bool HasDynamicMemoryAllocation = true;

    ///@}
    ///@name Public static operations
    ///@{

    static inline unsigned int Size(const ContainerType& rValue)
    {
        return (rValue.size() == 0 ? 0 : rValue.size() * DataTypeTraits<ValueType>::Size(rValue[0]));
    }

    static inline std::vector<unsigned int> Shape(const ContainerType& rValue)
    {
        std::vector<unsigned int> shape;
        if (rValue.size() == 0) {
            shape = DataTypeTraits<ValueType>::Shape(ValueType{});
        } else {
            shape = DataTypeTraits<ValueType>::Shape(rValue[0]);
        }
        shape.insert(shape.begin(), rValue.size());
        return shape;
    }

    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<unsigned int>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() >= 1) << "Invalid shape given for std::vector data type.";

        bool is_reshaped = false;

        if (rContainer.size() != rShape[0]) {
            rContainer.resize(rShape[0]);
            is_reshaped = true;
        }

        if constexpr(DataTypeTraits<ValueType>::HasDynamicMemoryAllocation) {
            std::vector<unsigned int> sub_data_type_shape(rShape.size() - 1);
            std::copy(rShape.begin() + 1, rShape.end(), sub_data_type_shape.begin());
            std::for_each(rContainer.begin(), rContainer.end(), [&is_reshaped, &sub_data_type_shape](auto& rValue) { is_reshaped = DataTypeTraits<ValueType>::Reshape(rValue, sub_data_type_shape) || is_reshaped; });
        }

        return is_reshaped;
    }

    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return rValue.data();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    inline static void FillToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        if (rContainer.size() != 0) {
            const auto stride = DataTypeTraits<ValueType>::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                DataTypeTraits<ValueType>::FillToContiguousData(pContiguousDataBegin + i * stride, rContainer[i]);
            }
        }
    }

    inline static void FillFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        if (rContainer.size() != 0) {
            const auto stride = DataTypeTraits<ValueType>::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                DataTypeTraits<ValueType>::FillFromContiguousData(rContainer[i], pContiguousDataBegin + i * stride);
            }
        }
    }

    ///@}
};

template<class TDataType>
class DataBuffer
{
public:
    ///@name Type definitions
    ///@{

    using BufferedDataTypeTrait = DataTypeTraits<TDataType>;

    using PrimitiveType = typename BufferedDataTypeTrait::PrimitiveType;

    static constexpr bool HasContiguousPrimitiveData = BufferedDataTypeTrait::HasContiguousPrimitiveData;

    ///@}
    ///@name Public operations
    ///@{

    inline PrimitiveType const * GetData(const TDataType& rContainer) const
    {
        if constexpr(HasContiguousPrimitiveData) {
            return BufferedDataTypeTrait::GetContiguousData(rContainer);
        } else {
            const auto size = BufferedDataTypeTrait::Size(rContainer);
            if (mData.size() != size) {
                mData.resize(size);
            }
            BufferedDataTypeTrait::FillToContiguousData(mData.begin(), rContainer);
            return mData.begin();
        }
    }

    inline PrimitiveType * GetData(TDataType& rContainer)
    {
        if constexpr(HasContiguousPrimitiveData) {
            return BufferedDataTypeTrait::GetContiguousData(rContainer);
        } else {
            const auto size = BufferedDataTypeTrait::Size(rContainer);
            if (mData.size() != size) {
                mData.resize(size);
            }
            BufferedDataTypeTrait::FillToContiguousData(mData.data(), rContainer);
            return mData.data();
        }
    }

    void UpdateValues(TDataType& rContainer) const
    {
        if constexpr(!HasContiguousPrimitiveData) {
            KRATOS_ERROR_IF_NOT(BufferedDataTypeTrait::Size(rContainer) == static_cast<unsigned int>(mData.size()))
                << "Size mismatch [ rContainer flat size = " << BufferedDataTypeTrait::Size(rContainer)
                << ", buffered data size = " << mData.size() << " ].\n";
            BufferedDataTypeTrait::FillFromContiguousData(rContainer, mData.data());
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::conditional_t<HasContiguousPrimitiveData, void *, std::vector<PrimitiveType>> mData;

    ///@}
};

}; // namespace Kratos
