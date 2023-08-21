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

/**
 * @brief Generic data type traits class for arithmetic types.
 *
 * @tparam TDataType        Arithmetic data type
 */
template<class TDataType> class DataTypeTraits
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = TDataType;

    using ValueType = TDataType;

    using PrimitiveType = TDataType;

    static constexpr bool IsContiguous = true;

    static constexpr bool IsDynamic = false;

    static constexpr unsigned int Dimension = 0;

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Returns whther the given @ref TCheckIndex is dynamic.
     *
     * This method returns true if the @ref TCheckIndex of the dimension
     * corresponds to a dynamic type. If not, this returns false.
     *
     * @tparam TCheckIndex          User input dimension index.
     * @return true                 If the input dimension index corresponds to dynamic data type.
     * @return false                If the input dimension index corresponds to static data type.
     */
    template<unsigned int TCheckIndex>
    static constexpr bool IsDimensionDynamic()
    {
        return IsDimensionDynamicImpl<TCheckIndex, 0>();
    }

    /**
     * @brief Returns the size of the value.
     *
     * @return TIndexType Size of the value.
     */
    template<class TIndexType = unsigned int>
    static inline TIndexType Size(const ContainerType&)
    {
        return 1;
    }

    /**
     * @brief Returns a vector with the shape of the value.
     *
     * Scalars have the shape of [] (an empty vector).
     *
     * @return std::vector<TIndexType>    Returns an empty vector as the shape.
     */
    template<class TIndexType = unsigned int>
    static inline std::vector<TIndexType> Shape(const ContainerType&)
    {
        return {};
    }

    /**
     * @brief Fills the given array with the shape values in each dimension.
     *
     * @throws If the array is not of the required size.
     *
     * @tparam int              Type of the shape values.
     * @param pShapeBegin       Begin of the shape array.
     * @param pShapeEnd         End of the shape array.
     */
    template<class TIndexType = unsigned int>
    static inline void Shape(
        const ContainerType&,
        TIndexType* pShapeBegin,
        TIndexType* pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) == 0)
            << "Invalid dimensions given to fill for primitive data type [ Expected dimension == 0, provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";
    }

    /**
     * @brief Reshapes the value.
     *
     * Scalar values are not reshaped. Hence, always returns false.
     * Return value is true if a change has been made to the passed value
     * due to prescribed shape.
     *
     * @param rShape        Required shape.
     * @return true         If the passed value is changed due to the required shape.
     * @return false        If the passed value is same as the given shape.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rValue,
        const std::vector<TIndexType>& rShape)
    {
        return Reshape(rValue, rShape.data(), rShape.data() + rShape.size());
    }

    /**
     * @brief Reshapes the value.
     *
     * Scalar values are not reshaped. Hence, always returns false.
     * Return value is true if a change has been made to the passed value
     * due to prescribed shape.
     *
     * @param rShape        Required shape.
     * @return true         If the passed value is changed due to the required shape.
     * @return false        If the passed value is same as the given shape.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType&,
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) == 0)
            << "Invalid shape/dimension given for primitive data type [ Expected shape = [], provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";
        return false;
    }

    /**
     * @brief Get the contiguous data pointer.
     *
     * @param rValue                    Value for which the contiguous data pointer is returned.
     * @return PrimitiveType const*     Contiguous data array pointer.
     */
    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        return &rValue;
    }

    /**
     * @brief Get the contiguous data pointer.
     *
     * @param rValue                    Value for which the contiguous data pointer is returned.
     * @return PrimitiveType const*     Contiguous data array pointer.
     */
    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        return &rValue;
    }

    /**
     * @brief Copies the given Value to contiguous array.
     *
     * This method copies content of the @ref rValue to a contiguous array
     * at @ref pContiguousDataBegin. The array pointed by @ref pContiguousDataBegin
     * should be sized correctly.
     *
     * @warning This may seg fault if the @ref pContiguousDataBegin is not correctly sized.
     *
     * @param pContiguousDataBegin      Starting value pointer of the contiguous array.
     * @param rValue                    Value to be copied to contiguous array.
     */
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rValue)
    {
        *pContiguousDataBegin = rValue;
    }

    /**
     * @brief Copies data from contiguous array to rValue.
     *
     * This method copies data from contiguous array to the passed @ref rValue.
     * The contiguous array is given by the @refpContiguousDataBegin.
     *
     * @warning This may seg fault if the @ref pContiguousDataBegin is not correctly sized.
     *
     * @param rValue                Output value which contains copied values.
     * @param pContiguousDataBegin  Array to copy data from.
     */
    inline static void CopyFromContiguousData(
        ContainerType& rValue,
        PrimitiveType const * pContiguousDataBegin)
    {
        rValue = *pContiguousDataBegin;
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TCheckIndex, unsigned int TCurrentIndex = 0>
    static constexpr bool IsDimensionDynamicImpl()
    {
        if constexpr(TCheckIndex == 0) {
            return false;
        } else {
            static_assert(TCheckIndex != TCheckIndex, "Invalid dimension index.");
        }
    }

    ///@}
    ///@name Friends
    ///@{

    template<class T>
    friend class DataTypeTraits;

    ///@}
};

/**
 * @brief Data type traits for array_1d data types
 *
 * @tparam TDataType    Data type of array_1d
 * @tparam TSize    Size of array_1d
 */
template<class TDataType, std::size_t TSize> class DataTypeTraits<array_1d<TDataType, TSize>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = array_1d<TDataType, TSize>;

    using ValueType = TDataType;

    using ValueTraits = DataTypeTraits<ValueType>;

    using PrimitiveType = typename ValueTraits::PrimitiveType;

    static constexpr bool IsDynamic = ValueTraits::IsDynamic;

    // boost ublas makes the underlying data structure contiguous for
    // any ValueType which are not dynamic recursively.
    static constexpr bool IsContiguous = !IsDynamic;

    static constexpr unsigned int Dimension = ValueTraits::Dimension + 1;

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Returns whther the given @ref TCheckIndex is dynamic.
     *
     * This method returns true if the @ref TCheckIndex of the dimension
     * corresponds to a dynamic type. If not, this returns false.
     *
     * @tparam TCheckIndex          User input dimension index.
     * @return true                 If the input dimension index corresponds to dynamic data type.
     * @return false                If the input dimension index corresponds to static data type.
     */
    template<unsigned int TCheckIndex>
    static constexpr bool IsDimensionDynamic()
    {
        return IsDimensionDynamicImpl<TCheckIndex, 0>();
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @ref PrimitiveType values contained in
     * the @ref rContainer recursively.
     *
     * @param rContainer        The container to calculate the size.
     * @return TIndexType       Number of primitive type values in rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline TIndexType Size(const ContainerType& rContainer)
    {
        if constexpr(TSize == 0) {
            return 0;
        } else {
            return TSize * ValueTraits::template Size<TIndexType>(rContainer[0]);
        }
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @ref rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline std::vector<TIndexType> Shape(const ContainerType& rContainer)
    {
        std::vector<TIndexType> shape(Dimension);
        Shape(rContainer, shape.data(), shape.data() + Dimension);
        return shape;
    }

    /**
     * @brief Fills the given array with the shape values in each dimension.
     *
     * @throws If the array is not of the required size.
     *
     * @tparam int              Type of the shape values.
     * @param pShapeBegin       Begin of the shape array.
     * @param pShapeEnd         End of the shape array.
     */
    template<class TIndexType = unsigned int>
    static inline void Shape(
        const ContainerType& rContainer,
        TIndexType* pShapeBegin,
        TIndexType* pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 1)
            << "Invalid dimensions given to fill for primitive data type [ Expected dimension >= 1, provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        if constexpr(TSize > 0) {
            ValueTraits::Shape(rContainer[0], pShapeBegin + 1, pShapeEnd);
        } else {
            ValueTraits::Shape(ValueType{}, pShapeBegin + 1, pShapeEnd);
        }
        pShapeBegin[0] = TSize;
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the @ref rShape
     * recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given
     * @ref rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<TIndexType>& rShape)
    {
        return Reshape(rContainer, rShape.data(), rShape.data() + rShape.size());
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the shape given by @ref pShapeBegin
     * and @ref pShapeEnd recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given shape
     * represented by @ref pShapeBegin and @ref pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 1 && *pShapeBegin == static_cast<TIndexType>(TSize))
            << "Invalid shape/dimension given for array_1d data type [ Expected shape = " << Shape(rContainer) << ", provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        bool is_reshaped = false;

        if constexpr(ValueTraits::IsDynamic) {
            std::for_each(rContainer.begin(), rContainer.end(), [&is_reshaped, pShapeBegin, pShapeEnd](auto& rValue) {
                is_reshaped = ValueTraits::Reshape(rValue, pShapeBegin + 1, pShapeEnd) || is_reshaped;
            });
        }

        return is_reshaped;
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            if constexpr(std::is_same_v<PrimitiveType, ValueType>) {
                return rValue.data().data();
            } else {
                // since the underlying data structure for recusive static data types
                // is contiguous in ublas types, we can do the following to get the
                // contiguous array.
                if constexpr(TSize > 0) {
                    return reinterpret_cast<PrimitiveType const *>(&rValue[0]);
                } else {
                    // not returning nullptr so, the return value can be subjected to
                    // arithmetic operations
                    return 0;
                }
            }
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            if constexpr(std::is_same_v<PrimitiveType, ValueType>) {
                return rValue.data().data();
            } else {
                // since the underlying data structure for recusive static data types
                // is contiguous in ublas types, we can do the following to get the
                // contiguous array.
                if constexpr(TSize > 0) {
                    return reinterpret_cast<PrimitiveType*>(&rValue[0]);
                } else {
                    // not returning nullptr so, the return value can be subjected to
                    // arithmetic operations
                    return 0;
                }
            }
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "This should be only called if the rValue is contiguous only.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @ref rContainer recursively to the
     * provided @ref pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Contaienr to copy data to the contiguous array.
     */
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        const auto stride = ValueTraits::Size(rContainer[0]);
        for (unsigned int i = 0; i < TSize; ++i) {
            ValueTraits::CopyToContiguousData(pContiguousDataBegin + i * stride, rContainer[i]);
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @ref pContiguousDataBegin pointer
     * to the @ref rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     */
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        const auto stride = ValueTraits::Size(rContainer[0]);
        for (unsigned int i = 0; i < TSize; ++i) {
            ValueTraits::CopyFromContiguousData(rContainer[i], pContiguousDataBegin + i * stride);
        }
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TCheckIndex, unsigned int TCurrentIndex>
    static constexpr bool IsDimensionDynamicImpl()
    {
        if constexpr(TCheckIndex == TCurrentIndex) {
            return false;
        } else {
            return ValueTraits::template IsDimensionDynamicImpl<TCheckIndex, TCurrentIndex + 1>();
        }
    }

    ///@}
    ///@name Friends
    ///@{

    template<class T>
    friend class DataTypeTraits;

    ///@}
};

/**
 * @brief Data type traits for DenseVector data types
 *
 * @tparam TDataType    Data type of DenseVector
 */
template<class TDataType> class DataTypeTraits<DenseVector<TDataType>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = DenseVector<TDataType>;

    using ValueType = TDataType;

    using ValueTraits = DataTypeTraits<ValueType>;

    using PrimitiveType = typename ValueTraits::PrimitiveType;

    static constexpr bool IsDynamic = true;

    // boost ublas makes the underlying data structure contiguous for
    // any ValueType which are not dynamic recursively.
    static constexpr bool IsContiguous = !ValueTraits::IsDynamic;

    static constexpr unsigned int Dimension = ValueTraits::Dimension + 1;

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Returns whther the given @ref TCheckIndex is dynamic.
     *
     * This method returns true if the @ref TCheckIndex of the dimension
     * corresponds to a dynamic type. If not, this returns false.
     *
     * @tparam TCheckIndex          User input dimension index.
     * @return true                 If the input dimension index corresponds to dynamic data type.
     * @return false                If the input dimension index corresponds to static data type.
     */
    template<unsigned int TCheckIndex>
    static constexpr bool IsDimensionDynamic()
    {
        return IsDimensionDynamicImpl<TCheckIndex, 0>();
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @ref PrimitiveType values contained in
     * the @ref rContainer recursively.
     *
     * @param rContainer        The container to calculate the size.
     * @return TIndexType       Number of primitive type values in rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline TIndexType Size(const ContainerType& rValue)
    {
        return (rValue.empty() ? 0 : rValue.size() * ValueTraits::template Size<TIndexType>(rValue[0]));
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @ref rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline std::vector<TIndexType> Shape(const ContainerType& rContainer)
    {
        std::vector<TIndexType> shape(Dimension);
        Shape(rContainer, shape.data(), shape.data() + Dimension);
        return shape;
    }

    /**
     * @brief Fills the given array with the shape values in each dimension.
     *
     * @throws If the array is not of the required size.
     *
     * @tparam int              Type of the shape values.
     * @param pShapeBegin       Begin of the shape array.
     * @param pShapeEnd         End of the shape array.
     */
    template<class TIndexType = unsigned int>
    static inline void Shape(
        const ContainerType& rContainer,
        TIndexType* pShapeBegin,
        TIndexType* pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 1)
            << "Invalid dimensions given to fill for primitive data type [ Expected dimension >= 1, provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        if (rContainer.empty()) {
            ValueTraits::Shape(ValueType{}, pShapeBegin + 1, pShapeEnd);
        } else {
            ValueTraits::Shape(rContainer[0], pShapeBegin + 1, pShapeEnd);
        }
        pShapeBegin[0] = rContainer.size();
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the @ref rShape
     * recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given
     * @ref rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<TIndexType>& rShape)
    {
        return Reshape(rContainer, rShape.data(), rShape.data() + rShape.size());
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the shape given by @ref pShapeBegin
     * and @ref pShapeEnd recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given shape
     * represented by @ref pShapeBegin and @ref pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 1)
            << "Invalid shape/dimension given for DenseVector data type [ Expected = " << Shape(rContainer) << ", provided = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        bool is_reshaped = false;

        if (rContainer.size() != *pShapeBegin) {
            rContainer.resize(*pShapeBegin, false);
            is_reshaped = true;
        }

        if constexpr(ValueTraits::IsDynamic) {
            std::for_each(rContainer.begin(), rContainer.end(), [&is_reshaped, pShapeBegin, pShapeEnd](auto& rValue) {
                is_reshaped = ValueTraits::Reshape(rValue, pShapeBegin + 1, pShapeEnd) || is_reshaped;
            });
        }

        return is_reshaped;
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            if constexpr(std::is_same_v<PrimitiveType, ValueType>) {
                return rValue.data().begin();
            } else {
                // since the underlying data structure for recusive static data types
                // is contiguous in ublas types, we can do the following to get the
                // contiguous array.
                if (rValue.size() > 0) {
                    return reinterpret_cast<PrimitiveType const *>(&rValue[0]);
                } else {
                    // not returning nullptr so, the return value can be subjected to
                    // arithmetic operations
                    return 0;
                }
            }
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            if constexpr(std::is_same_v<PrimitiveType, ValueType>) {
                return rValue.data().begin();
            } else {
                // since the underlying data structure for recusive static data types
                // is contiguous in ublas types, we can do the following to get the
                // contiguous array.
                if (rValue.size() > 0) {
                    return reinterpret_cast<PrimitiveType*>(&rValue[0]);
                } else {
                    // not returning nullptr so, the return value can be subjected to
                    // arithmetic operations
                    return 0;
                }
            }
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @ref rContainer recursively to the
     * provided @ref pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Contaienr to copy data to the contiguous array.
     */
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        if (rContainer.size() != 0) {
            const auto stride = ValueTraits::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                ValueTraits::CopyToContiguousData(pContiguousDataBegin + i * stride, rContainer[i]);
            }
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @ref pContiguousDataBegin pointer
     * to the @ref rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     */
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        if (rContainer.size() != 0) {
            const auto stride = ValueTraits::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                ValueTraits::CopyFromContiguousData(rContainer[i], pContiguousDataBegin + i * stride);
            }
        }
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TCheckIndex, unsigned int TCurrentIndex>
    static constexpr bool IsDimensionDynamicImpl()
    {
        if constexpr(TCheckIndex == TCurrentIndex) {
            return true;
        } else {
            return ValueTraits::template IsDimensionDynamicImpl<TCheckIndex, TCurrentIndex + 1>();
        }
    }

    ///@}
    ///@name Friends
    ///@{

    template<class T>
    friend class DataTypeTraits;

    ///@}
};

template<class TDataType> class DataTypeTraits<DenseMatrix<TDataType>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = DenseMatrix<TDataType>;

    using ValueType = TDataType;

    using ValueTraits = DataTypeTraits<ValueType>;

    using PrimitiveType = typename ValueTraits::PrimitiveType;

    static constexpr bool IsDynamic = true;

    // boost ublas makes the underlying data structure contiguous for
    // any ValueType which are not dynamic recursively.
    static constexpr bool IsContiguous = !ValueTraits::IsDynamic;

    static constexpr unsigned int Dimension = ValueTraits::Dimension + 2;

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Returns whther the given @ref TCheckIndex is dynamic.
     *
     * This method returns true if the @ref TCheckIndex of the dimension
     * corresponds to a dynamic type. If not, this returns false.
     *
     * @tparam TCheckIndex          User input dimension index.
     * @return true                 If the input dimension index corresponds to dynamic data type.
     * @return false                If the input dimension index corresponds to static data type.
     */
    template<unsigned int TCheckIndex>
    static constexpr bool IsDimensionDynamic()
    {
        return IsDimensionDynamicImpl<TCheckIndex, 0>();
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @ref PrimitiveType values contained in
     * the @ref rContainer recursively.
     *
     * @param rContainer        The container to calculate the size.
     * @return TIndexType       Number of primitive type values in rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline TIndexType Size(const ContainerType& rValue)
    {
        return (rValue.size1() == 0 || rValue.size2() == 0 ? 0 : rValue.size1() * rValue.size2() * ValueTraits::template Size<TIndexType>(rValue.data()[0]));
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @ref rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline std::vector<TIndexType> Shape(const ContainerType& rContainer)
    {
        std::vector<TIndexType> shape(Dimension);
        Shape(rContainer, shape.data(), shape.data() + Dimension);
        return shape;
    }

    /**
     * @brief Fills the given array with the shape values in each dimension.
     *
     * @throws If the array is not of the required size.
     *
     * @tparam int              Type of the shape values.
     * @param pShapeBegin       Begin of the shape array.
     * @param pShapeEnd         End of the shape array.
     */
    template<class TIndexType = unsigned int>
    static inline void Shape(
        const ContainerType& rContainer,
        TIndexType* pShapeBegin,
        TIndexType* pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 2)
            << "Invalid dimensions given to fill for primitive data type [ Expected dimension >= 2, provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        if (rContainer.size1() > 0 && rContainer.size2() > 0) {
            ValueTraits::Shape(rContainer(0, 0), pShapeBegin + 2, pShapeEnd);
        } else {
            ValueTraits::Shape(ValueType{}, pShapeBegin + 2, pShapeEnd);
        }
        pShapeBegin[0] = rContainer.size1();
        pShapeBegin[1] = rContainer.size2();
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the @ref rShape
     * recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given
     * @ref rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<TIndexType>& rShape)
    {
        return Reshape(rContainer, rShape.data(), rShape.data() + rShape.size());
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the shape given by @ref pShapeBegin
     * and @ref pShapeEnd recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given shape
     * represented by @ref pShapeBegin and @ref pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 2)
            << "Invalid shape/dimension given for DenseMatrix data type [ Expected = " << Shape(rContainer) << ", provided = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        bool is_reshaped = false;

        if (rContainer.size1() != pShapeBegin[0] || rContainer.size2() != pShapeBegin[1]) {
            rContainer.resize(pShapeBegin[0], pShapeBegin[1], false);
            is_reshaped = true;
        }

        if constexpr(ValueTraits::IsDynamic) {
            std::for_each(rContainer.data().begin(), rContainer.data().end(), [&is_reshaped, pShapeBegin, pShapeEnd](auto& rValue) {
                is_reshaped = ValueTraits::Reshape(rValue, pShapeBegin + 2, pShapeEnd) || is_reshaped;
            });
        }

        return is_reshaped;
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            if constexpr(std::is_same_v<PrimitiveType, ValueType>) {
                return rValue.data().begin();
            } else {
                // since the underlying data structure for recusive static data types
                // is contiguous in ublas types, we can do the following to get the
                // contiguous array.
                if (rValue.size1() > 0 && rValue.size2() > 0) {
                    return reinterpret_cast<PrimitiveType const*>(&rValue(0, 0));
                } else {
                    // not returning nullptr so, the return value can be subjected to
                    // arithmetic operations
                    return 0;
                }
            }
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            if constexpr(std::is_same_v<PrimitiveType, ValueType>) {
                return rValue.data().begin();
            } else {
                // since the underlying data structure for recusive static data types
                // is contiguous in ublas types, we can do the following to get the
                // contiguous array.
                if (rValue.size1() > 0 && rValue.size2() > 0) {
                    return reinterpret_cast<PrimitiveType*>(&rValue(0, 0));
                } else {
                    // not returning nullptr so, the return value can be subjected to
                    // arithmetic operations
                    return 0;
                }
            }
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @ref rContainer recursively to the
     * provided @ref pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Contaienr to copy data to the contiguous array.
     */
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        if (rContainer.size1() != 0 && rContainer.size2() != 0) {
            const auto stride = ValueTraits::Size(rContainer(0, 0));
            for (unsigned int i = 0; i < rContainer.size1() * rContainer.size2(); ++i) {
                ValueTraits::CopyToContiguousData(pContiguousDataBegin + i * stride, rContainer.data()[i]);
            }
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @ref pContiguousDataBegin pointer
     * to the @ref rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     */
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        if (rContainer.size1() != 0 && rContainer.size2() != 0) {
            const auto stride = ValueTraits::Size(rContainer(0, 0));
            for (unsigned int i = 0; i < rContainer.size1() * rContainer.size2(); ++i) {
                ValueTraits::CopyFromContiguousData(rContainer.data()[i], pContiguousDataBegin + i * stride);
            }
        }
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TCheckIndex, unsigned int TCurrentIndex>
    static constexpr bool IsDimensionDynamicImpl()
    {
        if constexpr(TCheckIndex == TCurrentIndex || TCheckIndex == TCurrentIndex + 1) {
            return true;
        } else {
            return ValueTraits::template IsDimensionDynamicImpl<TCheckIndex, TCurrentIndex + 2>();
        }
    }

    ///@}
    ///@name Friends
    ///@{

    template<class T>
    friend class DataTypeTraits;

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

    static constexpr bool IsContiguous = true;

    static constexpr bool IsDynamic = true;

    static constexpr unsigned int Dimension = 1;

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Returns whther the given @ref TCheckIndex is dynamic.
     *
     * This method returns true if the @ref TCheckIndex of the dimension
     * corresponds to a dynamic type. If not, this returns false.
     *
     * @tparam TCheckIndex          User input dimension index.
     * @return true                 If the input dimension index corresponds to dynamic data type.
     * @return false                If the input dimension index corresponds to static data type.
     */
    template<unsigned int TCheckIndex>
    static constexpr bool IsDimensionDynamic()
    {
        return IsDimensionDynamicImpl<TCheckIndex, 0>();
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @ref PrimitiveType values contained in
     * the @ref rContainer recursively.
     *
     * @param rContainer        The container to calculate the size.
     * @return TIndexType       Number of primitive type values in rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline TIndexType Size(const ContainerType& rValue)
    {
        return rValue.size();
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @ref rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline std::vector<TIndexType> Shape(const ContainerType& rContainer)
    {
        std::vector<TIndexType> shape(Dimension);
        Shape(rContainer, shape.data(), shape.data() + Dimension);
        return shape;
    }

    /**
     * @brief Fills the given array with the shape values in each dimension.
     *
     * @throws If the array is not of the required size.
     *
     * @tparam int              Type of the shape values.
     * @param pShapeBegin       Begin of the shape array.
     * @param pShapeEnd         End of the shape array.
     */
    template<class TIndexType = unsigned int>
    static inline void Shape(
        const ContainerType& rContainer,
        TIndexType* pShapeBegin,
        TIndexType* pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) == 1)
            << "Invalid dimensions given to fill for std::string data type [ Expected dimension == 1, provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";
        pShapeBegin[0] = rContainer.size();
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the @ref rShape
     * recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given
     * @ref rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<TIndexType>& rShape)
    {
        return Reshape(rContainer, rShape.data(), rShape.data() + rShape.size());
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the shape given by @ref pShapeBegin
     * and @ref pShapeEnd recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given shape
     * represented by @ref pShapeBegin and @ref pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) == 1)
            << "Invalid shape/dimension given for std::string data type [ Expected = "
            << Shape(rContainer) << ", provided = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        bool is_reshaped = false;

        if (rContainer.size() != pShapeBegin[0]) {
            rContainer.resize(pShapeBegin[0], false);
            is_reshaped = true;
        }

        return is_reshaped;
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        return rValue.data();
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        return rValue.data();
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @ref rContainer recursively to the
     * provided @ref pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Contaienr to copy data to the contiguous array.
     */
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        std::copy(rContainer.begin(), rContainer.end(), pContiguousDataBegin);
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @ref pContiguousDataBegin pointer
     * to the @ref rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     */
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        std::copy(pContiguousDataBegin, pContiguousDataBegin + rContainer.size(), rContainer.begin());
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TCheckIndex, unsigned int TCurrentIndex>
    static constexpr bool IsDimensionDynamicImpl()
    {
        if constexpr(TCheckIndex == TCurrentIndex) {
            return true;
        } else {
            static_assert(TCheckIndex != TCheckIndex, "Invalid dimension index.");
        }
    }

    ///@}
    ///@name Friends
    ///@{

    template<class T>
    friend class DataTypeTraits;

    ///@}
};

template<class TDataType> class DataTypeTraits<std::vector<TDataType>>
{
public:
    ///@name Type definitions
    ///@{

    using ContainerType = std::vector<TDataType>;

    using ValueType = TDataType;

    using ValueTraits = DataTypeTraits<ValueType>;

    using PrimitiveType = typename ValueTraits::PrimitiveType;

    static constexpr bool IsContiguous = std::is_same_v<PrimitiveType, ValueType>;

    static constexpr bool IsDynamic = true;

    static constexpr unsigned int Dimension = ValueTraits::Dimension + 1;

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Returns whther the given @ref TCheckIndex is dynamic.
     *
     * This method returns true if the @ref TCheckIndex of the dimension
     * corresponds to a dynamic type. If not, this returns false.
     *
     * @tparam TCheckIndex          User input dimension index.
     * @return true                 If the input dimension index corresponds to dynamic data type.
     * @return false                If the input dimension index corresponds to static data type.
     */
    template<unsigned int TCheckIndex>
    static constexpr bool IsDimensionDynamic()
    {
        return IsDimensionDynamicImpl<TCheckIndex, 0>();
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @ref PrimitiveType values contained in
     * the @ref rContainer recursively.
     *
     * @param rContainer        The container to calculate the size.
     * @return TIndexType       Number of primitive type values in rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline TIndexType Size(const ContainerType& rValue)
    {
        return (rValue.empty() ? 0 : rValue.size() * ValueTraits::template Size<TIndexType>(rValue[0]));
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @ref rContainer.
     */
    template<class TIndexType = unsigned int>
    static inline std::vector<TIndexType> Shape(const ContainerType& rContainer)
    {
        std::vector<TIndexType> shape(Dimension);
        Shape(rContainer, shape.data(), shape.data() + Dimension);
        return shape;
    }

    /**
     * @brief Fills the given array with the shape values in each dimension.
     *
     * @throws If the array is not of the required size.
     *
     * @tparam int              Type of the shape values.
     * @param pShapeBegin       Begin of the shape array.
     * @param pShapeEnd         End of the shape array.
     */
    template<class TIndexType = unsigned int>
    static inline void Shape(
        const ContainerType& rContainer,
        TIndexType* pShapeBegin,
        TIndexType* pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 1)
            << "Invalid dimensions given to fill for primitive data type [ Expected dimension >= 1, provided shape = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        if (rContainer.empty()) {
            ValueTraits::Shape(ValueType{}, pShapeBegin + 1, pShapeEnd);
        } else {
            ValueTraits::Shape(rContainer[0], pShapeBegin + 1, pShapeEnd);
        }
        pShapeBegin[0] = rContainer.size();
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the @ref rShape
     * recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given
     * @ref rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        const std::vector<TIndexType>& rShape)
    {
        return Reshape(rContainer, rShape.data(), rShape.data() + rShape.size());
    }

    /**
     * @brief Reshapes the given value to given shape.
     *
     * This method reshapes the given @ref rContainer to the shape given by @ref pShapeBegin
     * and @ref pShapeEnd recursively.
     *
     * If this method has changed size of @ref rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @ref rContainer to comply with the given shape
     * represented by @ref pShapeBegin and @ref pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContaienr has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 1)
            << "Invalid shape/dimension given for std::vector data type [ Expected = " << Shape(rContainer) << ", provided = "
            << std::vector<TIndexType>(pShapeBegin, pShapeEnd) << " ].\n";

        bool is_reshaped = false;

        if (rContainer.size() != pShapeBegin[0]) {
            rContainer.resize(pShapeBegin[0]);
            is_reshaped = true;
        }

        if constexpr(ValueTraits::IsDynamic) {
            std::for_each(rContainer.begin(), rContainer.end(), [&is_reshaped, pShapeBegin, pShapeEnd](auto& rValue) {
                is_reshaped = ValueTraits::Reshape(rValue, pShapeBegin + 1, pShapeEnd) || is_reshaped;
            });
        }

        return is_reshaped;
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType const * GetContiguousData(const ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            return rValue.data();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @ref rContainer.
     * If the underlying data structure of @ref rContainer is not contiguous, this will
     * throw a compiler time error.
     *
     * @param rValue                    Container to retireve the contiguous array pointer.
     * @return PrimitiveType const*     Contiguous array pointer.
     */
    inline static PrimitiveType * GetContiguousData(ContainerType& rValue)
    {
        if constexpr(IsContiguous) {
            return rValue.data();
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @ref rContainer recursively to the
     * provided @ref pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Contaienr to copy data to the contiguous array.
     */
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        if (rContainer.size() != 0) {
            const auto stride = ValueTraits::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                ValueTraits::CopyToContiguousData(pContiguousDataBegin + i * stride, rContainer[i]);
            }
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @ref pContiguousDataBegin pointer
     * to the @ref rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @ref pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     */
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin)
    {
        if (rContainer.size() != 0) {
            const auto stride = ValueTraits::Size(rContainer[0]);
            for (unsigned int i = 0; i < rContainer.size(); ++i) {
                ValueTraits::CopyFromContiguousData(rContainer[i], pContiguousDataBegin + i * stride);
            }
        }
    }

    ///@}

private:
    ///@name Private static operations
    ///@{

    template<unsigned int TCheckIndex, unsigned int TCurrentIndex>
    static constexpr bool IsDimensionDynamicImpl()
    {
        if constexpr(TCheckIndex == TCurrentIndex) {
            return true;
        } else {
            return ValueTraits::template IsDimensionDynamicImpl<TCheckIndex, TCurrentIndex + 1>();
        }
    }

    ///@}
    ///@name Friends
    ///@{

    template<class T>
    friend class DataTypeTraits;

    ///@}
};

}; // namespace Kratos
