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
#include <variant>
#include <type_traits>
#include <vector>

// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos {

template<class T, class... TList>
constexpr bool IsInList = std::disjunction<std::is_same<T, TList>...>::value;

template<class... TAlternativesList>
struct HoldsAlternative
{
    template<class... TVariantTypes>
    static bool Evaluate(std::variant<TVariantTypes...>&& rVariant) {
        return  (... || std::holds_alternative<TAlternativesList>(rVariant));
    }
};

template <class T>
struct BareTypeImpl {
    using type = std::remove_cv_t<std::remove_pointer_t<std::remove_cv_t<std::remove_reference_t<T>>>>;
};

template <typename T>
using BareType = typename BareTypeImpl<T>::type;

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
     * @brief Returns whther the given @p TCheckIndex is dynamic.
     *
     * This method returns true if the @p TCheckIndex of the dimension
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
     * @brief Checks if the given shape is valid.
     * @details This method checks if the given shape is valid in the sense
     *          that, dimensionality is correct as well as the number of components for each
     *          static dimension is less or equal to the actual static data types'
     *          number of components.
     * @param pShapeBegin           Beginning of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIntegerType>
    static bool inline IsValidShape(
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd)
    {
        // this class is reserved for primitive types such as int, double, bool,...
        // therefore, there should not be any dimensionality. Should only be an empty shape.
        return pShapeBegin == pShapeEnd;
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
     * @brief Get the size of the underlying container given the shape
     *
     * This method returns number of @p PrimitiveType values contained recursively
     * using the given shape with @p pShapeBegin indicating the pointer to the
     * start of the shape array and @p pShapeEnd indicating the pointer to the
     * end of the shape array.
     *
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     * @return TIndexType           Number of primitive type values in the shape.
     */
    template<class TIteratorType, class TIndexType = unsigned int>
    static inline TIndexType Size(
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        return 1;
    }

    /**
     * @brief Returns a vector with the shape of the value.
     * @details This return the shape of the value. Scalars have
     *          shape of [] (an empty vector).
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
     * This method copies content of the @p rValue to a contiguous array
     * at @p pContiguousDataBegin. The array pointed by @p pContiguousDataBegin
     * should be sized correctly.
     *
     * @warning This may seg fault if the @p pContiguousDataBegin is not correctly sized.
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
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may segfault if the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
     * @param pShapeBegin               Begining of the shape array.
     * @param pShapeEnd                 End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        *pContiguousDataBegin = rContainer;
    }

    /**
     * @brief Copies data from contiguous array to rValue.
     *
     * This method copies data from contiguous array to the passed @p rValue.
     * The contiguous array is given by the @refpContiguousDataBegin.
     *
     * @warning This may seg fault if the @p pContiguousDataBegin is not correctly sized.
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

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        rContainer = *pContiguousDataBegin;
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
            static_assert(sizeof(TCheckIndex) == 0 && false, "Invalid dimension index.");
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
     * @brief Returns whther the given @p TCheckIndex is dynamic.
     *
     * This method returns true if the @p TCheckIndex of the dimension
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
     * @brief Checks if the given shape is valid.
     * @details This method checks if the given shape is valid in the sense
     *          that, dimensionality is correct as well as the number of components for each
     *          static dimension is less or equal to the actual static data types'
     *          number of components.
     * @param pShapeBegin           Beginning of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIntegerType>
    static bool inline IsValidShape(
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd)
    {
        // this is an array_1d, so first dimension is used by the array_1d.
        // the static size check is done here.
        return pShapeBegin != pShapeEnd && (*pShapeBegin) <= TSize && ValueTraits::template IsValidShape<TIntegerType>(pShapeBegin + 1, pShapeEnd);
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @p PrimitiveType values contained in
     * the @p rContainer recursively.
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
     * @brief Get the size of the underlying container given the shape
     *
     * This method returns number of @p PrimitiveType values contained recursively
     * using the given shape with @p pShapeBegin indicating the pointer to the
     * start of the shape array and @p pShapeEnd indicating the pointer to the
     * end of the shape array.
     *
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     * @return TIndexType           Number of primitive type values in the shape.
     */
    template<class TIteratorType, class TIndexType = unsigned int>
    static inline TIndexType Size(
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        return (*pShapeBegin) * ValueTraits::template Size<TIteratorType, TIndexType>(pShapeBegin + 1, pShapeEnd);
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @p rContainer.
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
     * This method reshapes the given @p rContainer to the @p rShape
     * recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given
     * @p rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method reshapes the given @p rContainer to the shape given by @p pShapeBegin
     * and @p pShapeEnd recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given shape
     * represented by @p pShapeBegin and @p pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
     */
    template<class TIndexType = unsigned int>
    static inline bool Reshape(
        ContainerType& rContainer,
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        KRATOS_ERROR_IF_NOT(std::distance(pShapeBegin, pShapeEnd) >= 1 && *pShapeBegin <= static_cast<TIndexType>(TSize))
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
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "This should be only called if the rValue is contiguous only.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "This should be only called if the rValue is contiguous only.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
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
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
     * @param pShapeBegin               Begining of the shape array.
     * @param pShapeEnd                 End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        const auto stride = ValueTraits::Size(pShapeBegin + 1, pShapeEnd);
        for (unsigned int i = 0; i < *pShapeBegin; ++i) {
            ValueTraits::template CopyToContiguousData<TIteratorType>(pContiguousDataBegin + i * stride, rContainer[i], pShapeBegin + 1, pShapeEnd);
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
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

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        const auto stride = ValueTraits::Size(pShapeBegin + 1, pShapeEnd);
        for (unsigned int i = 0; i < *pShapeBegin; ++i) {
            ValueTraits::template CopyFromContiguousData<TIteratorType>(rContainer[i], pContiguousDataBegin + i * stride, pShapeBegin + 1, pShapeEnd);
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
     * @brief Returns whther the given @p TCheckIndex is dynamic.
     *
     * This method returns true if the @p TCheckIndex of the dimension
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
     * @brief Checks if the given shape is valid.
     * @details This method checks if the given shape is valid in the sense
     *          that, dimensionality is correct as well as the number of components for each
     *          static dimension is less or equal to the actual static data types'
     *          number of components.
     * @param pShapeBegin           Beginning of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIntegerType>
    static bool inline IsValidShape(
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd)
    {
        // this is a dense vector, hence the first dimension of the given shape
        // is used to identify the number of components in the dense vector.
        // Since the dense vector is dynamic, this does not check the given
        // number of components in the shape with what will be available.
        // it only checks whether the dimension is available in the shape.
        return pShapeBegin != pShapeEnd && ValueTraits::template IsValidShape<TIntegerType>(pShapeBegin + 1, pShapeEnd);
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @p PrimitiveType values contained in
     * the @p rContainer recursively.
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
     * @brief Get the size of the underlying container given the shape
     *
     * This method returns number of @p PrimitiveType values contained recursively
     * using the given shape with @p pShapeBegin indicating the pointer to the
     * start of the shape array and @p pShapeEnd indicating the pointer to the
     * end of the shape array.
     *
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     * @return TIndexType           Number of primitive type values in the shape.
     */
    template<class TIteratorType, class TIndexType = unsigned int>
    static inline TIndexType Size(
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        return (*pShapeBegin) * ValueTraits::template Size<TIteratorType, TIndexType>(pShapeBegin + 1, pShapeEnd);
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @p rContainer.
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
     * This method reshapes the given @p rContainer to the @p rShape
     * recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given
     * @p rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method reshapes the given @p rContainer to the shape given by @p pShapeBegin
     * and @p pShapeEnd recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given shape
     * represented by @p pShapeBegin and @p pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
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
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
     * @param pShapeBegin               Begining of the shape array.
     * @param pShapeEnd                 End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        // since the range [pShapeBegin, pShapeEnd) is checked using the IsValidShape method,
        // here we only check for the dynamic size.
        if (*pShapeBegin <= rContainer.size()) {
            const auto stride = ValueTraits::Size(pShapeBegin + 1, pShapeEnd);
            for (unsigned int i = 0; i < *pShapeBegin; ++i) {
                ValueTraits::template CopyToContiguousData<TIteratorType>(pContiguousDataBegin + i * stride, rContainer[i], pShapeBegin + 1, pShapeEnd);
            }
        } else {
            KRATOS_ERROR
                << "The given number of components are larger than the data size of DenseVector [ number of components in the dimension = "
                << *pShapeBegin << ", number of components available in the data = " << rContainer.size() << " ].\n";
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
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

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        if (rContainer.size() >= *pShapeBegin) {
            const auto stride = ValueTraits::Size(pShapeBegin + 1, pShapeEnd);
            for (unsigned int i = 0; i < *pShapeBegin; ++i) {
                ValueTraits::template CopyFromContiguousData<TIteratorType>(rContainer[i], pContiguousDataBegin + i * stride, pShapeBegin + 1, pShapeEnd);
            }
        } else {
            KRATOS_ERROR
                << "The given number of components are larger than the data size of DenseVector [ number of components in the dimension = "
                << *pShapeBegin << ", number of components available in the data = " << rContainer.size() << " ].\n";
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
     * @brief Returns whther the given @p TCheckIndex is dynamic.
     *
     * This method returns true if the @p TCheckIndex of the dimension
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
     * @brief Checks if the given shape is valid.
     * @details This method checks if the given shape is valid in the sense
     *          that, dimensionality is correct as well as the number of components for each
     *          static dimension is less or equal to the actual static data types'
     *          number of components.
     * @param pShapeBegin           Beginning of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIntegerType>
    static bool inline IsValidShape(
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd)
    {
        // this is for the matrix type.
        // first two dimensions of the shape is used to define the matrix row and column
        // number of components. Since matrix is dynamic, this
        // only check whether the required dimensions are available.
        return pShapeBegin != pShapeEnd && (pShapeBegin + 1) != pShapeEnd && ValueTraits::template IsValidShape<TIntegerType>(pShapeBegin + 2, pShapeEnd);
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @p PrimitiveType values contained in
     * the @p rContainer recursively.
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
     * @brief Get the size of the underlying container given the shape
     *
     * This method returns number of @p PrimitiveType values contained recursively
     * using the given shape with @p pShapeBegin indicating the pointer to the
     * start of the shape array and @p pShapeEnd indicating the pointer to the
     * end of the shape array.
     *
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     * @return TIndexType           Number of primitive type values in the shape.
     */
    template<class TIteratorType, class TIndexType = unsigned int>
    static inline TIndexType Size(
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        return (*pShapeBegin) * (*(pShapeBegin + 1)) * ValueTraits::template Size<TIteratorType, TIndexType>(pShapeBegin + 2, pShapeEnd);
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @p rContainer.
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
     * This method reshapes the given @p rContainer to the @p rShape
     * recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given
     * @p rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method reshapes the given @p rContainer to the shape given by @p pShapeBegin
     * and @p pShapeEnd recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given shape
     * represented by @p pShapeBegin and @p pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
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
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
     * @param pShapeBegin               Begining of the shape array.
     * @param pShapeEnd                 End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        // since the range [pShapeBegin, pShapeEnd) is checked using the IsValidShape method,
        // here we only check for the dynamic size.
        if (*pShapeBegin <= rContainer.size1() && *(pShapeBegin + 1) <= rContainer.size2()) {
            const auto stride = ValueTraits::Size(pShapeBegin + 2, pShapeEnd);
            for (unsigned int i = 0; i < *pShapeBegin; ++i) {
                for (unsigned int j = 0; j < *(pShapeBegin + 1); ++j) {
                    ValueTraits::template CopyToContiguousData<TIteratorType>(pContiguousDataBegin + i * stride * (*(pShapeBegin + 1)) + j * stride, rContainer(i, j), pShapeBegin + 2, pShapeEnd);
                }
            }
        } else {
            KRATOS_ERROR
                << "The given number of components are larger than the data size of DenseMatrix [ number of components in the dimensions = ("
                << *pShapeBegin << ", " << *(pShapeBegin + 1) << "), number of components available in the data = ("
                << rContainer.size1() << ", " << rContainer.size2() << ") ].\n";
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
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

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        // since the range [pShapeBegin, pShapeEnd) is checked using the IsValidShape method,
        // here we only check for the dynamic size.
        if (*pShapeBegin <= rContainer.size1() && *(pShapeBegin + 1) <= rContainer.size2()) {
            const auto stride = ValueTraits::Size(pShapeBegin + 2, pShapeEnd);
            for (unsigned int i = 0; i < *pShapeBegin; ++i) {
                for (unsigned int j = 0; j < *(pShapeBegin + 1); ++j) {
                    ValueTraits::template CopyFromContiguousData<TIteratorType>(rContainer(i, j), pContiguousDataBegin + i * stride * (*(pShapeBegin + 1)) + j * stride, pShapeBegin + 2, pShapeEnd);
                }
            }
        } else {
            KRATOS_ERROR
                << "The given number of components are larger than the data size of DenseMatrix [ number of components in the dimensions = ("
                << *pShapeBegin << ", " << *(pShapeBegin + 1) << "), number of components available in the data = ("
                << rContainer.size1() << ", " << rContainer.size2() << " ].\n";
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
     * @brief Returns whther the given @p TCheckIndex is dynamic.
     *
     * This method returns true if the @p TCheckIndex of the dimension
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
     * @brief Checks if the given shape is valid.
     * @details This method checks if the given shape is valid in the sense
     *          that, dimensionality is correct as well as the number of components for each
     *          static dimension is less or equal to the actual static data types'
     *          number of components.
     * @param pShapeBegin           Beginning of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIntegerType>
    static bool inline IsValidShape(
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd)
    {
        // this is also considered as an array of chars.
        // hence it should only have one dimension.
        // since string is dynamic, number of components
        // in the dimension is not checked. Only the
        // availability of the dimension is checked.
        return (pShapeBegin + 1) == pShapeEnd;
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @p PrimitiveType values contained in
     * the @p rContainer recursively.
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
     * @brief Get the size of the underlying container given the shape
     *
     * This method returns number of @p PrimitiveType values contained recursively
     * using the given shape with @p pShapeBegin indicating the pointer to the
     * start of the shape array and @p pShapeEnd indicating the pointer to the
     * end of the shape array.
     *
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     * @return TIndexType           Number of primitive type values in the shape.
     */
    template<class TIteratorType, class TIndexType = unsigned int>
    static inline TIndexType Size(
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        return *pShapeBegin;
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @p rContainer.
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
     * This method reshapes the given @p rContainer to the @p rShape
     * recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given
     * @p rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method reshapes the given @p rContainer to the shape given by @p pShapeBegin
     * and @p pShapeEnd recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given shape
     * represented by @p pShapeBegin and @p pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
     */
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer)
    {
        std::copy(rContainer.begin(), rContainer.end(), pContiguousDataBegin);
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
     * @param pShapeBegin               Begining of the shape array.
     * @param pShapeEnd                 End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        // since the range [pShapeBegin, pShapeEnd) is checked using the IsValidShape method,
        // here we only check for the dynamic size.
        if (*pShapeBegin <= rContainer.size()) {
            std::copy(rContainer.begin(), rContainer.begin() + (*pShapeBegin), pContiguousDataBegin);
        } else {
            KRATOS_ERROR
                << "The given number of components are larger than the data size of std::string [ number of components in the dimension = "
                << *pShapeBegin << ", number of components available in the data = " << rContainer.size() << " ].\n";
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
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

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        std::copy(pContiguousDataBegin, pContiguousDataBegin + *pShapeBegin, rContainer.begin());
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
            static_assert(sizeof(TCheckIndex) == 0 && false, "Invalid dimension index.");
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
     * @brief Returns whther the given @p TCheckIndex is dynamic.
     *
     * This method returns true if the @p TCheckIndex of the dimension
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
     * @brief Checks if the given shape is valid.
     * @details This method checks if the given shape is valid in the sense
     *          that, dimensionality is correct as well as the number of components for each
     *          static dimension is less or equal to the actual static data types'
     *          number of components.
     * @param pShapeBegin           Beginning of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIntegerType>
    static bool inline IsValidShape(
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd)
    {
        // This is a vector type, and first dimension represents
        // the number of components in the vector. Hence the availability
        // of the dimension is checked. Since std::vector is dynamic,
        // the number of component in the std::vector is not checked.
        return pShapeBegin != pShapeEnd && ValueTraits::template IsValidShape<TIntegerType>(pShapeBegin + 1, pShapeEnd);
    }

    /**
     * @brief Gets the size of underlying rContainer.
     *
     * This method returns number of @p PrimitiveType values contained in
     * the @p rContainer recursively.
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
     * @brief Get the size of the underlying container given the shape
     *
     * This method returns number of @p PrimitiveType values contained recursively
     * using the given shape with @p pShapeBegin indicating the pointer to the
     * start of the shape array and @p pShapeEnd indicating the pointer to the
     * end of the shape array.
     *
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     * @return TIndexType           Number of primitive type values in the shape.
     */
    template<class TIteratorType, class TIndexType = unsigned int>
    static inline TIndexType Size(
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        return (*pShapeBegin) * ValueTraits::template Size<TIteratorType, TIndexType>(pShapeBegin + 1, pShapeEnd);
    }

    /**
     * @brief Get the shape of the rContainer.
     *
     * This method returns the shape of the rContainer. Shape is calculated
     * in a recursive manner.
     *
     * @param rContainer                    Value to compute the shape.
     * @return std::vector<TIndexType>    Shape of the @p rContainer.
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
     * This method reshapes the given @p rContainer to the @p rShape
     * recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given
     * @p rShape.
     *
     * @param rContainer    Container to be resized recursively.
     * @param rShape        Shape to be used in resizing.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method reshapes the given @p rContainer to the shape given by @p pShapeBegin
     * and @p pShapeEnd recursively.
     *
     * If this method has changed size of @p rContainer or any of the elements
     * of the container, then true is returned. False is returned if no
     * change is required in @p rContainer to comply with the given shape
     * represented by @p pShapeBegin and @p pShapeEnd.
     *
     * @param rContainer    Container to be resized recursively.
     * @param pShapeBegin   Begin of the shape vector.
     * @param pShapeEnd     End of the shape vector.
     * @return true         If the rContainer or its elements has changed due to resizing.
     * @return false        If the rContainer has not changed.
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
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Get the Contiguous data pointer of the given container.
     *
     * This method returns the underlying contiguous data ppinter of the given @p rContainer.
     * If the underlying data structure of @p rContainer is not contiguous, this will
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
            static_assert(sizeof(TDataType) == 0 && false, "GetContiguousData should only be called if rValue is contiguous.");
        }
    }

    /**
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
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
     * @brief Copies the given container element value to contiguous array.
     *
     * This method copies all the elements of @p rContainer recursively to the
     * provided @p pContiguousDataBegin.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param pContiguousDataBegin      Contiguous array pointer.
     * @param rContainer                Container to copy data to the contiguous array.
     * @param pShapeBegin               Begining of the shape array.
     * @param pShapeEnd                 End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyToContiguousData(
        PrimitiveType* pContiguousDataBegin,
        const ContainerType& rContainer,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        // since the range [pShapeBegin, pShapeEnd) is checked using the IsValidShape method,
        // here we only check for the dynamic size.
        if (*pShapeBegin <= rContainer.size()) {
            const auto stride = ValueTraits::Size(pShapeBegin + 1, pShapeEnd);
            for (unsigned int i = 0; i < *pShapeBegin; ++i) {
                ValueTraits::template CopyToContiguousData<TIteratorType>(pContiguousDataBegin + i * stride, rContainer[i], pShapeBegin + 1, pShapeEnd);
            }
        } else {
            KRATOS_ERROR
                << "The given number of components are larger than the data size of std::vector [ number of components in the dimension = "
                << *pShapeBegin << ", number of components available in the data = " << rContainer.size() << " ].\n";
        }
    }

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
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

    /**
     * @brief Copies contiguous values to the container.
     *
     * This method copies all the contiguous values in the given @p pContiguousDataBegin pointer
     * to the @p rContainer.
     *
     * @warning This may seg-fault if the the contiguous array given by @p pContiguousDataBegin
     *          is not correctly sized.
     *
     * @param rContainer            Container to copy values to.
     * @param pContiguousDataBegin  Contiguous data container from which values are copied from.
     * @param pShapeBegin           Begining of the shape array.
     * @param pShapeEnd             End of the shape array.
     */
    template<class TIteratorType>
    inline static void CopyFromContiguousData(
        ContainerType& rContainer,
        PrimitiveType const * pContiguousDataBegin,
        TIteratorType pShapeBegin,
        TIteratorType pShapeEnd)
    {
        if (rContainer.size() >= *pShapeBegin) {
            const auto stride = ValueTraits::Size(pShapeBegin + 1, pShapeEnd);
            for (unsigned int i = 0; i < *pShapeBegin; ++i) {
                ValueTraits::template CopyFromContiguousData<TIteratorType>(rContainer[i], pContiguousDataBegin + i * stride, pShapeBegin + 1, pShapeEnd);
            }
        } else {
            KRATOS_ERROR
                << "The given number of components are larger than the data size of std::vector [ number of components in the dimension = "
                << *pShapeBegin << ", number of components available in the data = " << rContainer.size() << " ].\n";
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