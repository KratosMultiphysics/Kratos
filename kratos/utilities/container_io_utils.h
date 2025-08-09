//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"
#include "utilities/string_utilities.h"

// Application includes

namespace Kratos
{
///@name Kratos classes
///@{

class KRATOS_API(KRATOS_CORE) ContainerIOUtils
{
public:
    ///@name Public static operations
    ///@{

    /**
     * @brief Copies the values of each entity given by getter method in the container to contiguous array.
     * @details This method copies the values of each entity given by the getter method in the container to
     *          the given span.
     *              The @p TGetterType can be a lambda function or a function having the following signature
     *                  [](TDataType& rValue, const EntityType& rEntity) -> void.
     *                  The @p rGetter should get the value from @p rEntity and store it in the @p rValue input
     *                  variable. The full size value should be read, even if you want to have few components
     *                  of the @p rValue in the final contiguous array.
     *
     *              The @p pShapeBegin and @p pShapeEnd represents the required shape of the contiguous array
     *              represented by @p rDataSpan. Therefore, first dimension of the shape should represent the
     *              number of entities in the @p rContainer. The rest of the values of shape represent number of
     *              components in respective dimensions in the @p rValue of which is taken from the @p rGetter.
     *              This may be a subset of the @p TDataType.
     *                  Ex: if the @p TDataType is @ref array_1d<double_3> then, The vector representing
     *                      [ @p pShapeBegin + 1, @p pShapeEnd ) may be [1], [2] or [3].
     *                      if the @p TDataType is a dynamic type such as @ref DenseVector or @ref DenseMatrix,
     *                      then, the [ @p pShapeBegin + 1, @p pShapeEnd ) should represent values in each dimension
     *                      less than or equal to the values in each dimension of @p rValue.
     *
     *              The @p rDataSpan should represent a data span which can hold all the entities' data
     *              which are reshaped to the [ @p pShapeBegin, @p pShapeEnd ) shape.
     *
     * @throws If the first dimension of the shape represented by [ @p pShapeBegin, @p pShapeEnd ) is
     *         not equal to the number of entities in the @p rContainer.
     * @throws If the shape represented by [ @p pShapeBegin + 1, @p pShapeEnd ) is not a valid shape
     *         representing @p TDataType.
     * @throws If the size of the given @p rDataSpan does not hold enough memory allocated to
     *         put all the values.
     *
     * @tparam TDataType            The type of the value in each entity which is taken from the @p rGetter.
     * @tparam TContainerType       The type of the container containing all the entities.
     * @tparam TSpanType            The type of the span, to which all the retrieved and resized data from each entity will be copied to.
     * @tparam TIntegerType         The type of the integer used to represent the shape.
     * @tparam TGetterType          The type of the getter function.
     * @param rContainer            The container containing all the entities.
     * @param rDataSpan             The span, to which all the retrieved and resized data from each entity will be copied to.
     * @param pShapeBegin           Begining of the shape.
     * @param pShapeEnd             End of the shape.
     * @param rGetter               The getter function which retrieves values of @p TDataType from each entity in @p rContainer.
     */
    template<class TDataType, class TContainerType, class TSpanType, class TIntegerType, class TGetterType>
    static void CopyToContiguousArray(
        const TContainerType& rContainer,
        const TSpanType& rDataSpan,
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd,
        const TGetterType& rGetter)
    {
        KRATOS_TRY

        using value_type_traits = DataTypeTraits<TDataType>;

        KRATOS_ERROR_IF_NOT(pShapeBegin[0] == rContainer.size())
            << "First dimension of the  shape mismatch with the container size [ "
            << " container size = " << rContainer.size() << ", first dimension of the shape = " << pShapeBegin[0] << " ].\n";

        KRATOS_ERROR_IF_NOT(DataTypeTraits<TDataType>::IsValidShape(pShapeBegin + 1, pShapeEnd))
            << "Invalid data shape provided. [ data shape provided = ["
            << StringUtilities::JoinValues(pShapeBegin + 1, pShapeEnd, ",")
            << "], max possible sizes in each dimension  = "
            << DataTypeTraits<TDataType>::Shape(TDataType{}) << " ].\n";

        KRATOS_ERROR_IF_NOT(rDataSpan.size() == pShapeBegin[0] * DataTypeTraits<TDataType>::Size(pShapeBegin + 1, pShapeEnd))
            << "The span size mismatch with required size from the given shape [ "
            << "span size = " << rDataSpan.size() << ", shape = [" << StringUtilities::JoinValues(pShapeBegin, pShapeEnd, ",") << "] ].\n";

        const auto stride = value_type_traits::Size(pShapeBegin + 1, pShapeEnd);

        IndexPartition<unsigned int>(rContainer.size()).for_each(TDataType{}, [&rContainer,  &rGetter, pShapeBegin, pShapeEnd, stride, rDataSpan](const auto Index, auto& rTLS) {
            const auto& r_entity = *(rContainer.begin() + Index);
            rGetter(rTLS, r_entity);
            auto p_subrange_begin = rDataSpan.data() + Index * stride;
            value_type_traits::CopyToContiguousData(p_subrange_begin, rTLS, pShapeBegin + 1, pShapeEnd);
        });

        KRATOS_CATCH("");
    }

    /**
     * @brief Copy data from contiguous array to values given by the reference getter of each entity.
     * @details This method copies data from the contiguous array to values given by the reference getter
     *          of each entity.
     *              The @p rReferenceGetter lambda function should have the following signature.
     *                  []( EntityType& rEntity) -> TDataType&
     *              The reference getter should return a reference to the data which needs its' values
     *              copied from the contiguous array. The size checks for the dynamic types will be done
     *              in this method and errors will be thrown.
     *
     *              The @p pShapeBegin and @p pShapeEnd represents the shape of the contiguous array
     *              represented by @p rDataSpan. Therefore, first dimension of the shape should represent the
     *              number of entities in the @p rContainer. The rest of the values of shape represent number of
     *              components in respective dimensions in the @p rValue of which is set from the @p rSetter.
     *              This may be a subset of the @p TDataType.
     *                  Ex: if the @p TDataType is @ref array_1d<double_3> then, The vector representing
     *                      [ @p pShapeBegin + 1, @p pShapeEnd ) may be [1], [2] or [3].
     *                      if the @p TDataType is a dynamic type such as @ref DenseVector or @ref DenseMatrix,
     *                      then, the [ @p pShapeBegin + 1, @p pShapeEnd ) should represent values in each dimension
     *                      less than or equal to the values in each dimension of @p rValue.
     *
     *              The @p rDataSpan should represent a data span which can hold all the entities' data
     *
     *              which are reshaped to the [ @p pShapeBegin, @p pShapeEnd ) shape.
     * @throws If the first dimension of the shape represented by [ @p pShapeBegin, @p pShapeEnd ) is
     *         not equal to the number of entities in the @p rContainer.
     * @throws If the shape represented by [ @p pShapeBegin + 1, @p pShapeEnd ) is not a valid shape
     *         representing @p TDataType.
     * @throws If the size of the given @p rDataSpan does not hold enough memory allocated to
     *         put all the values.
     * @throws If the @p TDataType values is dynamic and they are not correctly sized in any of the
     *         entities.
     *
     * @tparam TDataType            The type of the value in each entity which is written by the @p rSetter.
     * @tparam TContainerType       The type of the container containing all the entities.
     * @tparam TSpanType            The type of the span, to which all the stored and resized data read to each entity.
     * @tparam TIntegerType         The type of the integer used to represent the shape.
     * @tparam TReferenceGetterType The type of the reference getter function.
     * @param rContainer            The container containing all the entities.
     * @param rDataSpan             The span, from which all the read and resized data to each entity will be written to.
     * @param pShapeBegin           Begining of the shape.
     * @param pShapeEnd             End of the shape.
     * @param rReferenceGetter      The reference getter function which gets the reference values of @p TDataType to each entity in @p rContainer.
     */
    template<class TDataType, class TContainerType, class TSpanType, class TIntegerType, class TReferenceGetterType>
    static void CopyFromContiguousDataArray(
        TContainerType& rContainer,
        const TSpanType& rDataSpan,
        TIntegerType const * pShapeBegin,
        TIntegerType const * pShapeEnd,
        const TReferenceGetterType& rReferenceGetter)
    {
        KRATOS_TRY

        using value_type_traits = DataTypeTraits<TDataType>;

        KRATOS_ERROR_IF_NOT(pShapeBegin[0] == rContainer.size())
            << "First dimension of the  shape mismatch with the container size [ "
            << " container size = " << rContainer.size() << ", first dimension of the shape = " << pShapeBegin[0] << " ].\n";

        KRATOS_ERROR_IF_NOT(DataTypeTraits<TDataType>::IsValidShape(pShapeBegin + 1, pShapeEnd))
            << "Invalid data shape provided. [ data shape provided = ["
            << StringUtilities::JoinValues(pShapeBegin + 1, pShapeEnd, ",")
            << "], max possible sizes in each dimension  = "
            << DataTypeTraits<TDataType>::Shape(TDataType{}) << " ].\n";

        KRATOS_ERROR_IF_NOT(rDataSpan.size() == pShapeBegin[0] * DataTypeTraits<TDataType>::Size(pShapeBegin + 1, pShapeEnd))
            << "The span size mismatch with required size from the given shape [ "
            << "span size = " << rDataSpan.size() << ", shape = [" << StringUtilities::JoinValues(pShapeBegin, pShapeEnd, ",") << "] ].\n";

        const auto stride = value_type_traits::Size(pShapeBegin + 1, pShapeEnd);

        IndexPartition<unsigned int>(rContainer.size()).for_each([&rContainer, &rReferenceGetter, rDataSpan, pShapeBegin, pShapeEnd, stride](const auto Index) {
            auto p_subrange_begin = rDataSpan.data() + Index * stride;
            value_type_traits::CopyFromContiguousData(rReferenceGetter(*(rContainer.begin() + Index)), p_subrange_begin, pShapeBegin + 1, pShapeEnd);
        });

        KRATOS_CATCH("");
    }

    ///@}
};

///@}
} // namespace Kratos