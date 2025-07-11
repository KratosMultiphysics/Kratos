//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <sstream>
#include <variant>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/data_type_traits.h"
#include "utilities/container_io_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) TensorAdaptorUtils
{
public:
    ///@name Type definitions
    ///@{

    using VariablePointerType = std::variant<
                                        Variable<double> const *,
                                        Variable<array_1d<double, 3>> const *,
                                        Variable<array_1d<double, 4>> const *,
                                        Variable<array_1d<double, 6>> const *,
                                        Variable<array_1d<double, 9>> const *,
                                        Variable<Vector> const *,
                                        Variable<Matrix> const *
                                    >;

    ///@}
    ///@name Public static operations
    ///@{

    template<class TDataType, class TContainerType, class TGetterType>
    static DenseVector<unsigned int> GetTensorShape(
        const TContainerType& rContainer,
        const TGetterType& rGetter)
    {
        KRATOS_TRY

        using value_traits = DataTypeTraits<TDataType>;

        KRATOS_ERROR_IF(rContainer.empty())
            << "The given container does not have any "
            << ModelPart::Container<TContainerType>::GetEntityName() << "s.\n";

        TDataType dummy{};
        rGetter(dummy, rContainer.front());

        DenseVector<unsigned int> tensor_shape(value_traits::Dimension + 1);
        tensor_shape[0] = rContainer.size();
        DataTypeTraits<TDataType>::Shape(dummy, tensor_shape.data().begin() + 1, tensor_shape.data().begin() + tensor_shape.size());

        return tensor_shape;

        KRATOS_CATCH("");
    }

    template<class TDataType, class TContainerType, class TIntegerType>
    static DenseVector<unsigned int> GetTensorShape(
        const TContainerType& rContainer,
        TIntegerType const * rDataShapeBegin,
        TIntegerType const * rDataShapeEnd)
    {
        KRATOS_TRY

        using value_traits = DataTypeTraits<TDataType>;

        KRATOS_ERROR_IF(rContainer.empty())
            << "The given container does not have any "
            << ModelPart::Container<TContainerType>::GetEntityName() << "s.\n";

        KRATOS_ERROR_IF_NOT(value_traits::IsValidShape(rDataShapeBegin, rDataShapeEnd))
            << "Invalid data shape provided. [ data shape provided = [" << StringUtilities::JoinValues(rDataShapeBegin, rDataShapeEnd, ",")
            << "], max possible sizes in each dimension  = "
            << value_traits::Shape(TDataType{}) << " ].\n";

        DenseVector<unsigned int> tensor_shape(value_traits::Dimension + 1);
        tensor_shape[0] = rContainer.size();
        std::copy(rDataShapeBegin, rDataShapeEnd, tensor_shape.begin() + 1);

        return tensor_shape;

        KRATOS_CATCH("");
    }

    template <class TContainerType, class TDataType, class TGetterType>
    static auto GetTensorShape(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TGetterType& rGetter)
    {
        return TensorAdaptorUtils::GetTensorShape<TDataType>(rContainer, rGetter);
    }

    template <class TContainerType, class TDataType, class TIntegerType>
    static auto GetTensorShape(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        TIntegerType const * pDataShapeBegin,
        TIntegerType const * pDataShapeEnd)
    {
        return TensorAdaptorUtils::GetTensorShape<TDataType>(rContainer, pDataShapeBegin, pDataShapeEnd);
    }

    template <class TContainerType, class TDataType, class TSpanType, class TGetterType>
    static void CollectVariableData(
        const DenseVector<unsigned int>& rTensorShape,
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TSpanType& Span,
        const TGetterType& rGetter)
    {
        KRATOS_TRY

        CopyToContiguousArrayNew<TDataType>(
            rContainer, Span, rTensorShape.data().begin(),
            rTensorShape.data().begin() + rTensorShape.size(), rGetter);

        KRATOS_CATCH("");
    }

    template <class TContainerType, class TDataType, class TSpanType, class TSetterType>
    static void StoreVariableData(
        const DenseVector<unsigned int>& rTensorShape,
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TSpanType& Span,
        const TSetterType& rSetter)
    {
        KRATOS_TRY

        CopyFromContiguousDataArrayNew<TDataType>(
            rContainer, Span, rTensorShape.data().begin(),
            rTensorShape.data().begin() + rTensorShape.size(), rSetter);

        KRATOS_CATCH("");
    }

    template<class TContainerType>
    static std::string Info(
        const std::string& rPrefix,
        const DenseVector<unsigned int>& rTensorShape,
        const TContainerType& rContainer)
    {
        std::stringstream info;
        info << rPrefix << "number of " << ModelPart::Container<TContainerType>::GetEntityName()
             << "(s) = " << rContainer.size() << ", shape = " << rTensorShape << ".\n";
        return info.str();
    }

    ///@}
};

/// @}
} // namespace Kratos