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

    template<class TDataType>
    [[nodiscard]] static TDataType GetZeroValue(
        const Variable<TDataType>& rVariable,
        const DenseVector<unsigned int>& rDataShape)
    {
        using data_type_traits = DataTypeTraits<TDataType>;

        using primitive_type = typename data_type_traits::PrimitiveType;

        // first create a zero value.
        TDataType zero{};

        // reshape it to correct size. This will not do anything if it is a static type
        // but in the case of dynamic type, this will reshape it to given size.
        data_type_traits::Reshape(zero, rDataShape.data().begin(), rDataShape.data().begin() + rDataShape.size());

        // now assign zeros to all the values.
        std::vector<primitive_type> zeros(data_type_traits::Size(zero), primitive_type{});
        data_type_traits::CopyFromContiguousData(zero, zeros.data());
        return zero;
    }

    ///@}
};

/// @}
} // namespace Kratos