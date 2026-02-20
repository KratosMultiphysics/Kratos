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
#include "includes/kernel.h"
#include "containers/variable.h"
#include "includes/model_part.h"
#include "utilities/data_type_traits.h"
#include "utilities/container_io_utils.h"
#include "tensor_adaptors/tensor_adaptor.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @ingroup TensorAdaptors
 * @brief Utility class providing static functions for tensor adaptors.
 * @author Suneth Warnakulasuriya
 */
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

    /**
     * @brief Computes the shape of a tensor represented by a container of data elements.
     *
     * @details This static utility function determines the shape of a tensor based on the provided @p rContainer
     *          and a getter function. The resulting shape is returned as a DenseVector of unsigned integers,
     *          where the first dimension corresponds to the size of the @p rContainer, and the remaining dimensions
     *          are determined by the shape of the data value in the first entity retrieved by @p rGetter.
     *
     * @tparam TContainerType The type of the @p rContainer holding the data elements.
     * @tparam TGetterType The type of the getter function used to extract data from the @p rContainer elements.
     * @tparam TDataType The type of the data elements contained in the @p rContainer.
     *
     * @param rContainer The @p rContainer holding the data elements.
     * @param rGetter A callable object or function that extracts the data element from a @p rContainer entry.
     *
     * @return DenseVector<unsigned int> A vector representing the shape of the tensor. The first entry is the
     *         number of elements in the @p rContainer, followed by the shape of the data type in the first entity in the @p rContainer
     *         retrieved through @p rGetter function.
     *
     * @note If the @p rContainer is empty, the first dimension will be zero. For static data types, the remaining
     *       dimensions will reflect the correct static shape. For dynamic types (e.g., Vector, Matrix), the number of
     *       dimensions will be correct, but all values will be zero.
     *
     * @warning In distributed environments, if the data type is dynamic and the @p rContainer is empty on a rank,
     *          the computed shape may be incorrect in that rank. Therefore, it is recommended to use a predefined data shape in such cases,
     *          communicated from ranks with non-empty containers.
     */
    template<class TDataType, class TContainerType, class TGetterType>
    static DenseVector<unsigned int> GetTensorShape(
        const TContainerType& rContainer,
        const TGetterType& rGetter)
    {
        KRATOS_TRY

        using value_traits = DataTypeTraits<TDataType>;

        if (!rContainer.empty()) {

            TDataType dummy{};
            rGetter(dummy, rContainer.front());

            DenseVector<unsigned int> tensor_shape(value_traits::Dimension + 1);
            tensor_shape[0] = rContainer.size();
            DataTypeTraits<TDataType>::Shape(dummy, tensor_shape.data().begin() + 1, tensor_shape.data().begin() + tensor_shape.size());

            return tensor_shape;
        } else {
            // if there are not entities in the rContainer, return a tensor shape which has
            // first dimension with zero values, and other dimensions having the exact values
            // it the TDataType is a static one. Otherwise if the TDataType is a dynamic type such
            // as Vector or Matrix, the tensor_shape will have correct number of dimensionality,
            // but the values in all the dimensions will be zeros.
            DenseVector<unsigned int> tensor_shape(value_traits::Dimension + 1);
            value_traits::Shape(TDataType{}, tensor_shape.data().begin() + 1, tensor_shape.data().begin() + tensor_shape.size());

            if constexpr(value_traits::IsDynamic) {
                KRATOS_WARNING_IF("TensorAdaptorUtils::GetTensorShape", Kernel::IsDistributedRun())
                    << "The " << TDataType{} << " is a dynamic value type, and the method is called in an distributed environment"
                    << " which may lead to wrong tensor shape in this rank since this rank has an empty container."
                    << " Therefore, it is suggested to use the method with the predefined data shape in the distributed environemnt where the"
                    << " data shape is computed on ranks where there are non-empty containers and communicated to the empty container ranks.\n";
            }

            tensor_shape[0] = 0;
            return tensor_shape;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Computes the shape of a tensor represented by a container and its data shape.
     *
     * @details This static utility function constructs a DenseVector representing the shape of a tensor.
     *          The first entry of the returned vector is the size of the input @p rContainer (number of elements),
     *          followed by the dimensions specified by the data shape iterators [ @p rDataShapeBegin, @p rDataShapeEnd ).
     *
     * @throws std::runtime_error If the provided data shape cannot be used to represent the @p TDataType.
     *
     * @tparam TContainerType The type of the container holding the tensor data.
     * @tparam TDataType The type of the data stored in the container.
     * @tparam TIntegerType The integer type used for shape specification.
     * @param rContainer The container holding the tensor data.
     * @param rDataShapeBegin Iterator pointing to the beginning of the data shape array.
     * @param rDataShapeEnd Iterator pointing to the end of the data shape array.
     * @return DenseVector<unsigned int> A vector containing the tensor shape: [container size, data shape...].
     */
    template<class TDataType, class TContainerType, class TIntegerType>
    static DenseVector<unsigned int> GetTensorShape(
        const TContainerType& rContainer,
        TIntegerType const * rDataShapeBegin,
        TIntegerType const * rDataShapeEnd)
    {
        KRATOS_TRY

        using value_traits = DataTypeTraits<TDataType>;

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


    /**
     * @brief Retrieves the shape of a tensor from a given container using a specified getter and a variable.
     *
     * @details This static function determines the shape (dimensions) of a tensor stored within a @p rContainer.
     *          It utilizes a getter function or functor to access the tensor data associated with a specific variable.
     *          The @p rVariable is merely used to identify the underlying @p TDataType.
     *
     * @tparam TContainerType The type of the container holding the tensor data.
     * @tparam TDataType The type of the tensor data.
     * @tparam TGetterType The type of the getter used to access the tensor data.
     * @param rContainer Reference to the container holding the tensor data.
     * @param rVariable Reference to the variable describing the tensor data.
     * @param rGetter The getter used to access the tensor data from the container.
     */
    template <class TContainerType, class TDataType, class TGetterType>
    static auto GetTensorShape(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TGetterType& rGetter)
    {
        return TensorAdaptorUtils::GetTensorShape<TDataType>(rContainer, rGetter);
    }

    /**
     * @brief Retrieves the shape of a tensor from the given container and shape information along with the variable.
     *
     * @details This static function obtains the shape of a tensor associated with a specific @p rVariable
     *          from the provided @p rContainer. The shape is determined using the range specified by
     *          the pointers @p pDataShapeBegin and @p pDataShapeEnd. The @p rVariable is merely used to
     *          determine the @p TDataType.
     *
     * @tparam TContainerType Type of the container holding the tensor data.
     * @tparam TDataType Type of the tensor data.
     * @tparam TIntegerType Type used for shape indices (e.g., int, std::size_t).
     * @param rContainer Reference to the container holding the tensor data.
     * @param rVariable Reference to the variable describing the tensor.
     * @param pDataShapeBegin Pointer to the beginning of the shape data.
     * @param pDataShapeEnd Pointer to the end of the shape data.
     */
    template <class TContainerType, class TDataType, class TIntegerType>
    static auto GetTensorShape(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        TIntegerType const * pDataShapeBegin,
        TIntegerType const * pDataShapeEnd)
    {
        return TensorAdaptorUtils::GetTensorShape<TDataType>(rContainer, pDataShapeBegin, pDataShapeEnd);
    }

    /**
     * @brief Returns a zero-initialized value of the specified data type and shape.
     *
     * @details This static utility function creates a zero value of type @p TDataType, reshapes it according to the provided shape,
     *          and fills all its elements with zeros. It supports both static and dynamic types, utilizing the corresponding
     *          data type traits for construction, reshaping, and assignment.
     *
     * @note    This is different to @ref Variable::Zero method because, the @ref Variable::Zero method returns a zero sized @p TDataType
     *          values if the @p TDataType is dynamic such as @ref Vector or @ref Matrix. This method returns correctly sized (based on @p rDataShape )
     *          zero valued @p TDataType values for both static and dynamic types.
     *
     * @tparam TDataType The data type of the value to be created.
     * @param rVariable The variable describing the type of the value.
     * @param rDataShape The desired shape of the zero value, as a dense vector of unsigned integers.
     * @return TDataType A zero-initialized value of the specified type and shape.
     */
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

    /**
     * @brief Get the Nodal Neighbours Count Tensor Adaptor for nodes in specified model part
     * @details This method returns number of conditions (if @p TContainerType is of @ref ModelPart::ConditionsContainerType ) or
     *          number of elements (if @p TContainerType is of @ref ModelPart::ElementsContainerType ) around nodes in the @p rModelPart .
     *
     *          Returning @ref TensorAdaptor will be having the @p rModelPart nodes as the container.
     *
     * @tparam TContainerType
     * @param rModelPart                    Model part on which the nodal neighbour count will be done.
     * @return TensorAdaptor<int>::Pointer  Tensor adaptor containing the number of neighbours for each node.
     */
    template<class TContainerType>
    static TensorAdaptor<int>::Pointer GetNodalNeighboursCountTensorAdaptor(ModelPart& rModelPart);

    ///@}
};

/// @}
} // namespace Kratos