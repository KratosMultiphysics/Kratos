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
#include <vector>
#include <type_traits>

// Project includes
#include "expression/expression.h"
#include "expression/literal_flat_expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Construct class to read into expressions from templated data values and write in to templated data values.
 *
 * @tparam TDataType
 */
template<class TDataType>
class KRATOS_API(KRATOS_CORE) VariableExpressionDataIO
{
public:
    ///@name Type definitions
    ///@{

    using Pointer = std::shared_ptr<VariableExpressionDataIO<TDataType>>;

    using IndexType = std::size_t;

    using RawType = std::conditional_t<
                            std::disjunction_v<
                                std::is_same<TDataType, int>,
                                std::is_same<TDataType, std::vector<int>>
                            >, int, double>;

    using RawLiteralFlatExpression = LiteralFlatExpression<RawType>;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Variable Expression Data IO object from a given SampleValue
     *
     * @details This will assign correct shape information for static data types. For dynamic
     *          data types, it is required to have the correct shape in the provided SampleValue.
     *
     * @param SampleValue       Sample value to provide the shape information.
     */
    VariableExpressionDataIO(const TDataType& SampleValue);

    /**
     * @brief Construct a new Variable Expression Data IO object from a given shape
     *
     * @param rShape            Shape of the data to be transferred.
     */
    VariableExpressionDataIO(const std::vector<IndexType>& rShape);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Create a VariableExpressionDataIO from a given Sample value.
     *
     * @param SampleValue       Sample value.
     * @return Pointer          Pointer to a VariableExpressionDataIO.
     */
    static Pointer Create(const TDataType& SampleValue);

    /**
     * @brief Create a VariableExpressionDataIO from a given shape.
     *
     * @param rShape            Shape of the data type.
     * @return Pointer          Pointer to a VariableExpressionDataIO.
     */
    static Pointer Create(const std::vector<IndexType>& rShape);

    /**
     * @brief Evaluate the expression and assign the values to rOutput at the given index.
     *
     * @details This assigns the evaluated expression values to passed rOutput. The shape of the
     *          @ref rExpression and @ref rOutput should be checked before using this method.
     *
     * @throws  If the expression item shape and output shape mismatches in full debug.
     *
     * @param rOutput           Output containing evaluated rExpression values.
     * @param rExpression       Expression to be evaluated.
     * @param EntityIndex       Entity index at which the expression will be evaluated.
     */
    void Assign(
        TDataType& rOutput,
        const Expression& rExpression,
        const IndexType EntityIndex) const;

    /**
     * @brief Read into expression the values from given value.
     *
     * @details This assigns values at the given EntityIndex with the values from Value. The shape of the
     *          @ref rExpression and @ref Value should be checked before using this method.
     *
     * @throws  If the expression item shape and output shape mismatches if full debug.
     *
     * @param rExpression       Expression to be assigned.
     * @param EntityIndex       Entity of the expression to be assigned.
     * @param Value             Value to be assigned to expression.
     */
    void Read(
        RawLiteralFlatExpression& rExpression,
        const IndexType EntityIndex,
        const TDataType& Value) const;

    /**
     * @brief Get the shape of the data type
     *
     * @return const std::vector<IndexType>     Shape of the data type.
     */
    const std::vector<IndexType> GetItemShape() const { return mShape; }

    ///@}
private:
    ///@name Private members
    ///@{

    std::vector<IndexType> mShape;

    ///@}
};

} // namespace Kratos