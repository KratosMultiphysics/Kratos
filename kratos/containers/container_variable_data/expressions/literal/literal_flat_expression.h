//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/container_variable_data/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Expression to hold a literal with a flattened data structure
 *
 * This expression holds data of a literal with rShape by
 * flattening the data structure to a vector. This allocates
 * a Vector with large sizes.
 *
 */
class LiteralFlatExpression : public Expression {
public:
    ///@name Type definition
    ///@{

    using Pointer = Kratos::intrusive_ptr<LiteralFlatExpression>;

    ///@}
    ///@name Life cycle
    ///@{

    LiteralFlatExpression(
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    LiteralFlatExpression(const LiteralFlatExpression& rOther) = delete;

    ~LiteralFlatExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief A specific create method is used in here to return a pointer to the LiteralFlatExpression
     *
     * This method is used to return the derrived class LiteralFlatExpression pointer
     * so that additional SetData method can be accessed to initialize
     * the vector data.
     *
     * @param NumberOfEntities                  Number of entities in the container.
     * @param rShape                            Shape of the data in each entitiy.
     * @return LiteralFlatExpression::Pointer   Returns an intrusive pointer to LiteralFlatExpression.
     */
    static LiteralFlatExpression::Pointer Create(
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    void SetData(
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex,
        const double Value);

    double Evaluate(
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::vector<IndexType> mShape;

    Vector mData;

    ///@}
};

} // namespace Kratos