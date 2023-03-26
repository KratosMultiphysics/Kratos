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

// Application includes
#include "custom_utilities/container_variable_data/expressions/expression.h"

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
    ///@name Life cycle
    ///@{

    LiteralFlatExpression(const IndexType NumberOfEntities,
                          const std::vector<IndexType>& rShape);

    LiteralFlatExpression(const LiteralFlatExpression& rOther) = delete;

    ~LiteralFlatExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

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