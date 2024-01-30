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
#include <iterator>

// Project includes
#include "expression/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Used to create lazy expression to reshape given expression to a new shape.
 *
 */
class UnaryReshapeExpression : public Expression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TIteratorType>
    UnaryReshapeExpression(
        Expression::ConstPointer pExpression,
        TIteratorType Begin,
        TIteratorType End)
        : Expression(pExpression->NumberOfEntities()),
        mpSourceExpression(pExpression),
        mShape(Begin, End)
    {
        KRATOS_ERROR_IF_NOT(this->GetItemComponentCount() == mpSourceExpression->GetItemComponentCount())
            << "The source and destination should have same stride. [ source stride = "
            << mpSourceExpression->GetItemComponentCount() << ", destination stride = "
            << this->GetItemComponentCount() << " ].\n"
            << "Source expression:\n"
            << *mpSourceExpression << "\n";
    }

    ///@}
    ///@name Public operations
    ///@{

    template<class TIteratorType>
    static Expression::Pointer Create(
        Expression::ConstPointer pExpression,
        TIteratorType Begin,
        TIteratorType End)
    {
        return Kratos::make_intrusive<UnaryReshapeExpression>(std::move(pExpression), Begin, End);
    }

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override
    {
        return mpSourceExpression->Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex);
    }

    const std::vector<IndexType> GetItemShape() const override
    {
        return mShape;
    }

    IndexType GetMaxDepth() const override
    {
        return mpSourceExpression->GetMaxDepth() + 1;
    }

    std::string Info() const override
    {
        std::stringstream msg;
        msg << mpSourceExpression << "-> [";
        for (const auto v : mShape) {
            msg << " " << v;
        }
        msg << " ]";
        return msg.str();
    }

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const Expression::ConstPointer mpSourceExpression;

    const std::vector<IndexType> mShape;

    ///@}
};

} // namespace Kratos