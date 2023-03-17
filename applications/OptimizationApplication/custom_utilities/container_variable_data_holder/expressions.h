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
#include <cmath>
#include <atomic>
#include <vector>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the expression types.
 *
 * This is used to represent an expression for arithmatic evaluation
 *
 */
class Expression {
public:
    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<Expression>;

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    virtual ~Expression() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Evalute the expression for the given entity index and component index and returns the value
     *
     * @param EntityIndex       Entity index
     * @param ComponentIndex    Component index
     * @return double           Evaluated expression
     */
    virtual double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const = 0;

    /**
     * @brief Get the Dimension of the expression
     *
     * @return IndexType    Dimension of the expression
     */
    virtual IndexType GetDimension() const = 0;

    ///@}

private:
    ///@name Private operations
    ///@{

    //*********************************************
    // this block is needed for refcounting in th eintrusive ptr
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const Expression* x);

    friend void intrusive_ptr_release(const Expression* x);
    //*********************************************

    ///@}
};

/**
 * @brief Expression to hold a double value
 *
 */
class LiteralDoubleExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    LiteralDoubleExpression(const double Value);

    LiteralDoubleExpression(const LiteralDoubleExpression& rOther) = delete;

    ~LiteralDoubleExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(const double Value);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mValue;

    ///@}
};

/**
 * @brief Expression to hold a 3D array value
 *
 */
class LiteralArray3Expression : public Expression {
public:
    ///@name Life cycle
    ///@{

    LiteralArray3Expression(
        const array_1d<double, 3>& Value,
        const IndexType Dimension);

    LiteralArray3Expression(const LiteralArray3Expression& rOther) = delete;

    ~LiteralArray3Expression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        const array_1d<double, 3>& Value,
        const IndexType Dimension);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const array_1d<double, 3> mValue;

    const IndexType mDimension;

    ///@}
};

/**
 * @brief Expression to hold a vector
 *
 */
class LiteralVectorExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    LiteralVectorExpression(
        Kratos::shared_ptr<const Vector> pValue,
        const IndexType Dimension);

    LiteralVectorExpression(const LiteralVectorExpression& rOther) = delete;

    ~LiteralVectorExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Kratos::shared_ptr<Vector> pValue,
        const IndexType Dimension);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Kratos::shared_ptr<const Vector> mpValue;

    const IndexType mDimension;

    ///@}
};

/**
 * @brief Expression to add two given expressions.
 *
 */
class BinaryAddExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    BinaryAddExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinaryAddExpression(const BinaryAddExpression& rOther) = delete;

    ~BinaryAddExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    ///@}
};

/**
 * @brief Expression to substract two given expressions.
 *
 */
class BinarySubstractExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    BinarySubstractExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinarySubstractExpression(const BinarySubstractExpression& rOther) = delete;

    ~BinarySubstractExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    ///@}
};

/**
 * @brief Expression to multiply two given expressions.
 *
 */
class BinaryMultiplyExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    BinaryMultiplyExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinaryMultiplyExpression(const BinaryMultiplyExpression& rOther) = delete;

    ~BinaryMultiplyExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    ///@}
};

/**
 * @brief Expression to divide two given expressions.
 *
 */
class BinaryDivideExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    BinaryDivideExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinaryDivideExpression(const BinaryDivideExpression& rOther) = delete;

    ~BinaryDivideExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    ///@}
};

/**
 * @brief Expression to raise one expression to the power of other expression.
 *
 */
class BinaryPowerExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    BinaryPowerExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinaryPowerExpression(const BinaryPowerExpression& rOther) = delete;

    ~BinaryPowerExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    ///@}
};

} // namespace Kratos