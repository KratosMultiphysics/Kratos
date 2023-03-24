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
     * @brief Get the Shape of the expression
     *
     * @return const std::vector<IndexType>     Size of each dimension is in the vector elements.
     */
    virtual const std::vector<IndexType> GetShape() const = 0;

    /**
     * @brief Get the Local Size of the expression
     *
     * @return IndexType
     */
    IndexType GetLocalSize() const;

    ///@}
    ///@name Input and output
    ///@{

    virtual std::string Info() const = 0;

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

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

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

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

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
        const std::vector<IndexType>& rShape);

    LiteralVectorExpression(const LiteralVectorExpression& rOther) = delete;

    ~LiteralVectorExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Kratos::shared_ptr<Vector> pValue,
        const std::vector<IndexType>& rShape);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const Kratos::shared_ptr<const Vector> mpValue;

    const std::vector<IndexType> mShape;

    IndexType mLocalSize;

    ///@}
};

class BinaryExpression: public Expression
{
public:
    ///@name Life cycle
    ///@{

    BinaryExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinaryExpression(const BinaryExpression& rOther) = delete;

    ~BinaryExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    const std::vector<IndexType> GetShape() const override;

    virtual std::string Operation() const = 0;

    std::string Info() const override;

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    ///@}
};

/**
 * @brief Expression to add two given expressions.
 *
 */
class BinaryAddExpression : public BinaryExpression {
public:
    ///@name Base class exposures
    ///@{

    using BinaryExpression::BinaryExpression;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    std::string Operation() const override { return "+"; }

    ///@}
};

/**
 * @brief Expression to substract two given expressions.
 *
 */
class BinarySubstractExpression : public BinaryExpression {
public:
    ///@name Base class exposures
    ///@{

    using BinaryExpression::BinaryExpression;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    std::string Operation() const override { return "-"; }

    ///@}
};

/**
 * @brief Expression to multiply two given expressions.
 *
 */
class BinaryMultiplyExpression : public BinaryExpression {
public:
    ///@name Base class exposures
    ///@{

    using BinaryExpression::BinaryExpression;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    std::string Operation() const override { return "*"; }

    ///@}
};

/**
 * @brief Expression to divide two given expressions.
 *
 */
class BinaryDivideExpression : public BinaryExpression {
public:
    ///@name Base class exposures
    ///@{

    using BinaryExpression::BinaryExpression;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    std::string Operation() const override { return "/"; }

    ///@}
};

/**
 * @brief Expression to raise one expression to the power of other expression.
 *
 */
class BinaryPowerExpression : public BinaryExpression {
public:
    ///@name Base class exposures
    ///@{

    using BinaryExpression::BinaryExpression;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    std::string Operation() const override { return "^"; }

    ///@}
};

/**
 * @brief Expression to compute weighted multiplication.
 *
 */
class BinaryWeightedMultiplicationExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    BinaryWeightedMultiplicationExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinaryWeightedMultiplicationExpression(const BinaryWeightedMultiplicationExpression& rOther) = delete;

    ~BinaryWeightedMultiplicationExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    ///@}
};

/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const Expression& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos