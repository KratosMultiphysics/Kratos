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
#include <string>
#include <variant>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

class HistoricalContainerDataIO
{
public:
    ///@name Public operations
    ///@{

    template<class TDataType>
    static TDataType& GetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable)
    {
        return rNode.FastGetSolutionStepValue(rVariable);
    }

    template<class TDataType>
    static void SetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rNode.FastGetSolutionStepValue(rVariable) = rValue;
    }

    ///@}
};

class NonHistoricalContainerDataIO
{
public:
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        return rEntity.GetValue(rVariable);
    }

    template<class TEntityType, class TDataType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.SetValue(rVariable, rValue);
    }

    ///@}
};

class PropertiesContainerDataIO
{
public:
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    static TDataType& GetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable)
    {
        return rEntity.GetProperties().GetValue(rVariable);
    }

    template<class TEntityType, class TDataType>
    static void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rEntity.GetProperties().SetValue(rVariable, rValue);
    }

    ///@}
};


class Expression
{
public:
    using Pointer = Kratos::intrusive_ptr<Expression>;

    virtual ~Expression() = default;

    virtual double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const = 0;

    virtual IndexType GetDimension() const = 0;

private:
    //*********************************************
    //this block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const Expression* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const Expression* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
        }
    }
    //*********************************************
};

class LiteralDoubleExpression : public Expression
{
public:
    LiteralDoubleExpression(const double Value);

    LiteralDoubleExpression(const LiteralDoubleExpression& rOther) = delete;

    ~LiteralDoubleExpression() override = default;

    static Expression::Pointer Create(const double Value);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const double mValue;
};

class LiteralArray3Expression : public Expression
{
public:
    LiteralArray3Expression(const array_1d<double, 3>& Value, const IndexType Dimension);

    LiteralArray3Expression(const LiteralArray3Expression& rOther) = delete;

    ~LiteralArray3Expression() override = default;

    static Expression::Pointer Create(const array_1d<double, 3>& Value, const IndexType Dimension);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const array_1d<double, 3> mValue;

    const IndexType mDimension;
};

class LiteralVectorExpression : public Expression
{
public:
    LiteralVectorExpression(Kratos::shared_ptr<Vector> pValue, const IndexType Dimension);

    LiteralVectorExpression(const LiteralVectorExpression& rOther) = delete;

    ~LiteralVectorExpression() override = default;

    static Expression::Pointer Create(Kratos::shared_ptr<Vector> pValue, const IndexType Dimension);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const Kratos::shared_ptr<Vector> mpValue;

    const IndexType mDimension;
};

class BinaryAddExpression : public Expression
{
public:
    BinaryAddExpression(Expression::Pointer pLeft, Expression::Pointer pRight);

    BinaryAddExpression(const BinaryAddExpression& rOther) = delete;

    ~BinaryAddExpression() override = default;

    static Expression::Pointer Create(Expression::Pointer pLeft, Expression::Pointer pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;
};

class BinarySubstractExpression : public Expression
{
public:
    BinarySubstractExpression(Expression::Pointer pLeft, Expression::Pointer pRight);

    BinarySubstractExpression(const BinarySubstractExpression& rOther) = delete;

    ~BinarySubstractExpression() override = default;

    static Expression::Pointer Create(Expression::Pointer pLeft, Expression::Pointer pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;
};

class BinaryMultiplyExpression : public Expression
{
public:
    BinaryMultiplyExpression(Expression::Pointer pLeft, Expression::Pointer pRight);

    BinaryMultiplyExpression(const BinaryMultiplyExpression& rOther) = delete;

    ~BinaryMultiplyExpression() override = default;

    static Expression::Pointer Create(Expression::Pointer pLeft, Expression::Pointer pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;
};

class BinaryDivideExpression : public Expression
{
public:
    BinaryDivideExpression(Expression::Pointer pLeft, Expression::Pointer pRight);

    BinaryDivideExpression(const BinaryDivideExpression& rOther) = delete;

    ~BinaryDivideExpression() override = default;

    static Expression::Pointer Create(Expression::Pointer pLeft, Expression::Pointer pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;
};

class BinaryPowerExpression : public Expression
{
public:
    BinaryPowerExpression(Expression::Pointer pLeft, Expression::Pointer pRight);

    BinaryPowerExpression(const BinaryPowerExpression& rOther) = delete;

    ~BinaryPowerExpression() override = default;

    static Expression::Pointer Create(Expression::Pointer pLeft, Expression::Pointer pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;
};

template <class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataHolderBase {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(ContainerVariableDataHolderBase);

    ///@}
    ///@name Life cycle
    ///#{

    virtual ~ContainerVariableDataHolderBase() = default;

    ///@}
    ///@name Public operations
    ///@{

    void CopyDataFrom(const ContainerVariableDataHolderBase<TContainerType>& rOther);

    void SetDataToZero();

    ///@}
    ///@name Input and output
    ///@{

    void SetExpression(Expression::Pointer pExpression);

    const Expression& GetExpression() const;

    IndexType GetDataDimension() const;

    ModelPart& GetModelPart();

    const ModelPart& GetModelPart() const;

    TContainerType& GetContainer();

    const TContainerType& GetContainer() const;

    virtual std::string Info() const;

    std::string PrintData() const;

    ///@}
protected:
    ///@name Life cycle
    ///@{

    /// Constructor with the model part
    ContainerVariableDataHolderBase(ModelPart& rModelPart);

    /// Copy constructor
    ContainerVariableDataHolderBase(const ContainerVariableDataHolderBase& rOther);

    ///@}
    ///@name Protected member variables
    ///@{

    Expression::Pointer mpExpression = nullptr;

    ModelPart* const mpModelPart;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerVariableDataHolderBase<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}

} // namespace Kratos