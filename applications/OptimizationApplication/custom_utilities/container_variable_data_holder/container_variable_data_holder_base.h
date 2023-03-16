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
    virtual double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const = 0;

    virtual IndexType GetDimension() const = 0;
};

class LiteralDoubleExpression : public Expression
{
public:
    LiteralDoubleExpression(const double Value);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const double mValue;
};

class LiteralArray3Expression : public Expression
{
public:
    LiteralArray3Expression(const array_1d<double, 3>& Value, const IndexType Dimension);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const array_1d<double, 3> mValue;

    const IndexType mDimension;
};

class LiteralVectorExpression : public Expression
{
public:
    LiteralVectorExpression(std::shared_ptr<Vector> pValue, const IndexType Dimension);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const std::shared_ptr<Vector> mpValue;

    const IndexType mDimension;
};

class BinaryAddExpression : public Expression
{
public:
    BinaryAddExpression(std::shared_ptr<Expression> pLeft, std::shared_ptr<Expression> pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const std::shared_ptr<Expression> mpLeft;

    const std::shared_ptr<Expression> mpRight;
};

class BinarySubstractExpression : public Expression
{
public:
    BinarySubstractExpression(std::shared_ptr<Expression> pLeft, std::shared_ptr<Expression> pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const std::shared_ptr<Expression> mpLeft;

    const std::shared_ptr<Expression> mpRight;
};

class BinaryMultiplyExpression : public Expression
{
public:
    BinaryMultiplyExpression(std::shared_ptr<Expression> pLeft, std::shared_ptr<Expression> pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const std::shared_ptr<Expression> mpLeft;

    const std::shared_ptr<Expression> mpRight;
};

class BinaryDivideExpression : public Expression
{
public:
    BinaryDivideExpression(std::shared_ptr<Expression> pLeft, std::shared_ptr<Expression> pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const std::shared_ptr<Expression> mpLeft;

    const std::shared_ptr<Expression> mpRight;
};

class BinaryPowerExpression : public Expression
{
public:
    BinaryPowerExpression(std::shared_ptr<Expression> pLeft, std::shared_ptr<Expression> pRight);

    double Evaluate(const IndexType EntityIndex, const IndexType ComponentIndex) const override;

    IndexType GetDimension() const override;
private:
    const std::shared_ptr<Expression> mpLeft;

    const std::shared_ptr<Expression> mpRight;
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

    void SetExpression(std::shared_ptr<Expression> pExpression);

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

    std::shared_ptr<Expression> mpExpression = nullptr;

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