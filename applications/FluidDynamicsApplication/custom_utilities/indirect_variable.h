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

#if !defined(KRATOS_FLUID_INDIRECT_VARIABLE_H)
#define KRATOS_FLUID_INDIRECT_VARIABLE_H

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template<class TDataType>
class IndirectData
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using TVariableType = Variable<TDataType>;

    using IndirectDataType = IndirectData<TDataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    IndirectData()
        : mSet([](const TDataType) {}),
          mGet([]() -> TDataType { return TDataType{}; })
    {};

    IndirectData(
        const TVariableType& rVariable,
        NodeType& rNode,
        const IndexType Step)
        : mSet([&rVariable, &rNode, Step](const TDataType Value) { rNode.FastGetSolutionStepValue(rVariable, Step) = Value;}),
          mGet([&rVariable, &rNode, Step]() -> TDataType { return rNode.FastGetSolutionStepValue(rVariable, Step); })
    {};

    ///@}
    ///@name Operators
    ///@{

    inline TDataType operator+(const TDataType T)
    {
        return mGet() + T;
    }

    inline TDataType operator-(const TDataType T)
    {
        return mGet() - T;
    }

    inline TDataType operator*(const TDataType T)
    {
        return mGet() * T;
    }

    inline TDataType operator/(const TDataType T)
    {
        return mGet() / T;
    }

    inline TDataType operator+(const IndirectDataType& T)
    {
        return mGet() + T.mGet();
    }

    inline TDataType operator-(const IndirectDataType& T)
    {
        return mGet() - T.mGet();
    }

    inline TDataType operator*(const IndirectDataType& T)
    {
        return mGet() * T.mGet();
    }

    inline TDataType operator/(const IndirectDataType& T)
    {
        return mGet() / T.mGet();
    }

    inline IndirectDataType& operator=(const TDataType T)
    {
        mSet(T);
        return *this;
    }

    inline IndirectDataType& operator+=(const TDataType T)
    {
        mSet(mGet() + T);
        return *this;
    }

    inline IndirectDataType& operator-=(const TDataType T)
    {
        mSet(mGet() - T);
        return *this;
    }

    inline IndirectDataType& operator*=(const TDataType T)
    {
        mSet(mGet() * T);
        return *this;
    }

    inline IndirectDataType& operator/=(const TDataType T)
    {
        mSet(mGet() / T);
        return *this;
    }

    inline IndirectDataType& operator=(const IndirectDataType& T)
    {
        mSet(T.mGet());
        return *this;
    }

    inline IndirectDataType& operator+=(const IndirectDataType& T)
    {
        mSet(mGet() + T.mGet());
        return *this;
    }

    inline IndirectDataType& operator-=(const IndirectDataType& T)
    {
        mSet(mGet() - T.mGet());
        return *this;
    }

    inline IndirectDataType& operator*=(const IndirectDataType& T)
    {
        mSet(mGet() * T.mGet());
        return *this;
    }

    inline IndirectDataType& operator/=(const IndirectDataType& T)
    {
        mSet(mGet() / T.mGet());
        return *this;
    }

    inline bool operator==(const TDataType T)
    {
        return mGet() == T;
    }

    inline bool operator!=(const TDataType T)
    {
        return mGet() != T;
    }

    inline bool operator>=(const TDataType T)
    {
        return mGet() >= T;
    }

    inline bool operator<=(const TDataType T)
    {
        return mGet() <= T;
    }

    inline bool operator>(const TDataType T)
    {
        return mGet() > T;
    }

    inline bool operator<(const TDataType T)
    {
        return mGet() < T;
    }

    inline bool operator==(const IndirectDataType& T)
    {
        return mGet() == T.mGet();
    }

    inline bool operator!=(const IndirectDataType& T)
    {
        return mGet() != T.mGet();
    }

    inline bool operator>=(const IndirectDataType& T)
    {
        return mGet() >= T.mGet();
    }

    inline bool operator<=(const IndirectDataType& T)
    {
        return mGet() <= T.mGet();
    }

    inline bool operator>(const IndirectDataType& T)
    {
        return mGet() > T.mGet();
    }

    inline bool operator<(const IndirectDataType& T)
    {
        return mGet() < T.mGet();
    }

    inline std::ostream& print(std::ostream& os) const
    {
        return os << mGet();
    }

    ///@}
private:
    ///@name Private Members
    ///@{

    const std::function<void(const TDataType)> mSet;
    const std::function<TDataType()> mGet;

    ///@}
};

template <class TDataType>
std::ostream& operator<<(std::ostream& os, const IndirectData<TDataType>& s)
{
    return s.print(os);
}

template<class TVariableDataType>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) IndirectVariable
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using IndirectDataType = IndirectData<TVariableDataType>;

    using TVariableType = Variable<TVariableDataType>;

    KRATOS_CLASS_POINTER_DEFINITION(IndirectVariable);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Indirect Scalar Variable object
     *
     * This constructor constructs an indirect scalar variable which
     * will return 0.0 for all the nodes, and will not set
     * any variable in the node if setter is used.
     *
     * This is usefull in the case where the value
     * is not required and wants to avoid branching
     * using if blocks and if wants to have generalized
     * formulations.
     *
     */
    IndirectVariable()
        : mrVariable(TVariableType::StaticObject()),
          mConstGetter(&IndirectVariable<TVariableDataType>::DefaultValueGetter),
          mIndirectData(&IndirectVariable<TVariableDataType>::DefaultIndirectDataGetter)
    {
    }

    /**
     * @brief Construct a new Indirect Scalar Variable object
     *
     * This constructor is used to construct an indirect variable
     * which can modify/retrieve nodal solution step data
     *
     * @param rVariable
     */

    IndirectVariable(
        const TVariableType& rVariable)
        : mrVariable(rVariable),
          mConstGetter(&IndirectVariable<TVariableDataType>::NodalValueGetter),
          mIndirectData(&IndirectVariable<TVariableDataType>::NodalIndirectDataGetter)
    {
    }

    ///@}
    ///@name Public Member Operations
    ///@{

    inline const TVariableType& GetVariable() const
    {
        return mrVariable;
    }

    template <class TUnaryFunction>
    inline void ExecuteForVariable(TUnaryFunction &&f) const
    {
        if (mrVariable != TVariableType::StaticObject()) {
            f(mrVariable);
        }
    }

    ///@}
    ///@name Operators
    ///@{

    inline IndirectDataType operator()(NodeType& rNode) const
    {
        return  (this->*(this->mIndirectData))(rNode, 0);
    }

    inline TVariableDataType operator()(const NodeType& rNode) const
    {
        return  (this->*(this->mConstGetter))(rNode, 0);
    }

    inline IndirectDataType operator()(NodeType& rNode, const IndexType Step) const
    {
        return  (this->*(this->mIndirectData))(rNode, Step);
    }

    inline TVariableDataType operator()(const NodeType& rNode, const IndexType Step) const
    {
        return  (this->*(this->mConstGetter))(rNode, Step);
    }

    ///@}
private:
    ///@name Private Members
    ///@{

    const TVariableType& mrVariable;

    TVariableDataType (IndirectVariable<TVariableDataType>::*mConstGetter)(
        const NodeType&,
        const IndexType) const;

    IndirectDataType (IndirectVariable<TVariableDataType>::*mIndirectData)(
        NodeType&,
        const IndexType) const;

    // nodal value methods
    inline TVariableDataType NodalValueGetter(
        const NodeType& rNode,
        const IndexType Step) const
    {
        return rNode.FastGetSolutionStepValue(mrVariable, Step);
    }

    inline IndirectDataType NodalIndirectDataGetter(
        NodeType& rNode,
        const IndexType Step) const
    {
        return IndirectDataType(mrVariable, rNode, Step);
    }

    // default value methos
    inline TVariableDataType DefaultValueGetter(
        const NodeType& rNode,
        const IndexType Step) const
    {
        return TVariableDataType{};
    }

    inline IndirectDataType DefaultIndirectDataGetter(
        NodeType& rNode,
        const IndexType Step) const
    {
        return IndirectDataType();
    }


    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_INDIRECT_VARIABLE_H
