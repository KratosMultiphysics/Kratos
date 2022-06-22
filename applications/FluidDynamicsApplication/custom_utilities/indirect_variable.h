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

template<class TVariableType>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) IndirectVariable
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

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
        : mNonConstGetter([&](NodeType&, const IndexType) -> TVariableType& { mTemp = TVariableType{}; return mTemp;}),
          mConstGetter([&](const NodeType&, const IndexType) -> TVariableType { mTemp = TVariableType{}; return mTemp;})
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
        const Variable<TVariableType>& rVariable)
        : mNonConstGetter([&rVariable](NodeType& rNode, const IndexType Step) -> TVariableType& { return rNode.FastGetSolutionStepValue(rVariable, Step);}),
          mConstGetter([&rVariable](const NodeType& rNode, const IndexType Step) -> TVariableType { return rNode.FastGetSolutionStepValue(rVariable, Step);})
    {
    }

    ///@}
    ///@name Operators
    ///@{

    TVariableType& operator()(NodeType& rNode)
    {
        return mNonConstGetter(rNode, 0);
    }

    TVariableType operator()(const NodeType& rNode) const
    {
        return mConstGetter(rNode, 0);
    }

    TVariableType& operator()(NodeType& rNode, const IndexType Step)
    {
        return mNonConstGetter(rNode, Step);
    }

    TVariableType operator()(const NodeType& rNode, const IndexType Step) const
    {
        return mConstGetter(rNode, Step);
    }

    ///@}
private:
    ///@name Private Members
    ///@{

    TVariableType mTemp{};

    const std::function<TVariableType&(NodeType&, const IndexType)> mNonConstGetter;
    const std::function<TVariableType(const NodeType&, const IndexType)> mConstGetter;

    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_INDIRECT_VARIABLE_H
