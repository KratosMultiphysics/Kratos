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
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"

// Application includes

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos classes
///@{

class KRATOS_API(KRATOS_CORE) IndirectScalarVariable
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    KRATOS_CLASS_POINTER_DEFINITION(IndirectScalarVariable);

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
    IndirectScalarVariable()
        : mrVariable(Variable<double>::StaticObject()),
          mIsZeroVariable(true)
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
    IndirectScalarVariable(
        const Variable<double>& rVariable)
        : mrVariable(rVariable),
          mIsZeroVariable(false)
    {
    }

    ///@}
    ///@name Public Member Operations
    ///@{

    inline const Variable<double>& GetVariable() const
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mIsZeroVariable) << "Calling IndirectScalarVariable::GetVariable for Zero variable is not allowed.\n";
        return mrVariable;

        KRATOS_CATCH("");
    }

    inline bool IsZeroVariable() const
    {
        return mIsZeroVariable;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        if (mIsZeroVariable) {
            buffer << "Indirect variable which is always zero.";
        } else {
            buffer << "Indirect variable for " << mrVariable;
        }
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    ///@}
    ///@name Operators
    ///@{

    inline double& operator()(NodeType& rNode) const
    {
        mDefaultValue = 0.0;
        return  mIsZeroVariable ? mDefaultValue : rNode.FastGetSolutionStepValue(mrVariable);
    }

    inline double operator()(const NodeType& rNode) const
    {
        return  mIsZeroVariable ? 0.0 : rNode.FastGetSolutionStepValue(mrVariable);
    }

    inline double operator()(const NodeType& rNode, const IndexType Step) const
    {
        return  mIsZeroVariable ? 0.0 : rNode.FastGetSolutionStepValue(mrVariable, Step);
    }

    inline double& operator()(NodeType& rNode, const IndexType Step) const
    {
        mDefaultValue = 0.0;
        return  mIsZeroVariable ? mDefaultValue : rNode.FastGetSolutionStepValue(mrVariable, Step);
    }

    ///@}
private:
    ///@name Private static members
    ///@{

    // these members are required for objects created with the default constructor
#ifdef KRATOS_SMP_OPENMP
    #ifdef KRATOS_COMPILED_IN_WINDOWS
        __declspec(thread) static double mDefaultValue;
    #else
        static double mDefaultValue;
        #pragma omp threadprivate(mDefaultValue)
    #endif
#elif defined(KRATOS_SMP_CXX11)
    static thread_local mDefaultValue;
#else
    static double mDefaultValue;
#endif

    ///@}
    ///@name Private local members
    ///@{

    const Variable<double>& mrVariable;

    const bool mIsZeroVariable;

    ///@}
};

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  IndirectScalarVariable& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const IndirectScalarVariable& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

///@}

///@}

} // namespace Kratos