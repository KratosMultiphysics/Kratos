//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#if !defined(KRATOS_VARIABLE_CHECK_UTILITIES_H_INCLUDED)
#define KRATOS_VARIABLE_CHECK_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/node.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Methods for common operations in Element and Condition Check functions.
class CheckUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CheckUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CheckUtilities);

    ///@}
    ///@name Operations
    ///@{

    /** @brief Check if a variable is registered in Kratos.
    *   @param rVar The variable to be checked.
    *   Raises Kratos::Exception if the check fails.
    */
    static void CheckVariableKey(const Kratos::VariableData& rVar)
    {
        KRATOS_ERROR_IF(rVar.Key() == 0)
            << rVar.Name() << " Key is 0." << std::endl
            << "Check that Kratos variables have been correctly registered and "
               "all required applications have been imported."
            << std::endl;
    }

    /** @brief Check if a variable has been added to a node's solution step data
    * container (regular variable version).
    *   @param rVar The variable to be checked.
    *   @param rNode The node to be checked.
    *   Raises Kratos::Exception if the check fails.
    */
    template <class TValueType>
    static void CheckVariableInNodalData(const Kratos::Variable<TValueType>& rVar,
                                         Node<3>& rNode)
    {
        KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rVar))
            << "Missing " << rVar.Name()
            << " variable in solution step data for node " << rNode.Id() << "."
            << std::endl;
    }

    /** @brief Check if a variable has been added to a node's solution step data
    * container (variable component version).
    *   @param rVar The variable to be checked.
    *   @param rNode The node to be checked.
    *   Raises Kratos::Exception if the check fails.
    */
    template <class TAdaptorType>
    static void CheckVariableInNodalData(const Kratos::VariableComponent<TAdaptorType>& rVar,
                                         Node<3>& rNode)
    {
        KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rVar))
            << "Missing " << rVar.Name()
            << " variable in solution step data for node " << rNode.Id() << "."
            << std::endl;
    }

    /** @brief Check if a variable has been defined as Degree of Freedom for a
    * given node.
    *   @param rVar The variable to be checked.
    *   @param rNode The node to be checked.
    *   Raises Kratos::Exception if the check fails.
    */
    static void CheckDofInNode(const Kratos::VariableData& rVar, Node<3>& rNode)
    {
        KRATOS_ERROR_IF_NOT(rNode.HasDofFor(rVar))
            << "Missing Degree of Freedom for " << rVar.Name() << " in node "
            << rNode.Id() << "." << std::endl;
    }

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    CheckUtilities();

    /// Assignment operator.
    CheckUtilities& operator=(CheckUtilities const& rOther);

    /// Copy constructor.
    CheckUtilities(CheckUtilities const& rOther);

    ///@}

}; // Class CheckUtilities

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_VARIABLE_CHECK_UTILITIES_H_INCLUDED  defined
