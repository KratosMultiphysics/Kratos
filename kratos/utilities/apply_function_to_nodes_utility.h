//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_APPLY_FUNCTION_TO_NODES_UTILITY_H_INCLUDED)
#define  KRATOS_APPLY_FUNCTION_TO_NODES_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ApplyFunctionToNodesUtility
 * @ingroup KratosCore
 * @brief This function applies a givn function to its nodes calling GenericFunctionUtility
 * @details The functions can be constructed by providing a python-defined method of the type
 * @author Riccardo Rossi
 */
class ApplyFunctionToNodesUtility
{
public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// Counted pointer of ApplyFunctionToNodesUtility
    KRATOS_CLASS_POINTER_DEFINITION(ApplyFunctionToNodesUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rNodes The nodes where to set the values
     * @param pFunction The function to set
     */

    ApplyFunctionToNodesUtility(
        ModelPart::NodesContainerType& rNodes,
        GenericFunctionUtility::Pointer pFunction
        ): mrNodes(rNodes),
           mpFunction(pFunction)
    {

    };

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method applies the function for a given variable at a certain time
     * @param rVariable The variable type
     * @param t The time variable
     */
    void ApplyFunction(
        const Variable<double>& rVariable,
        const double t
        )
    {
        // The first node iterator
        const auto it_node_begin = mrNodes.begin();

        if(!mpFunction->UseLocalSystem()) {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->CallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                it_node->FastGetSolutionStepValue(rVariable) = value;
            }
        } else {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->RotateAndCallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                it_node->FastGetSolutionStepValue(rVariable) = value;
            }
        }
    }

    /**
     * @brief This method returns all the evaluated values in a given time
     * @param t The time variable
     */
    std::vector<double> ReturnFunction(const double t)
    {
        // The first node iterator
        const auto it_node_begin = mrNodes.begin();

        // The vector containing the values
        std::vector<double> values(mrNodes.size());

        //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
        if(!mpFunction->UseLocalSystem()) {
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->CallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                values[k] = value;
            }
        } else {
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->RotateAndCallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                values[k] = value;
            }
        }

        return values;
    }

private:

    ///@name Protected member Variables
    ///@{

    ModelPart::NodesContainerType& mrNodes;           /// The nodes where to set the function
    GenericFunctionUtility::Pointer mpFunction; /// The function to set

    ///@}
}; /// ApplyFunctionToNodesUtility

} /// namespace Kratos

#endif // KRATOS_APPLY_FUNCTION_TO_NODES_UTILITY_H_INCLUDED
