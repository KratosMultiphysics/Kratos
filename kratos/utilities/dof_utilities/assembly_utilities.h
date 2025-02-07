//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/sparse_graph.h"
#include "includes/define.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class AssemblyUtilities
 * @ingroup KratosCore
 * @brief Utility class for handling the build
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
class AssemblyUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Build type definition
    enum class BuildType
    {
        Block,
        Elimination
    };

    /// Pointer definition of AssemblyUtilities
    KRATOS_CLASS_POINTER_DEFINITION(AssemblyUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssemblyUtilities() = delete;

    ///@}
    ///@name Operations
    ///@{

    static void SetUpSparseGraph(
        const ModelPart& rModelPart,
        SparseGraph<typename ModelPart::IndexType>& rSparseGraph)
    {
        if (!rSparseGraph.IsEmpty()) {
            KRATOS_WARNING("AssemblyUtilities") << "Provided sparse graph is not empty and will be cleared." << std::endl;
            rSparseGraph.Clear();
        }

        Element::EquationIdVectorType eq_ids;
        for (auto& r_elem : rModelPart.Elements()) {
            r_elem.EquationIdVector(eq_ids, rModelPart.GetProcessInfo());
            rSparseGraph.AddEntries(eq_ids);
        }
        for (auto& r_cond : rModelPart.Conditions()) {
            r_cond.EquationIdVector(eq_ids, rModelPart.GetProcessInfo());
            rSparseGraph.AddEntries(eq_ids);
        }
    }

    ///@}
    ///@name Input and output
    ///@{


    ///@}
private:
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
}; // Class AssemblyUtilities

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
