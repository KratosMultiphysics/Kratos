//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_DERIVATIVES_RECOVERY_UTILITY_H_INCLUDED
#define KRATOS_DERIVATIVES_RECOVERY_UTILITY_H_INCLUDED

// System includes
#include <unordered_set>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// @brief Forward declaration of ModelPart
class ModelPart;

/**
 * @ingroup ShallowWaterApplication
 * @class DerivativesRecoveryUtility
 * @brief Superconvergent patch recovery for linear meshes using quadratic polynomials
 * @details Zhimin Zhangand and Ahmed Naga (2006)
 * "A new finite element gradient recovery method: superconvergence property"
 * SIAM J. Sci. Comput., 26(4), 1192–1213. (22 pages)
 *
 * The polynomial follows the next convention:
 *    p_2(x, y, z, Node) = trans(P)*a
 *    trans(P) = (1, x, y, z, x^2, y^2, z^2, xy, xz, yz)
 *    trans(a) = a_i
 */
template<std::size_t TDim>
class DerivativesRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DerivativesRecoveryUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static void Check(ModelPart& rModelPart);

    static void ExtendNeighborsPatch(ModelPart& rModelPart);

    static void CalculatePolynomialWeights(ModelPart& rModelPart);

    static void RecoverDivergence(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    static void RecoverGradient(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable,
        const Variable<array_1d<double,3>>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    static void RecoverLaplacian(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    static void RecoverLaplacian(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rOriginVariable,
        const Variable<array_1d<double,3>>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static double CalculateMaximumDistance(
        const NodeType& rNode,
        GlobalPointersVector<NodeType>& rNeighbors);

    static void FindExtendedNeighbors(
        NodeType& rNode,
        GlobalPointersVector<NodeType>& rNeighbors,
        std::unordered_set<int>& rExtendedNeighborsId);

    static void AppendExtendedNeighbors(
        ModelPart& rModelPart,
        GlobalPointersVector<NodeType>& rNeighbors,
        std::unordered_set<int>& rExtendedNeighborsId);

    static bool CalculateNodalPolynomialWeights(
        NodeType& rNode);

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class DerivativesRecoveryUtility

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_DERIVATIVES_RECOVERY_UTILITY_H_INCLUDED defined
