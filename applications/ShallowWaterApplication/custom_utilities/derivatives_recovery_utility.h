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
class KRATOS_API(SHALLOW_WATER_APPLICATION) DerivativesRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node NodeType;

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

    /**
     * @brief Check that all the required variables are present in the nodal database.
     * @param rModelPart 
     */
    static void Check(ModelPart& rModelPart);

    /**
     * @brief Check that all the nodes have the required number of neighbors (6 in 2D and 10 in 3D)
     * @param rModelPart 
     */
    static void CheckRequiredNeighborsPatch(ModelPart& rModelPart);

    /**
     * @brief Same as CheckRequiredNeighborsPatch with a custom number of required neighbors
     * @param rModelPart 
     * @param RequiredNeighbors 
     */
    static void ExtendNeighborsPatch(ModelPart& rModelPart, const std::size_t RequiredNeighbors);

    /**
     * @brief Fit a polynomial of degree=2 and store its derivative coefficients in FIRST_DERIVATIVE_WEIGHTS and SECOND_DERIVATIVE_WEIGHTS
     * @details It calls CheckRequiredNeighborsPatch
     * @see CheckRequiredNeighborsPatch
     * @param rModelPart 
     */
    static void CalculatePolynomialWeights(ModelPart& rModelPart);

    /**
     * @brief Recover the divergence of a vector variable
     * @param rModelPart 
     * @param rOriginVariable 
     * @param rDestinationVariable 
     * @param BufferStep 
     */
    static void RecoverDivergence(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    /**
     * @brief Recover the gradient of a scalar variable
     * @param rModelPart 
     * @param rOriginVariable 
     * @param rDestinationVariable 
     * @param BufferStep 
     */
    static void RecoverGradient(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable,
        const Variable<array_1d<double,3>>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    /**
     * @brief Recover the laplacian of a scalar variable
     * @param rModelPart 
     * @param rOriginVariable 
     * @param rDestinationVariable 
     * @param BufferStep 
     */
    static void RecoverLaplacian(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    /**
     * @brief Recover the laplacian of a vector variable as the gradient of the divergence
     * @param rModelPart 
     * @param rOriginVariable 
     * @param rDestinationVariable 
     * @param BufferStep 
     */
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

    static bool GeneralizedInvertMatrix(
        Matrix& rInputMatrix,
        Matrix& rResult);

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
