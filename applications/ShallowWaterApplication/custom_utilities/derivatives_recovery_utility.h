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
 * @details Zhimin Zhangand Ahmed Naga (2006)
 * "A new finite element gradient recovery method: superconvergence property"
 * SIAM J. Sci. Comput., 26(4), 1192–1213. (22 pages)
 *
 * The polynomial follows the next convention:
 *    p_2(x, y, z, Node) = trans(P)*a
 *    trans(P) = (1, x, y, z, x^2, y^2, z^2, xy, xz, yz)
 *    trans(a) = a_i
 */
class DerivativesRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef Node<3> NodeType;

    static constexpr std::size_t TDim = 2;

    ///@}
    ///@name Pointer Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DerivativesRecoveryUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    // /// Default constructor.
    // DerivativesRecoveryUtility();

    // /// Destructor.
    // virtual ~DerivativesRecoveryUtility();

    ///@}
    ///@name Operations
    ///@{

    static void CalculateDivergence(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    static void CalculateGradient(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable,
        const Variable<array_1d<double,3>>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    static void CalculateLaplacian(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rOriginVariable,
        const Variable<array_1d<double,3>>& rDestinationVariable,
        const Variable<double>& rIntermediateVariable,
        const std::size_t BufferStep = 0);

    static void CalculateSuperconvergentDivergence(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    static void CalculateSuperconvergentGradient(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable,
        const Variable<array_1d<double,3>>& rDestinationVariable,
        const std::size_t BufferStep = 0);

    static void CalculateSuperconvergentLaplacian(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3>>& rOriginVariable,
        const Variable<array_1d<double,3>>& rDestinationVariable,
        const Variable<double>& rIntermediateVariable,
        const std::size_t BufferStep = 0);

    static void ExtendNeighborsPatch(ModelPart& rModelPart);

    static void CalculatePolynomialWeights(ModelPart& rModelPart);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
    // virtual std::string Info() const;

    // /// Print information about this object.
    // virtual void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    static double CalculateMaximumDistance(
        const NodeType& rNode,
        GlobalPointersVector<NodeType>& rNeighbors);

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
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


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    // /// Assignment operator.
    // DerivativesRecoveryUtility& operator=(DerivativesRecoveryUtility const& rOther);

    // /// Copy constructor.
    // DerivativesRecoveryUtility(DerivativesRecoveryUtility const& rOther);


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
