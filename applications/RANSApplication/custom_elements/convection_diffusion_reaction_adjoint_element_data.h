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

#if !defined(KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief Adjoint data holder for convection diffusion reaction element
 *
 * This class holds interfaces to compute derivatives of the convection diffusion reaction
 * equation coefficients. The equation is illustrated below.
 *
 * \[
 *      \frac{\partial \phi}{\partial t} + u_i\frac{\partial \phi}{\partial x_i} + \nu\frac{\partial^2 \phi}{\partial x_i^2} + s\phi = f
 * \]
 *
 * @tparam TDerivativeVariableDim       Dimensionality of the derivative variable
 * @tparam TNumNodes                    Number of nodes in the geometry
 * @tparam TDataContainerType           Data holder
 */
template <unsigned int TDerivativeVariableDim, unsigned int TNumNodes, class TDataContainerType>
class ConvectionDiffusionReactionAdjointElementData
{
public:
    ///@name Public type definitions
    ///@{

    static constexpr unsigned int TDerivativesSize = TDerivativeVariableDim * TNumNodes;

    ///@}
    ///@name Life cycle
    ///@{

    ConvectionDiffusionReactionAdjointElementData(
        const TDataContainerType& rElementData)
        : mrElementData(rElementData)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    const TDataContainerType& GetElementData() const
    {
        return mrElementData;
    }

    ///@}
    ///@name Public abstract interface
    ///@{

    /**
     * @brief Calculates effective velocity derivatives
     *
     * \[
     *      M_a^i = \frac{\partial u_i}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     *
     * @param rOutput                       Output derivatives
     * @param rShapeFunctions               Shape functions
     * @param rShapeFunctionDerivatives     Shape function derivatives
     */
    virtual void CalculateEffectiveVelocityDerivatives(
        BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates effective viscosity derivatives
     *
     * \[
     *      V_a = \frac{\partial \nu}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     *
     * @param rOutput                       Output derivatives
     * @param rShapeFunctions               Shape functions
     * @param rShapeFunctionDerivatives     Shape function derivatives
     */
    virtual void CalculateEffectiveKinematicViscosityDerivatives(
        BoundedVector<double, TDerivativesSize>& rOutput,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates reaction term derivatives
     *
     * \[
     *      V_a = \frac{\partial s}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     *
     * @param rOutput                       Output derivatives
     * @param rShapeFunctions               Shape functions
     * @param rShapeFunctionDerivatives     Shape function derivatives
     */
    virtual void CalculateReactionTermDerivatives(
        BoundedVector<double, TDerivativesSize>& rOutput,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates source term derivatives
     *
     * \[
     *      V_a = \frac{\partial f}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     *
     * @param rOutput                       Output derivatives
     * @param rShapeFunctions               Shape functions
     * @param rShapeFunctionDerivatives     Shape function derivatives
     */
    virtual void CalculateSourceTermDerivatives(
        BoundedVector<double, TDerivativesSize>& rOutput,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    ///@}

private:
    ///@name Private members
    ///@{

    const TDataContainerType& mrElementData;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT_DATA_H_INCLUDED