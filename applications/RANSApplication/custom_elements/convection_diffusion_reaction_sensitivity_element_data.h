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

#if !defined(KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_SENSITIVITY_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_SENSITIVITY_ELEMENT_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "utilities/geometrical_sensitivity_utility.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief Adjoint data holder for convection diffusion reaction element
 *
 * This class holds interfaces to compute shape derivatives of the convection diffusion reaction
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
class ConvectionDiffusionReactionSensitivityElementData
{
public:
    ///@name Life cycle
    ///@{

    ConvectionDiffusionReactionSensitivityElementData(
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
     * @brief Calculate effective velocity shape derivatives
     *
     * \[
     *      output_i = \frac{\partial u_i}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     * $c$ and $k$ is stored inside rShapeParameters.
     *
     * @param rShapeParameters
     * @param rShapeFunctions
     * @param rShapeFunctionDerivatives
     * @param detJ_deriv
     * @param rDN_Dx_deriv
     * @return array_1d<double, 3>
     */
    virtual array_1d<double, 3> CalculateEffectiveVelocityDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

    /**
     * @brief Calculate effective viscosity shape derivatives
     *
     * \[
     *      output = \frac{\partial \nu}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     * $c$ and $k$ is stored inside rShapeParameters.
     *
     * @param rShapeParameters
     * @param rShapeFunctions
     * @param rShapeFunctionDerivatives
     * @param detJ_deriv
     * @param rDN_Dx_deriv
     * @return array_1d<double, 3>
     */
    virtual double CalculateEffectiveKinematicViscosityDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

    /**
     * @brief Calculate reaction term shape derivatives
     *
     * \[
     *      output = \frac{\partial s}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     * $c$ and $k$ is stored inside rShapeParameters.
     *
     * @param rShapeParameters
     * @param rShapeFunctions
     * @param rShapeFunctionDerivatives
     * @param detJ_deriv
     * @param rDN_Dx_deriv
     * @return array_1d<double, 3>
     */
    virtual double CalculateReactionTermDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

    /**
     * @brief Calculate source term shape derivatives
     *
     * \[
     *      output = \frac{\partial f}{\partial \omega_a}
     * \]
     *
     * In here $u_i$ derivative is calculated w.r.t. $\omega_a$.
     * If $\omega$ is a not a scalar then $a = c * TDim + k$, where
     * $c$ is the node index, $k$ is the direction index
     *
     * $c$ and $k$ is stored inside rShapeParameters.
     *
     * @param rShapeParameters
     * @param rShapeFunctions
     * @param rShapeFunctionDerivatives
     * @param detJ_deriv
     * @param rDN_Dx_deriv
     * @return array_1d<double, 3>
     */
    virtual double CalculateSourceTermDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

    ///@}

private:
    ///@name Private members
    ///@{

    const TDataContainerType& mrElementData;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_SENSITIVITY_ELEMENT_DATA_H_INCLUDED