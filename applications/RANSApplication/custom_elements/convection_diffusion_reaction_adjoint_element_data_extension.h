//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT_DATA_EXTENSION_H_INCLUDED)
#define KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT_DATA_EXTENSION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/geometrical_sensitivity_utility.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class ConvectionDiffusionReactionAdjointElementDataExtension
{
public:
    /**
     * @brief Calculates effective velocity scalar derivatives
     *
     * This method is used to get scalar derivatives of effective velocity at gauss points in scalar transport equation.
     * \[
     *      M^{ci} = \frac{\partial u_i}{\partial \phi^c}
     * \]
     * Where $i$ is the effective velocity direction, $\phi^c$ is the scalar variables nodal index and $M^{ia}$ is the
     * rOutput matrix's $c^{th}$ row and $i^{th}$ column.
     *
     * @param rOutput                       Output matrix containing scalar derivatives. Rows corresponds to scalar derivative nodal index, columns corresponds to direction of effective velocity
     * @param rDerivativeVariable           Scalar derivative variable
     * @param rShapeFunctions               Shape function values at gauss point
     * @param rShapeFunctionDerivatives     Shape function derivatives at gauss point
     */
    virtual void CalculateEffectiveVelocityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates effective velocity vector derivatives
     *
     * This method is used to get vector derivatives of effective velocity at gauss points in scalar transport equation.
     * \[
     *      v = \frac{\partial u_i}{\partial \phi^c_j}
     * \]
     * Where $i$ is the effective velocity direction, $\phi^c_j$ is the vector variables $c^{th}$ node's $j^{th}$ direction
     * and v is the returned output. $c$ and $j$ information is passed through rDerivativeParameters variable.
     *
     * @param rDerivativeVariable           Vector derivative variable
     * @param rDerivativeParameters         Derivative parameters containing derivative nodal index and derivative direction
     * @param rShapeFunctions               Shape function values at gauss point
     * @param rShapeFunctionDerivatives     Shape function derivatives at gauss point
     * @return double                       The derivative of effective velocity
     */
    virtual array_1d<double, 3> CalculateEffectiveVelocityDerivatives(
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const ShapeParameter& rDerivativeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates shape derivative of effective velocity
     *
     * This method is used to get shape derivative of effective velocity at gauss points in scalar transport equation
     * \[
     *      v = \frac{\partial u_i}{\partial x^c_j}
     * \]
     * Where $i$ is the effective velocity direction, $x^c_j$ is the shape locations $c^{th}$ node's $j^{th}$ direction
     * and v is the returned output. $c$ and $j$ information is passed through rShapeParameters variable.
     *
     * @param rShapeParameters              Derivative shape parameters containing derivative nodal index and derivative direction
     * @param rShapeFunctions               Shape function values at gauss points
     * @param rShapeFunctionDerivatives     Shape function derivatives at gauss points
     * @param detJ_deriv                    Jacobian determinant shape derivative
     * @param rDN_Dx_deriv                  Shape function derivatives' shape derivative
     * @return double                       The shape derivative of effective velocity
     */
    virtual array_1d<double, 3> CalculateEffectiveVelocityShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

    /**
     * @brief Calculates effective kinematic viscosity scalar derivatives
     *
     * This method is used to calculate scalar derivative of effective kinematic viscosity at gauss points in scalar transport equations
     * \[
     *      V^c = \frac{\partial \nu}{\phi^c}
     * \]
     * Where $\nu$ is the effective kinematic viscosity, $\phi^c$ is the scalar variables nodal index and $V^c$ is the
     * rOutput vector's $c^{th}$ row.
     *
     * @param rOutput                       Output matrix containing $\nu$ derivatives for each nodal scalar derivative
     * @param rDerivativeVariable           Scalar derivative variable
     * @param rShapeFunctions               Shape function values at gauss point
     * @param rShapeFunctionDerivatives     Shape function derivatives at gauss point
     */
    virtual void CalculateEffectiveKinematicViscosityDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates effective kinematic viscosity vector derivatives
     *
     * This method is used to calculate vector derivative of effective kinematic viscosity at gauss points in scalar transport equations
     * \[
     *      M^{ci} = \frac{\partial \nu}{\phi^c_i}
     * \]
     * Where $\nu$ is the effective kinematic viscosity, $\phi^c_i$ is the scalar variables nodal index and direction and $M^{ci}$ is the
     * rOutput matrix's $c^{th}$ row and $i^{th}$ column.
     *
     * @param rOutput                       Output matrix containing $\nu$ derivatives for each nodal scalar derivative and its directions
     * @param rDerivativeVariable           Vector derivative variable
     * @param rShapeFunctions               Shape function values at gauss point
     * @param rShapeFunctionDerivatives     Shape function derivatives at gauss point
     */
    virtual void CalculateEffectiveKinematicViscosityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    virtual double CalculateEffectiveKinematicViscosityShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

    virtual void CalculateReactionTermDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                                  const Variable<double>& rDerivativeVariable,
                                                  const Vector& rShapeFunctions,
                                                  const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates reaction term vector derivatives
     *
     * This method is used to calculate vector derivative of reaction term at gauss points in scalar transport equations
     * \[
     *      M^{ci} = \frac{\partial s}{\phi^c_i}
     * \]
     * Where $s$ is the reaction term, $\phi^c_i$ is the scalar variables nodal index and direction and $M^{ci}$ is the
     * rOutput matrix's $c^{th}$ row and $i^{th}$ column.
     *
     * @param rOutput                       Output matrix containing reaction term derivatives for each nodal scalar derivative and its directions
     * @param rDerivativeVariable           Vector derivative variable
     * @param rShapeFunctions               Shape function values at gauss point
     * @param rShapeFunctionDerivatives     Shape function derivatives at gauss point
     */
    virtual void CalculateReactionTermDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    virtual double CalculateReactionTermShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

    virtual void CalculateSourceTermDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                                const Variable<double>& rDerivativeVariable,
                                                const Vector& rShapeFunctions,
                                                const Matrix& rShapeFunctionDerivatives) const = 0;

    /**
     * @brief Calculates source term vector derivatives
     *
     * This method is used to calculate vector derivative of source term at gauss points in scalar transport equations
     * \[
     *      M^{ci} = \frac{\partial f}{\phi^c_i}
     * \]
     * Where $f$ is the source term, $\phi^c_i$ is the scalar variables nodal index and direction and $M^{ci}$ is the
     * rOutput matrix's $c^{th}$ row and $i^{th}$ column.
     *
     * @param rOutput                       Output matrix containing source term derivatives for each nodal scalar derivative and its directions
     * @param rDerivativeVariable           Vector derivative variable
     * @param rShapeFunctions               Shape function values at gauss point
     * @param rShapeFunctionDerivatives     Shape function derivatives at gauss point
     */
    virtual void CalculateSourceTermDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    virtual double CalculateSourceTermShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;
};

///@}
} // namespace Kratos

#endif