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
    class ScalarDerivative
    {
    public:
        virtual void CalculateEffectiveVelocityDerivatives(
            BoundedMatrix<double, TNumNodes, TDim>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const = 0;

        virtual void CalculateEffectiveKinematicViscosityDerivatives(
            BoundedVector<double, TNumNodes>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const = 0;

        virtual void CalculateReactionTermDerivatives(
            BoundedVector<double, TNumNodes>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const = 0;

        virtual void CalculateSourceTermDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                                    const Vector& rShapeFunctions,
                                                    const Matrix& rShapeFunctionDerivatives) const = 0;
    };

    class VectorDerivative
    {
    public:
        virtual void CalculateEffectiveVelocityDerivatives(
            BoundedMatrix<double, TNumNodes * TDim, TDim>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const = 0;

        virtual void CalculateEffectiveKinematicViscosityDerivatives(
            BoundedMatrix<double, TNumNodes, TDim>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const = 0;

        virtual void CalculateReactionTermDerivatives(
            BoundedMatrix<double, TNumNodes, TDim>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const = 0;

        virtual void CalculateSourceTermDerivatives(
            BoundedMatrix<double, TNumNodes, TDim>& rOutput,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives) const = 0;
    };

    class ShapeDerivative
    {
    public:
        virtual array_1d<double, 3> CalculateEffectiveVelocityShapeDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

        virtual double CalculateEffectiveKinematicViscosityShapeDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

        virtual double CalculateReactionTermShapeDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;

        virtual double CalculateSourceTermShapeDerivatives(
            const ShapeParameter& rShapeParameters,
            const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives,
            const double detJ_deriv,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const = 0;
    };
};

///@}
} // namespace Kratos

#endif