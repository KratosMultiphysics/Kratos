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

#if !defined(KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT_DATA_H_INCLUDED

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
#include "custom_elements/convection_diffusion_reaction_element_data.h"

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

template <unsigned int TDim, unsigned int TNumNodes>
class ScalarConvectionDiffusionReactionAdjointElementData
    : public ScalarConvectionDiffusionReactionElementData
{
public:
    using BaseType = ScalarConvectionDiffusionReactionElementData;
    using GeometryType = BaseType::GeometryType;

    ScalarConvectionDiffusionReactionAdjointElementData(const GeometryType& rGeometry)
        : BaseType(rGeometry)
    {
    }

    virtual double CalculateEffectiveKinematicViscosityDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual void CalculateEffectiveKinematicViscosityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual double CalculateReactionTermDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual void CalculateReactionTermDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual double CalculateSourceTermDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual void CalculateSourceTermDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual double CalculateEffectiveKinematicViscosityShapeSensitivity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual double CalculateReactionTermShapeSensitivity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual double CalculateSourceTermShapeSensitivity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const = 0;
};
} // namespace Kratos

#endif