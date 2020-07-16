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

#if !defined(KRATOS_K_EPSILON_ADJOINT_ELEMENT_DATA_EPSILON_ADJOINT_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_K_EPSILON_ADJOINT_ELEMENT_DATA_EPSILON_ADJOINT_ELEMENT_DATA_H_INCLUDED

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
#include "custom_elements/convection_diffusion_reaction_adjoint_element_data_extension.h"
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/epsilon_element_data.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{
namespace KEpsilonElementData
{
template <unsigned int TDim, unsigned int TNumNodes>
class EpsilonAdjointElementData
    : public ConvectionDiffusionReactionAdjointElementDataExtension<TDim, TNumNodes>,
      public EpsilonElementData<TDim>
{
public:
    using NodeType = Node<3>;
    using GeometryType = Geometry<NodeType>;
    using PrimalBaseType = EpsilonElementData<TDim>;
    using AdjointBaseType =
        ConvectionDiffusionReactionAdjointElementDataExtension<TDim, TNumNodes>;

    EpsilonAdjointElementData(const GeometryType& rGeometry)
        : PrimalBaseType(rGeometry)
    {
    }

    static const Variable<double>& GetAdjointScalarVariable()
    {
        return RANS_SCALAR_2_ADJOINT_1;
    }
    static const Variable<double>& GetAdjointScalarRateVariable()
    {
        return RANS_SCALAR_2_ADJOINT_3;
    }

    static const std::string GetName()
    {
        return "KEpsilonEpsilonAdjointElementData";
    }

    static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGaussPointData(const Vector& rShapeFunctions,
                                 const Matrix& rShapeFunctionDerivatives,
                                 const int Step = 0) override;

    void CalculateEffectiveVelocityDerivatives(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                               const Variable<double>& rDerivativeVariable,
                                               const Vector& rShapeFunctions,
                                               const Matrix& rShapeFunctionDerivatives) const override;

    void CalculateEffectiveVelocityDerivatives(
        BoundedMatrix<double, TNumNodes * TDim, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const override;

    array_1d<double, 3> CalculateEffectiveVelocityShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

    void CalculateEffectiveKinematicViscosityDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const override;

    void CalculateEffectiveKinematicViscosityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rDerivativeVariable,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const override;

    double CalculateEffectiveKinematicViscosityShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

    void CalculateReactionTermDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                          const Variable<double>& rDerivativeVariable,
                                          const Vector& rShapeFunctions,
                                          const Matrix& rShapeFunctionDerivatives) const override;

    void CalculateReactionTermDerivatives(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                          const Variable<array_1d<double, 3>>& rDerivativeVariable,
                                          const Vector& rShapeFunctions,
                                          const Matrix& rShapeFunctionDerivatives) const override;

    double CalculateReactionTermShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

    void CalculateSourceTermDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                        const Variable<double>& rDerivativeVariable,
                                        const Vector& rShapeFunctions,
                                        const Matrix& rShapeFunctionDerivatives) const override;

    void CalculateSourceTermDerivatives(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                        const Variable<array_1d<double, 3>>& rDerivativeVariable,
                                        const Vector& rShapeFunctions,
                                        const Matrix& rShapeFunctionDerivatives) const override;

    double CalculateSourceTermShapeDerivatives(
        const ShapeParameter& rShapeParameters,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

private:
    BoundedMatrix<double, TNumNodes, TDim> mNodalVelocity;
    BoundedVector<double, TNumNodes> mGaussTurbulentKinematicViscositySensitivitiesK;
    BoundedVector<double, TNumNodes> mGaussTurbulentKinematicViscositySensitivitiesEpsilon;
};
} // namespace KEpsilonElementData
///@}
} // namespace Kratos

#endif