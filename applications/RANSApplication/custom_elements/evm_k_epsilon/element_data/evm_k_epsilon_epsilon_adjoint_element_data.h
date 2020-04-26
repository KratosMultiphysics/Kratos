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

#if !defined(KRATOS_EVM_K_EPSILON_ELEMENT_DATA_UTILITIES_EPSILON_ADJOINT_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_EVM_K_EPSILON_ELEMENT_DATA_UTILITIES_EPSILON_ADJOINT_ELEMENT_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_adjoint_element_data_extension.h"
#include "custom_elements/evm_k_epsilon/element_data/evm_k_epsilon_epsilon_element_data.h"

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

namespace EvmKEpsilonAdjointElementDataUtilities
{
template <unsigned int TDim, unsigned int TNumNodes>
class EpsilonAdjointElementData
    : public EvmKEpsilonElementDataUtilities::EpsilonElementData<TDim>,
      public ScalarConvectionDiffusionReactionAdjointElementDataExtension<TDim, TNumNodes>
{
public:
    using PrimalBaseType = EvmKEpsilonElementDataUtilities::EpsilonElementData<TDim>;
    using AdjointBaseType =
        ScalarConvectionDiffusionReactionAdjointElementDataExtension<TDim, TNumNodes>;
    using NodeType = Node<3>;
    using GeometryType = typename PrimalBaseType::GeometryType;

    static const Variable<double>& GetAdjointScalarVariable();
    static const Variable<double>& GetAdjointScalarRateVariable();

    static void Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

    EpsilonAdjointElementData(const GeometryType& rGeometry)
        : PrimalBaseType(rGeometry)
    {
    }

    void CalculateGaussPointData(const Vector& rShapeFunctions,
                                 const Matrix& rShapeFunctionDerivatives,
                                 const int Step = 0) override;

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

    void CalculateReactionTermDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                          const Variable<double>& rDerivativeVariable,
                                          const Vector& rShapeFunctions,
                                          const Matrix& rShapeFunctionDerivatives) const override;

    void CalculateReactionTermDerivatives(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                          const Variable<array_1d<double, 3>>& rDerivativeVariable,
                                          const Vector& rShapeFunctions,
                                          const Matrix& rShapeFunctionDerivatives) const override;

    void CalculateSourceTermDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                        const Variable<double>& rDerivativeVariable,
                                        const Vector& rShapeFunctions,
                                        const Matrix& rShapeFunctionDerivatives) const override;

    void CalculateSourceTermDerivatives(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                        const Variable<array_1d<double, 3>>& rDerivativeVariable,
                                        const Vector& rShapeFunctions,
                                        const Matrix& rShapeFunctionDerivatives) const override;

    double CalculateEffectiveKinematicViscosityShapeSensitivity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

    double CalculateReactionTermShapeSensitivity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

    double CalculateSourceTermShapeSensitivity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const override;

protected:
    BoundedMatrix<double, TNumNodes, TDim> mNodalVelocity;
    BoundedVector<double, TNumNodes> mGaussTurbulentKinematicViscositySensitivitiesK;
    BoundedVector<double, TNumNodes> mGaussTurbulentKinematicViscositySensitivitiesEpsilon;
};

} // namespace EvmKEpsilonAdjointElementDataUtilities
} // namespace Kratos

#endif