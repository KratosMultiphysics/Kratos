//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_RESIDUAL_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_RESIDUAL_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// A response function for drag.
/**
 * The response function is defined as:
 *
 * \f[
 * \bar{D} = \Sigma_{n=1}^N D^n \Delta t
 * \f]
 *
 * if "integrate_in_time" is true.
 */
template <unsigned int TDim>
class ResidualResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using IndexType = std::size_t;

    using ArrayD = array_1d<double, TDim>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    KRATOS_CLASS_POINTER_DEFINITION(ResidualResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ResidualResponseFunction(
        Parameters Settings,
        ModelPart& rModelPart);

    /// Destructor.
    ~ResidualResponseFunction() override = default;

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep() override;

    void CalculateGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Element& rAdjointElement,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    double CalculateValue(ModelPart& rModelPart) override;

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    ModelPart& mrModelPart;
    double mMomentumResidualWeight;
    double mContinuityResidualWeight;
    bool mIsElementErrorSaved;

    ///@}

private:
    ///@name Private Operations
    ///@{

    void CalculateStabilizationParameters(
        double& TauU,
        double& TauP,
        const double VelocityMagnitude,
        const double ElementSize,
        const double KinematicViscosity,
        const double TauDynamicMultiplier,
        const ProcessInfo& rCurrentProcessInfo) const;

    void CalculateStabilizationParameterDerivatives(
        double& TauUDerivative,
        double& TauPDerivative,
        const double TauU,
        const double TauP,
        const double VelocityMagnitude,
        const double VelocityMagnitudeDerivative,
        const double ElementSize,
        const double ElementSizeDerivative,
        const double KinematicViscosity,
        const ProcessInfo& rCurrentProcessInfo) const;

    void CalculateGeometryData(
        const Element& rElement,
        Vector& rGaussWeights,
        Matrix& rNContainer,
        ShapeFunctionDerivativesArrayType& rDN_DX) const;

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_RESPONSE_FUNCTION_H_INCLUDED defined */
