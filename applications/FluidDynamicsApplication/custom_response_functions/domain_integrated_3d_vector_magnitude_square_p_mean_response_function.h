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

#if !defined(KRATOS_DOMAIN_INTEGRATED_MAGNITUDE_P_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_DOMAIN_INTEGRATED_MAGNITUDE_P_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
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
template<unsigned int TDim>
class DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using ElementType = ModelPart::ElementType;

    using ConditionType = ModelPart::ConditionType;

    using GeometryType = typename ConditionType::GeometryType;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction(
        Parameters Settings,
        ModelPart& rModelPart);

    /// Destructor.
    ~DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

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

    double CalculateValue(ModelPart&) override;

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mIntegrationDomainModelPartName;

    const Variable<array_1d<double, 3>> *mpVariable;
    const Flags *mpDomainFlag;

    bool mIsElementsConsidered = false;
    bool mIsConditionsConsidered = false;

    double mIntegrationDomainSize;
    double mDomainIntegratedSquareMean;
    double mStartTime;

    IndexType mPower;
    IndexType mDofPosition;

    std::vector<Matrix> mShapeFunctions;

    ///@}
    ///@name Protected Operations
    ///@{

    template<class TEntityType>
    void CalculateResponseGradientContribution(
        const TEntityType& rEntity,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo);

    template<class TEntityType>
    void CalculateResponsePartialSensitivity(
        const TEntityType& rEntity,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo);

    template<class TEntityType>
    double CalculateDomainSizeDerivative(
        const TEntityType& rEntity,
        const IndexType DerivativeNodeIndex,
        const IndexType DerivativeDirectionIndex) const;

    double CalculateIntegrationDomainSize() const;

    double CalculateGeometryValueContribution(
        const GeometryType& rGeometry,
        Matrix& rShapeFunctions) const;

    void CalculateGeometryData(
        const GeometryType& rGeometry,
        Matrix& rNContainer) const;

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_DOMAIN_INTEGRATED_MAGNITUDE_P_RESPONSE_FUNCTION_H_INCLUDED defined */
