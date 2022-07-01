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

#if !defined(KRATOS_GAUSS_POINT_KREISSELMEIER_AGGREGATION_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_GAUSS_POINT_KREISSELMEIER_AGGREGATION_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

class GaussPointKreisselmeierAggregationResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(GaussPointKreisselmeierAggregationResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    GaussPointKreisselmeierAggregationResponseFunction(
        Parameters Settings,
        ModelPart& rModelPart);

    /// Destructor.
    ~GaussPointKreisselmeierAggregationResponseFunction() override = default;
    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    double CalculateValue(ModelPart& rModelPart) override;

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

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mCriticalModelPartName;
    std::string mPerturbationVariableName;
    double mGaussPointValueScalingFactor;
    double mStepSize;

    const Variable<double>* mpGaussPointValueScalarVariable;
    std::vector<const Variable<double>*> mPrimalStateScalarVariablePointersList;
    std::vector<const Variable<double>*> mPrimalStateFirstDerivativeScalarVariablePointersList;
    std::vector<const Variable<double>*> mPrimalStateSecondDerivativeScalarVariablePointersList;

    std::map<int,double> mKSPrefactors;
    double mSumKSPrefactors;
    double mRho;

    int mEchoLevel;

    bool mAreKSPrefactorsInitialized;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateFiniteDifferenceStateVariableSensitivities(
        Vector& rOutput,
        ModelPart::ElementType& rElement,
        const std::vector<const Variable<double>*>& rDerivativeVariablePointersList,
        const ProcessInfo& rProcessInfo);

    void CalculateFiniteDifferenceShapeVariableSensitivities(
        Vector& rOutput,
        ModelPart::ElementType& rElement,
        const ProcessInfo& rProcessInfo);

    void CalculateStateDerivative(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const std::vector<const Variable<double>*>& rDerivativeVariablesList,
        const ProcessInfo& rProcessInfo);

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_GAUSS_POINT_KREISSELMEIER_AGGREGATION_RESPONSE_FUNCTION_H_INCLUDED defined */
