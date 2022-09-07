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

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) GaussPointKreisselmeierAggregationResponseFunction : public AdjointResponseFunction
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

    enum GradientMode
    {
        SEMI_ANALITIC,
        ANALYTIC
    };

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mCriticalModelPartName;
    double mAggregationPenalty;
    int mEchoLevel;

    GradientMode mGradientMode;

    // semi-analytic gradient mode settings
    const Variable<double>* mpPerturbationVariable = nullptr;
    double mPerturbationSize = 0.0;
    const Variable<std::string>* mpDeisgnVariableNameStorageVariable;

    // scaling factor for gauss point evaluations
    double mGaussPointValueScalingFactor = 0.0;

    // gauss point variables
    const Variable<double>* mpGaussPointValueScalarVariable;
    const Variable<Matrix>* mpGaussPointValueGradientVariable = nullptr;
    const Variable<Matrix>* mpGaussPointValueFirstDerivativeVariable = nullptr;
    const Variable<Matrix>* mpGaussPointValueSecondDerivativeVariable = nullptr;
    const Variable<Matrix>* mpGaussPointValueShapeDerivativeVariable = nullptr;

    // variables required for gauss point value and derivative computations.
    std::map<int, double> mKSPrefactors;
    double mSumKSPrefactors;
    bool mAreKSPrefactorsInitialized;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TDataType>
    void SetGaussPointVariable(
        const Variable<TDataType>*& pVariable,
        const std::string& rVariableName,
        const std::string& rMsg);

    template<class TEntityType>
    void CalculateGaussPointDerivatives(
        Vector& rOutput,
        TEntityType& rEntity,
        const double Factor,
        const Variable<Matrix>* pVariable,
        const Matrix& rResidualGradient,
        const ProcessInfo& rProcessInfo) const;

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_GAUSS_POINT_KREISSELMEIER_AGGREGATION_RESPONSE_FUNCTION_H_INCLUDED defined */
