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

#if !defined(KRATOS_FREQUENCY_BIN_COMPONENT_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_FREQUENCY_BIN_COMPONENT_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"

// Application includes

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

template<unsigned int TDim>
class WindowedFrequencyBinComponentResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(WindowedFrequencyBinComponentResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    WindowedFrequencyBinComponentResponseFunction(Parameters Settings, ModelPart& rModelPart);

    /// Destructor.
    ~WindowedFrequencyBinComponentResponseFunction()  = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

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

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    int mEchoLevel;

    Vector mPointShapeFunctionValues;
    IndexType mPointElementId = 0;
    array_1d<double, 3> mPointCoordinates;
    array_1d<double, 3> mVelocityDirection;

    int mFrequencyBinIndex;
    int mTotalNumberOfTimeSteps;
    int mWindowSize;

    bool mIsRealComponentRequested;

    ///@}
    ///@name Private Operations
    ///@{

    void FindPointElement();

    double CalculateHannWindowValue();

    ///@}
};

///@} // Kratos Classes

///@} // RANSApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_FREQUENCY_BIN_COMPONENT_RESPONSE_FUNCTION_H_INCLUDED defined */
