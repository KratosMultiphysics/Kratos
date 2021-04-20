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

#if !defined(KRATOS_DRAG_FREQUENCY_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_DRAG_FREQUENCY_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <functional>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"

// Application includes
#include "drag_response_function.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

template <unsigned int TDim>
class DragFrequencyResponseFunction : public DragResponseFunction<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = DragResponseFunction<TDim>;

    KRATOS_CLASS_POINTER_DEFINITION(DragFrequencyResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DragFrequencyResponseFunction(
        Parameters Settings,
        ModelPart& rModelPart);

    /// Destructor.
    ~DragFrequencyResponseFunction() override = default;

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


    void CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;


    void CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Element& rAdjointElement,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    double CalculateValue(ModelPart& rModelPart) override;

    ///@}

private:
    ///@name Private Member Variables
    ///@{

    using BaseType::mrModelPart;
    using BaseType::mStructureModelPartName;
    using BaseType::mDragDirection;

    int mEchoLevel;

    int mFrequencyBinIndex;
    int mTotalNumberOfTimeSteps;
    int mWindowSize;

    bool mIsRealComponentRequested;

    std::function<double(double)> mComponentFunction;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateDragFrequencyContribution(
        const Matrix& rDerivativesOfResidual,
        const Element::NodesArrayType& rNodes,
        const int Step,
        Vector& rDerivativesOfDrag) const;

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_DRAG_FREQUENCY_RESPONSE_FUNCTION_H_INCLUDED defined */
