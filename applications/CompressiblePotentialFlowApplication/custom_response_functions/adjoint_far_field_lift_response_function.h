//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nunez, Inigo Lopez
//

#ifndef ADJOINT_LIFT_FAR_FIELD_RESPONSE_FUNCTION_H
#define ADJOINT_LIFT_FAR_FIELD_RESPONSE_FUNCTION_H


// System includes

// External includes

// Project includes
#include "adjoint_potential_response_function.h"


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

/** \brief AdjointLiftFarFieldResponseFunction
*
* This is a response function that computes the lift sensitivity using the far field integral
*/
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) AdjointLiftFarFieldResponseFunction : public AdjointPotentialResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointLiftFarFieldResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointLiftFarFieldResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointLiftFarFieldResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    ~AdjointLiftFarFieldResponseFunction() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    void InitializeSolutionStep() override;

    double CalculateValue(ModelPart& rModelPart) override;

    void CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo) override;

    void CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mReferenceChord;
    double mDynamicPressure;
    double mStepSize;
    std::string mFarFieldModelPartName;
    array_1d<double, 3> mFreeStreamVelocity;
    array_1d<double, 3> mWakeNormal;


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    double ComputeLiftContribution(Element& rElement, const array_1d<double, 3> rNormal, const ProcessInfo& rProcessInfo);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class AdjointLiftFarFieldResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_LIFT_FAR_FIELD_RESPONSE_FUNCTION_H_INCLUDED
