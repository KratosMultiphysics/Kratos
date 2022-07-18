// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kateryna Pindak
//

#ifndef ADJOINT_EIGENFREQUENCY_RESPONSE_FUNCTION_H
#define ADJOINT_EIGENFREQUENCY_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "adjoint_structural_response_function.h"

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

/** \brief
*
* It is designed to be used in adjoint sensitivity analysis.
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointEigenfrequencyResponseFunction : public AdjointStructuralResponseFunction
{

public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointEigenfrequencyResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointEigenfrequencyResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointEigenfrequencyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    ~AdjointEigenfrequencyResponseFunction() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    double CalculateValue(ModelPart& rPrimalModelPart) override;

    using AdjointStructuralResponseFunction::CalculateGradient;

    void CalculateGradient(const Element& rAdjointElement,
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
    double mWeightedSumOfEigenfrequencies;
    double mWeightedSumOfEigenvalues;

    ModelPart& mrAdjointModelPart;
    Element::Pointer mpElementInAdjointPart = nullptr;
    double mDelta;
    std::vector<int> mTracedEigenfrequencyIds;
    std::vector<double> mWeightingFactors;
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                    const std::string& rVariableName,
                                    const Matrix& rSensitivityMatrix,
                                    Vector& rSensitivityGradient,
                                    const ProcessInfo& rProcessInfo);

    void GetWeightedSumOfEigenfrequencies(const ProcessInfo& rProcessInfo);

    void GetWeightedSumOfEigenvalues(const ProcessInfo& rProcessInfo);

    void CheckSettingsForGradientAnalysis(Parameters rResponseSettings);

    void DetermineTracedEigenfrequencies(Parameters rResponseSettings);
    
    bool AreSeveralEigenfrequenciesTraced();

    void CalculateLinearWeights(Parameters rResponseSettings);

    void UseDefaultWeight();
    
    void CheckIfAllNecessaryEigenvaluesAreComputed(const ProcessInfo& rProcessInfo);

    double GetEigenvalue(const int eigenfrequency_id, const ProcessInfo& rProcessInfo);

    void GetWeightedSumOfEigenvectors(Element& rAdjointElement, Vector& rWeightedSumOfEigenvectors, const ProcessInfo& rProcessInfo);

    void DetermineEigenvectorOfElement(ModelPart::ElementType& rElement,
                                        const int eigenfrequency_id,
                                        Vector& rEigenvectorOfElement,
                                        const ProcessInfo& CurrentProcessInfo);

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
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_EIGENFREQUENCY_RESPONSE_FUNCTION_H_INCLUDED