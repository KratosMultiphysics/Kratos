// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner, modified by E. Wehrle for KS
//

#ifndef ADJOINT_KS_MAX_STRESS_RESPONSE_FUNCTION_H
#define ADJOINT_KS_MAX_STRESS_RESPONSE_FUNCTION_H


// System includes

// External includes

// Project includes
#include "adjoint_structural_response_function.h"
#include "stress_response_definitions.h"

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

/** \brief AdjointKSMaxStressResponseFunction
 *
 * This is a response function which approximates the max stress via
 * modified Kreisselmeier-Steinhauser function aggregation
 * It is designed to be used in adjoint sensitivity analysis.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointKSMaxStressResponseFunction : public AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointKSMaxStressResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointKSMaxStressResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointKSMaxStressResponseFunction(ModelPart& rModelPart,
                                       Parameters ResponseSettings);

    /// Destructor.
    ~AdjointKSMaxStressResponseFunction() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
                                     const Variable<array_1d<double,
                                     3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                     const Variable<array_1d<double,
                                     3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override;

    double CalculateValue(ModelPart& rModelPart) override;

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

    ModelPart& mrAdjointModelPart;
    std::string mCriticalPartName;
    Element::Pointer mpTracedElementInAdjointPart = nullptr;
    std::vector<int> mAggregatedElementIds;
    StressTreatment mStressTreatment;
    TracedStressType mTracedStressType;
    SizeType mEchoLevel = 0;
    double max_mean_stress = 0.0;
    double KS_exp_sum = 0.0;
    std::map<int,double> mean_stress_vector;
    //double pKS = 100.0;  //this should be a variable to set
    //pKS mpKS;
    double pKS;

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

    void ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix,
                                     Vector& rResponseGradient);

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

}; // Class AdjointKSMaxStressResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_KS_MAX_STRESS_RESPONSE_FUNCTION_H_INCLUDED
