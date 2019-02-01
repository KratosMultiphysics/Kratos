// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kevin Braun, https://github.com/MFusseder
//

#ifndef DIRECT_SENSITIVITY_LOCAL_STRESS_RESPONSE_FUNCTION_H
#define DIRECT_SENSITIVITY_LOCAL_STRESS_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "direct_sensitivity_response_function.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/** \brief AdjointStructuralResponseFunction
*
* This is the response base class for responses in structural mechanics.
* It is designed to be used in adjoint sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DirectSensitivityLocalStressResponseFunction : public DirectSensitivityResponseFunction 
{
public:
    ///@name Type Definitions
    ///@{

    typedef DirectSensitivityResponseFunction BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(DirectSensitivityLocalStressResponseFunction);

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Vector VectorType;

    ///@}
    ///@name Pointer Definitions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DirectSensitivityLocalStressResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
     ~DirectSensitivityLocalStressResponseFunction();    

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    void CalculateGradient(Element& rDirectElement, 
                            const Matrix& rLHS,
                            Matrix& rResponseGradientMatrix,                                     
                            const ProcessInfo& rProcessInfo) override;

       
    
    void CalculatePartialSensitivity(Element& rDirectElement, 
                            DirectSensitivityVariable& DesignVariable,
                            Matrix& rSensitivityGradient,
                            const ProcessInfo& rProcessInfo) override;

    int GetNumberOfTracedGaussPoints() override;

    


    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{
        
        StressTreatment mStressTreatment;
        std::vector<TracedStressType> mTracedStressesVector;
        std::vector<unsigned int> mIdOfLocationVector;
              
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    void CalculateElementContributionToGradient(Element& rDirectElement,
                            const Matrix& rLHS,
                            Vector& rResponseGradientVector,
                            unsigned int& rIdOfGaussPoint,
                            ProcessInfo& rProcessInfo); 

    void CalculateElementContributionToPartialSensitivity(Element& rDirectElement,
                            const std::string& rVariableName,
                            Vector& rSensitivityGradient,
                            unsigned int& rIdOfGaussPoint,
                            ProcessInfo& rProcessInfo);

    void ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient);

    void ExtractGaussPointStressDerivative(const Matrix& rStressDerivativesMatrix, 
                            unsigned int& rIdOfGaussPoint,
                            Vector& rResponseGradient);
    ///@}
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* DIRECT_SENSITIVITY_LOCAL_STRESS_RESPONSE_FUNCTION_H defined */
