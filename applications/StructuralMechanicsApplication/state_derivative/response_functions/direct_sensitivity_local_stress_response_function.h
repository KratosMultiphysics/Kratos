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
#include "derivative_builder.h"

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

    typedef Element::DofsVectorType DofsVectorType;

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
                            Variable<array_1d<double, 3>> const& rStressVariable,
                            std::vector<std::vector<array_1d<double, 3>>>& rOutput, 
                            const ProcessInfo& rProcessInfo) override;

       
    
    void CalculatePartialSensitivity(Element& rDirectElement, 
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<array_1d<double, 3>> const& rStressVariable, 
                            std::vector<array_1d<double, 3>>& rOutput, 
                            const ProcessInfo& rProcessInfo) override;

    std::vector<std::string> GetResponseSensitivityVariableVector() override;

    


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
        std::vector<std::string> mTracedStressesVector;
        std::vector<unsigned int> mIdOfLocationVector;
              
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    void CalculateElementContributionToPartialSensitivity(Element& rDirectElement,
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<array_1d<double, 3>> const& rStressVariable,
                            std::vector<std::vector<array_1d<double, 3>>>& rOutput,
                            const ProcessInfo& rProcessInfo);

    void ExtractMeanStressDerivative(const std::vector<array_1d<double, 3>>& rStressDerivativesMatrix,
                            std::vector<array_1d<double, 3>>& rOutput);

    void ExtractGaussPointStressDerivative(const std::vector<array_1d<double, 3>>& rStressDerivativesMatrix,
                            std::vector<array_1d<double, 3>>& rOutput);

    
    ///@}
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* DIRECT_SENSITIVITY_LOCAL_STRESS_RESPONSE_FUNCTION_H defined */
