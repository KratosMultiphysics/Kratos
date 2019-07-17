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

/** \brief DirectSensitivityLocalStressResponseFunction
*
* This is the response class for local stresses used in the direct sensitivity analysis.
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
    DirectSensitivityLocalStressResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings, std::string& ResponseVariableName);

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
    
    void CalculateGradient(Element& rDirectElement,                            
                            Variable<Matrix> const& rStressVariable,
                            std::vector<std::vector<Matrix>>& rOutput, 
                            const ProcessInfo& rProcessInfo) override;

       
    
    void CalculatePartialSensitivity(Element& rDirectElement, 
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<array_1d<double, 3>> const& rStressVariable, 
                            std::vector<array_1d<double, 3>>& rOutput, 
                            const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<Matrix> const& rStressVariable, 
                                    std::vector<Matrix>& rOutput, 
                                    const ProcessInfo& rProcessInfo) override;

    std::string GetEvaluationFlag() override;

    

    


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

    template <typename TDataType>
    void PartialSensitivityBuilder(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<TDataType> const& rStressVariable, 
                                    std::vector<TDataType>& rOutput, 
                                    const ProcessInfo& rProcessInfo);
    
    template <typename TDataType>
    void CalculateElementContributionToPartialSensitivity(Element& rDirectElement,
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<TDataType> const& rStressVariable,
                                    std::vector<std::vector<TDataType>>& rOutput,
                                    const ProcessInfo& rProcessInfo);

    template <typename TDataType>
    void CalculateElementContributionToGradient(Element& rDirectElement,
                                    Variable<TDataType> const& rStressVariable,
                                    std::vector<std::vector<TDataType>>& rOutput,
                                    const ProcessInfo& rProcessInfo);

    template <typename TDataType>
    void ExtractMeanStressDerivative(const std::vector<TDataType>& rStressDerivativesMatrix,
                            std::vector<TDataType>& rOutput);

    template <typename TDataType>
    void ExtractGaussPointStressDerivative(const std::vector<TDataType>& rStressDerivativesMatrix,
                            std::vector<TDataType>& rOutput);

    
    ///@}
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* DIRECT_SENSITIVITY_LOCAL_STRESS_RESPONSE_FUNCTION_H defined */
