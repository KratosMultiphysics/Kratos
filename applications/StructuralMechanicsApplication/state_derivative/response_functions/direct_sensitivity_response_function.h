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

#ifndef DIRECT_SENSITIVITY_RESPONSE_FUNCTION_H
#define DIRECT_SENSITIVITY_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "state_derivative/variable_utilities/direct_sensitivity_variable.h"
#include "derivative_builder.h"
#include "state_derivative/math_functions/vector_math.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/** \brief DirectSensitivityResponseFunction
*
* This is the response base class for responses in structural mechanics.
* It is designed to be used in direct sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DirectSensitivityResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DirectSensitivityResponseFunction);

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    typedef Variable<array_1d<double, 3>> VariableWithComponentsType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef Vector VectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Pointer Definitions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DirectSensitivityResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings, std::string& ResponseVariableName);


    /// Destructor.
    ~DirectSensitivityResponseFunction() {}

    ///@}
    ///@name Operators
    ///@{
    
    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize();

    virtual void CalculateGradient(Element& rDirectElement,                            
                            Variable<array_1d<double, 3>> const& rStressVariable,
                            std::vector<std::vector<array_1d<double, 3>>>& rOutput, 
                            const ProcessInfo& rProcessInfo);

    virtual void CalculateGradient(Element& rDirectElement,                            
                            Variable<Matrix> const& rStressVariable,
                            std::vector<std::vector<Matrix>>& rOutput, 
                            const ProcessInfo& rProcessInfo);
    
    virtual void CalculateGradient(Node<3>& rNode,                            
                            Variable<array_1d<double, 3>> const& rStressVariable,
                            std::vector<array_1d<double, 3>>& rOutput, 
                            const ProcessInfo& rProcessInfo); 

    virtual void CalculateGradient(Node<3>& rNode,                            
                            Variable<Matrix> const& rStressVariable,
                            std::vector<Matrix>& rOutput, 
                            const ProcessInfo& rProcessInfo);   

    virtual void CalculatePartialSensitivity(Element& rDirectElement, 
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<array_1d<double, 3>> const& rStressVariable, 
                            std::vector<array_1d<double, 3>>& rOutput, 
                            const ProcessInfo& rProcessInfo);
    
    virtual void CalculatePartialSensitivity(Element& rDirectElement, 
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<Matrix> const& rStressVariable, 
                            std::vector<Matrix>& rOutput, 
                            const ProcessInfo& rProcessInfo);

    virtual void CalculatePartialSensitivity(Node<3>& rNode, 
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<array_1d<double, 3>> const& rStressVariable, 
                            array_1d<double, 3>& rOutput, 
                            const ProcessInfo& rProcessInfo);
    
    virtual void CalculatePartialSensitivity(Node<3>& rNode, 
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<Matrix> const& rStressVariable, 
                            Matrix& rOutput, 
                            const ProcessInfo& rProcessInfo);
    
    virtual std::string GetResponseVariableName();

    virtual std::string GetEvaluationFlag(); 

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;
    unsigned int mGradientMode;
    std::string mResponseVariableName;
        
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* DIRECT_SENSITIVITY_RESPONSE_FUNCTION_H defined */