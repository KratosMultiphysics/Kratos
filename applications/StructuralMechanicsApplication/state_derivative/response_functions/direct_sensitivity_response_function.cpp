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


// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "direct_sensitivity_response_function.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

    /// Constructor.
    DirectSensitivityResponseFunction::DirectSensitivityResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings, std::string& ResponseVariableName)
      : mrModelPart(rModelPart), mResponseVariableName(ResponseVariableName)
    {
        KRATOS_TRY;

        // Set gradient mode
        const std::string& gradient_mode = ResponseSettings["gradient_mode"].GetString();

        KRATOS_CATCH("");
    }

    void DirectSensitivityResponseFunction::Initialize()
    {
        KRATOS_TRY; 

        KRATOS_CATCH("");
    }

    
    void DirectSensitivityResponseFunction::CalculateGradient(Element& rDirectElement,                            
                                    Variable<array_1d<double, 3>> const& rResponseVariable,
                                    std::vector<std::vector<array_1d<double, 3>>>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculateGradient() needs to be implemented by the derived class.\n";  
    }

    void DirectSensitivityResponseFunction::CalculateGradient(Element& rDirectElement,                            
                            Variable<Matrix> const& rStressVariable,
                            std::vector<std::vector<Matrix>>& rOutput, 
                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculateGradient() needs to be implemented by the derived class.\n";  
    }

    void DirectSensitivityResponseFunction::CalculateGradient(Node<3>& rNode,                            
                                    Variable<array_1d<double, 3>> const& rStressVariable,
                                    std::vector<array_1d<double, 3>>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculateGradient() needs to be implemented by the derived class.\n";  
    } 

    void DirectSensitivityResponseFunction::CalculateGradient(Node<3>& rNode,                            
                                    Variable<Matrix> const& rStressVariable,
                                    std::vector<Matrix>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculateGradient() needs to be implemented by the derived class.\n";  
    }       
    
    void DirectSensitivityResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<array_1d<double, 3>> const& rResponseVariable, 
                                    std::vector<array_1d<double, 3>>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculatePartialSensitivity() needs to be implemented by the derived class.\n";  
    }

    void DirectSensitivityResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<Matrix> const& rResponseVariable, 
                                    std::vector<Matrix>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculatePartialSensitivity() needs to be implemented by the derived class.\n";  
    }

    void DirectSensitivityResponseFunction::CalculatePartialSensitivity(Node<3>& rNode, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<array_1d<double, 3>> const& rStressVariable, 
                                    array_1d<double, 3>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculatePartialSensitivity() needs to be implemented by the derived class.\n";  
    }

    void DirectSensitivityResponseFunction::CalculatePartialSensitivity(Node<3>& rNode, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<Matrix> const& rStressVariable, 
                                    Matrix& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculatePartialSensitivity() needs to be implemented by the derived class.\n";  
    }
    
        
    std::string DirectSensitivityResponseFunction::GetResponseVariableName()
    {
        return mResponseVariableName;
    }

    std::string DirectSensitivityResponseFunction::GetEvaluationFlag()
    {
        KRATOS_ERROR << "GetEvaluationFlag() needs to be implemented by the derived class.\n";
        std::string flag = "undefined";
        return flag;
    } 
    
};

