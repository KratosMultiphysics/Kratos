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
#include "direct_sensitivity_element_data_variable.h"
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/variable_utils.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
    DirectSensitivityElementDataVariable::DirectSensitivityElementDataVariable(ModelPart& rModelPart, Parameters VariableSettings)
    : DirectSensitivityVariable(rModelPart, VariableSettings)
    {
        // Get traced element
        mIdTracedElement = VariableSettings["traced_element_id"].GetInt();
        mpTracedElement = rModelPart.pGetElement(mIdTracedElement);
                
        //Get perturbation size
        mDelta = VariableSettings["delta"].GetDouble();
        
        // Get the name of the design variable 
        mDesignVariableName = VariableSettings["design_variable_name"].GetString();
                
        std::cout << "ch Construction finished" << std::endl;
    }


    DirectSensitivityElementDataVariable::~DirectSensitivityElementDataVariable()
    {
        std::cout << "[DirectSensitivityElementDataVariable] :: Destruction finished" << std::endl;
    }


    void DirectSensitivityElementDataVariable::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

   
    void DirectSensitivityElementDataVariable::CalculatePseudoLoadVector(Element& rDirectElement, const Matrix& rLHS, Vector& rPseudoLoadVector, 
                                            const ProcessInfo& rProcessInfo)
    {
        
        KRATOS_TRY;
        
        // Define working variable rOutput
        Matrix rOutput;
        rOutput.resize(rLHS.size1(), rLHS.size2(), false);

        rOutput.clear();

        // Sizing and clearing of the pseudo load vector
        if (rPseudoLoadVector.size() != rLHS.size1())
            rPseudoLoadVector.resize(rLHS.size1(), false);

        rPseudoLoadVector.clear();

        if (rDirectElement.Id() == mpTracedElement->Id()) 
        {   
            // Set perturbation size         
            rDirectElement.SetValue(PERTURBATION_SIZE, mDelta);

            // Calculate Sensitivity Matrix
            if (KratosComponents<Variable<double>>::Has(mDesignVariableName))                                          
                rDirectElement.CalculateSensitivityMatrix(ReadScalarDesignVariables(mDesignVariableName), rOutput, rProcessInfo);
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName))
                rDirectElement.CalculateSensitivityMatrix(ReadVectorDesignVariables(mDesignVariableName), rOutput, rProcessInfo);
            else
                KRATOS_ERROR << "Unsupported variable: " << mDesignVariableName << "." << std::endl; 
                
            // Transform rOutput from CalculateSensitivityMatrix(matrix) into PseudoLoadVector(vector)
            for(IndexType i = 0; i < rLHS.size1(); ++i)                    
                rPseudoLoadVector[i] = rOutput(0, i);           
        }
        else
        {
            // rPseudoLoadVector = zero vector    
            for(IndexType i = 0; i < rLHS.size1(); ++i)            
                rPseudoLoadVector[i] = 0;                          
        }  
            
        KRATOS_CATCH("");           
    }

    void DirectSensitivityElementDataVariable::CalculatePseudoLoadVector(Condition& rDirectCondition, const Matrix& rLHS, Vector& rPseudoLoadVector, 
                                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // Sizing and clearing of the pseudo load vector
        if (rPseudoLoadVector.size() != rLHS.size1())
            rPseudoLoadVector.resize(rLHS.size1(), false);

        rPseudoLoadVector.clear();

        // rPseudoLoadVector = zero vector    
        for(IndexType i = 0; i < rLHS.size1(); ++i)            
            rPseudoLoadVector[i] = 0;

        KRATOS_CATCH("");
    }


    Variable<double> DirectSensitivityElementDataVariable::ReadScalarDesignVariables(std::string const& rVariableName)
    {
        const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
            return r_variable;                      
    }

    Variable<array_1d<double,3>> DirectSensitivityElementDataVariable::ReadVectorDesignVariables(std::string const& rVariableName)
    {
        const Variable<array_1d<double,3>>& r_variable =
                KratosComponents<Variable<array_1d<double,3>>>::Get(rVariableName);
            return r_variable;                  
    }

    void DirectSensitivityElementDataVariable::PerturbDesignVariable(Element& rDirectElement, Variable<double>& rDesignVariable)                                                 
    {              
        const double current_property_value = rDirectElement.GetProperties()[rDesignVariable];
        rDirectElement.SetValue(rDesignVariable, (current_property_value + mDelta)); 
    }

    void DirectSensitivityElementDataVariable::UnperturbDesignVariable(Element& rDirectElement, Variable<double>& rDesignVariable)                                                  
    {              
        const double current_property_value = rDirectElement.GetProperties()[rDesignVariable];
        rDirectElement.SetValue(rDesignVariable, (current_property_value - mDelta)); 
    }
}    