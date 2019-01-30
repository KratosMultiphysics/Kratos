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


namespace Kratos
{
    // Constructor
    DirectSensitivityElementDataVariable::DirectSensitivityElementDataVariable(ModelPart& rModelPart, Parameters VariableSettings)
    : DirectSensitivityVariable(rModelPart, VariableSettings)
    {
        // Get traced element
        mIdTracedElement = VariableSettings["traced_element_id"].GetInt();
        mpTracedElement = rModelPart.pGetElement(mIdTracedElement);
    }

    // Destructor
    DirectSensitivityElementDataVariable::~DirectSensitivityElementDataVariable(){}


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
   
            // Check if the sizes of the pseudo load vector and the sensitivity matrix match 
            KRATOS_ERROR_IF_NOT( rPseudoLoadVector.size() == rOutput.size2() ) << 
                "Sizes of the pseudo load vector and the sensitivity matrix do not match" << std::endl;
            
            // Transform rOutput from CalculateSensitivityMatrix(matrix) into rPseudoLoadVector(vector)
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

        // Sizing of the pseudo load vector
        if (rPseudoLoadVector.size() != rLHS.size2())
            rPseudoLoadVector.resize(rLHS.size2(), false);

        // rPseudoLoadVector = zero vector    
        for(IndexType i = 0; i < rLHS.size2(); ++i)            
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

    
    void DirectSensitivityElementDataVariable::PerturbDesignVariable(Element& rDirectElement)                                                 
    {              
        if (KratosComponents<Variable<double>>::Has(mDesignVariableName))
        {
            // Save property pointer
            Properties::Pointer p_global_properties = rDirectElement.pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rDirectElement.SetProperties(p_local_property);

            // perturb the design variable
            const double current_property_value = rDirectElement.GetProperties()[ ReadScalarDesignVariables(mDesignVariableName) ];
            p_local_property->SetValue( ReadScalarDesignVariables(mDesignVariableName) , (current_property_value + mDelta) );
        }
        else if (KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName))                      
            KRATOS_ERROR << "No element_data_variable of type Variable<array_1d<double,3>> is implemented yet." << std::endl;
        else
            KRATOS_ERROR << "Unsupported variable: " << mDesignVariableName << "." << std::endl;         
    }    

    
    void DirectSensitivityElementDataVariable::UnperturbDesignVariable(Element& rDirectElement)                                                  
    {              
        if (KratosComponents<Variable<double>>::Has(mDesignVariableName))
        {
            // Save property pointer
            Properties::Pointer p_global_properties = rDirectElement.pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rDirectElement.SetProperties(p_local_property);

            // perturb the design variable
            const double current_property_value = rDirectElement.GetProperties()[ ReadScalarDesignVariables(mDesignVariableName) ];
            p_local_property->SetValue( ReadScalarDesignVariables(mDesignVariableName) , (current_property_value - mDelta) );
        }
        else if (KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName))                      
            KRATOS_ERROR << "No element_data_variable of type Variable<array_1d<double,3>> is implemented yet." << std::endl;
        else
            KRATOS_ERROR << "Unsupported variable: " << mDesignVariableName << "." << std::endl; 
    }  
    
    unsigned int  DirectSensitivityElementDataVariable::GetTracedElementId() 
    {              
        return mIdTracedElement;
    } 
}    