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
        if (mTracedElementIdVector.size() != 0)
            mTracedElementIdVector.resize(0);

        for (auto& elem_i : mpTracedModelPart->Elements())
            mTracedElementIdVector.push_back(elem_i.Id()); 
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

        if (std::find(mTracedElementIdVector.begin(), mTracedElementIdVector.end(),rDirectElement.Id() ) != mTracedElementIdVector.end())    
        {   
            // Calculate Sensitivity Matrix
            if (KratosComponents<Variable<double>>::Has(mDesignVariableName))
            {
                const Variable<double>& r_scalar_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);                                          
                rDirectElement.CalculateSensitivityMatrix(r_scalar_variable, rOutput, rProcessInfo);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName))
            {   
                const Variable<array_1d<double,3>>& r_vector_variable =  KratosComponents<Variable<array_1d<double,3>>>::Get(mDesignVariableName);
                rDirectElement.CalculateSensitivityMatrix(r_vector_variable, rOutput, rProcessInfo);
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << mDesignVariableName << "." << std::endl;
   
            // Check if the sizes of the pseudo load vector and the sensitivity matrix match 
            KRATOS_ERROR_IF_NOT( rPseudoLoadVector.size() == rOutput.size2() ) << 
                "Sizes of the pseudo load vector and the sensitivity matrix do not match" << std::endl;
            
            // Transform rOutput(matrix) from CalculateSensitivityMatrix into rPseudoLoadVector(vector)
            for(IndexType i = 0; i < rLHS.size1(); ++i)                    
                rPseudoLoadVector[i] = rOutput(0, i);           
        }
        else
            rPseudoLoadVector = ZeroVector(rPseudoLoadVector.size()); 

        KRATOS_CATCH("");           
    }

    void DirectSensitivityElementDataVariable::CalculatePseudoLoadVector(Condition& rDirectCondition, const Matrix& rLHS, Vector& rPseudoLoadVector, 
                                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // Sizing of the pseudo load vector
        if (rPseudoLoadVector.size() != rLHS.size2())
            rPseudoLoadVector.resize(rLHS.size2(), false);

        // Set pseudo load to zero    
        rPseudoLoadVector = ZeroVector(rPseudoLoadVector.size());

        KRATOS_CATCH("");
    }
    
    std::vector<unsigned int>  DirectSensitivityElementDataVariable::GetTracedElementId() 
    {              
        return mTracedElementIdVector;
    } 
}    