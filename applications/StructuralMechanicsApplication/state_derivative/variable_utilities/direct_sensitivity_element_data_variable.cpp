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
#include "custom_response_functions/response_utilities/element_finite_difference_utility.h"

namespace Kratos
{
    DirectSensitivityElementDataVariable::DirectSensitivityElementDataVariable(ModelPart& rModelPart, Parameters VariableSettings)
    : DirectSensitivityVariable(rModelPart, VariableSettings)
    {
        KRATOS_TRY;

        // Get traced element
        /*const int mIdTracedElement = VariableSettings["traced_element_id"].GetInt();
        mTracedElement = rModelPart.GetElement(mIdTracedElement);

        //Get perturbation size
        const double mDelta = VariableSettings["delta"].GetDouble();

        // Get design variable
        const Variable<double> mDesignVariable= VariableSettings["design_variable"].GetDouble();*/
        

        std::cout << "DirectSensitivityElementDataVariable wurde konstruiert!!!!!!!!!!" << std::endl;

        KRATOS_CATCH("");
    }

    DirectSensitivityElementDataVariable::~DirectSensitivityElementDataVariable()
    {
        std::cout << "DirectSensitivityElementDataVariable wurde zerstört!!!!!!!!!!" << std::endl;
    }

    void DirectSensitivityElementDataVariable::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }
   
    void DirectSensitivityElementDataVariable::CalculatePseudoLoadVector(Matrix& rPseudoLoadVector, ProcessInfo& rProcessInfo)
    {
        /*KRATOS_TRY;

        // Get design variable
        Variable<double> rDesignVariable = mDesignVariable;
        // Get perturbation size
        const double rPerturbationSize = mDelta;
        // Get process info
        ProcessInfo process_info = rProcessInfo;
        // Get traced Element 
        Element rElem = mTracedElement;
        
        // Get pseudo-load vector using utility
        ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(rElem, rDesignVariable, rPerturbationSize, rPseudoLoadVector, process_info);            
        
        KRATOS_CATCH("");*/
           
    }
 
}    