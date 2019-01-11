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
#include "direct_sensitivity_variable.h"
#include "utilities/variable_utils.h"
#include "custom_response_functions/response_utilities/element_finite_difference_utility.h"

namespace Kratos
{

    /// Constructor.
    DirectSensitivityVariable::DirectSensitivityVariable(ModelPart& rModelPart, Parameters VariableSettings)
      : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        //Get perturbation size
        //mDelta = VariableSettings["delta"].GetDouble();

        std::cout << "[DirectSensitivityDataVariable] :: Construction finished" << std::endl;

        KRATOS_CATCH("");
    }

    /// Destructor.
    DirectSensitivityVariable::~DirectSensitivityVariable()
    {
        std::cout << "[DirectSensitivityDataVariable] :: Destruction finished" << std::endl;
    }

    ///@}
    ///@name Operations
    ///@{
        
    void DirectSensitivityVariable::Initialize()
    {
        KRATOS_TRY;
        
        KRATOS_CATCH("");
    }

      
    void DirectSensitivityVariable::CalculatePseudoLoadVector(Element& rDirectElement, const Matrix& rRHS, 
                                    Vector& rPseudoLoadVector, const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "CalculatePseudoLoadVector should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }  

    void DirectSensitivityVariable::CalculatePseudoLoadVector(Condition& rDirectCondition, const Matrix& rLHS, 
                                    Vector& rPseudoLoadVector, const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "CalculatePseudoLoadVector should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }  

    void DirectSensitivityVariable::PerturbDesignVariable(Element& rElement, Variable<double>& rDesignVariable)
    {
        KRATOS_ERROR << "PerturbDesignVariable should be implemented in the derived class." << std::endl;
    }
    
    void DirectSensitivityVariable::UnperturbDesignVariable(Element& rElement, Variable<double>& rDesignVariable)
    {
        KRATOS_ERROR << "UnperturbDesignVariable should be implemented in the derived class." << std::endl;
    }

    Variable<double> DirectSensitivityVariable::ReadScalarDesignVariables(std::string const& rVariableName)
    {
        KRATOS_ERROR << "ReadScalarDesignVariables should be implemented in the derived class." << std::endl;
    }

    Variable<array_1d<double,3>> DirectSensitivityVariable::ReadVectorDesignVariables(std::string const& rVariableName)
    {
        KRATOS_ERROR << "ReadVectorDesignVariables should be implemented in the derived class." << std::endl;
    }
}   