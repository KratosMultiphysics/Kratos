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


namespace Kratos
{

    /// Constructor.
    DirectSensitivityVariable::DirectSensitivityVariable(ModelPart& rModelPart, Parameters VariableSettings)
      : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        //Get perturbation size
        //const double mDelta = VariableSettings["delta"].GetDouble();

        std::cout << "DirectSensitivityVariable wurde konstruiert!!!!!!!!!" << std::endl;

        KRATOS_CATCH("");
    }

    /// Destructor.
    DirectSensitivityVariable::~DirectSensitivityVariable()
    {
    }

    ///@}
    ///@name Operations
    ///@{
        
    void DirectSensitivityVariable::Initialize()
    {
        KRATOS_TRY;
        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;
        KRATOS_CATCH("");
    }

      
    void DirectSensitivityVariable::CalculatePseudoLoadVector( Matrix& rPseudoLoadVector, ProcessInfo& rProcessInfo)
    {
        /*KRATOS_TRY;
        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;
        KRATOS_CATCH("");*/
    }  
}   