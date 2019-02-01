// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
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
    DirectSensitivityResponseFunction::DirectSensitivityResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
      : mrModelPart(rModelPart)
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
                                            const Matrix& rLHS,
                                            Matrix& rResponseGradientMatrix,                                            
                                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculateGradient needs to be implemented by the derived class.\n";
    }
    
    void DirectSensitivityResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                            DirectSensitivityVariable& DesignVariable, 
                                            Matrix& rSensitivityGradient,
                                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "CalculatePartialGradient needs to be implemented by the derived class.\n";
    }

    int DirectSensitivityResponseFunction::GetNumberOfTracedGaussPoints()
    {
        KRATOS_ERROR << "GetNumberOfTracedGaussPoints() needs to be implemented by the derived class.\n";
        return 0;
    }
};

