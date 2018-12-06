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
#include "adjoint_structural_response_function.h"
#include "utilities/variable_utils.h"


namespace Kratos
{

    /// Constructor.
    AdjointStructuralResponseFunction::AdjointStructuralResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
      : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        // Set gradient mode
        const std::string& gradient_mode = ResponseSettings["gradient_mode"].GetString();

        // Mode 1: semi-analytic sensitivities
        if (gradient_mode == "semi_analytic")
        {
            mGradientMode = 1;
            double delta = ResponseSettings["step_size"].GetDouble();
            mDelta = delta;
        }
        else
            KRATOS_ERROR << "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: " <<  gradient_mode << std::endl;


        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::Initialize()
    {
        KRATOS_TRY;

        if(mGradientMode == 1)
        {
            VariableUtils().SetNonHistoricalVariable(PERTURBATION_SIZE, mDelta, mrModelPart.Elements());
            VariableUtils().SetNonHistoricalVariable(PERTURBATION_SIZE, mDelta, mrModelPart.Conditions());
        }

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    double AdjointStructuralResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_ERROR << "CalculateValue needs to be implemented by the derived class.\n";
    }
};


