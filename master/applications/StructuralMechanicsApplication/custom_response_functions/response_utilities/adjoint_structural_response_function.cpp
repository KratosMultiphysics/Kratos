// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
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
      : mrModelPart(rModelPart),
        mResponseSettings(ResponseSettings)
    {
        KRATOS_TRY;

        // Set gradient mode
        const std::string& gradient_mode = ResponseSettings["gradient_mode"].GetString();

        // Mode 1: semi-analytic sensitivities
        if (gradient_mode == "semi_analytic")
            mGradientMode = 1;
        else
            KRATOS_ERROR << "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: " <<  gradient_mode << std::endl;


        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::Initialize()
    {
        KRATOS_TRY;

        if(mGradientMode == 1)
        {
            double perturbation_size = mResponseSettings["step_size"].GetDouble();
            mrModelPart.GetProcessInfo()[PERTURBATION_SIZE] = perturbation_size;

            bool adapt_perturbation_size = false;
            if(mResponseSettings.Has("adapt_step_size"))
                adapt_perturbation_size = mResponseSettings["adapt_step_size"].GetBool();
            mrModelPart.GetProcessInfo()[ADAPT_PERTURBATION_SIZE] = adapt_perturbation_size;
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


