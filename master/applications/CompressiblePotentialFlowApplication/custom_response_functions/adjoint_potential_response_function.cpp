//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nunez, based on Martin Fusseder work, https://github.com/MFusseder
//


// System includes

// External includes

// Project includes
#include "adjoint_potential_response_function.h"
#include "utilities/variable_utils.h"


namespace Kratos
{

    /// Constructor.
    AdjointPotentialResponseFunction::AdjointPotentialResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
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
        else if (gradient_mode == "analytic")
        {
            mGradientMode = 2;
        }
        else
            KRATOS_ERROR << "Specified gradient_mode not recognized. The only "
                            "options are: semi_analytic and analytic. "
                            "Specified gradient_mode: "
                         << gradient_mode << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointPotentialResponseFunction::Initialize()
    {
        KRATOS_TRY;

        if(mGradientMode == 1)
        {
            VariableUtils().SetNonHistoricalVariable(SCALE_FACTOR, mDelta, mrModelPart.Elements());
            VariableUtils().SetNonHistoricalVariable(SCALE_FACTOR, mDelta, mrModelPart.Conditions());
        }

        KRATOS_CATCH("");
    }

    void AdjointPotentialResponseFunction::InitializeSolutionStep()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    void AdjointPotentialResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
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

    double AdjointPotentialResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_ERROR << "CalculateValue needs to be implemented by the derived class.\n";
    }
};


