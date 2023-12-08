// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "newmark_quasistatic_U_Pw_scheme.hpp"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerQuasistaticUPwScheme
    : public NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerQuasistaticUPwScheme);

    BackwardEulerQuasistaticUPwScheme()
        : NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>(1.0, 1.0, 1.0)
    {
    }

protected:
    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        this->SetDeltaTime(rModelPart.GetProcessInfo()[DELTA_TIME]);
        rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] = 1.0 / this->GetDeltaTime();
        rModelPart.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0 / this->GetDeltaTime();

        KRATOS_CATCH("")
    }

    template <class T>
    void SetDerivative(const Variable<T>& derivative_variable,
                       const Variable<T>& instance_variable,
                       Node& rNode)
    {
        rNode.FastGetSolutionStepValue(derivative_variable) =
            (rNode.FastGetSolutionStepValue(instance_variable, 0) -
             rNode.FastGetSolutionStepValue(instance_variable, 1)) /
            this->GetDeltaTime();
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        // Update Acceleration, Velocity and DtPressure
        block_for_each(rModelPart.Nodes(), [this](Node& rNode)
        {
            for (const auto& r_variable_with_derivative : this->GetVariableDerivatives())
            {
                SetDerivative(r_variable_with_derivative.first_time_derivative,
                              r_variable_with_derivative.instance, rNode);

                // Make sure that setting the second_time_derivative is done
                // after setting the first_time_derivative.
                SetDerivative(r_variable_with_derivative.second_time_derivative,
                              r_variable_with_derivative.first_time_derivative, rNode);
            }

            SetDerivative(DT_WATER_PRESSURE, WATER_PRESSURE, rNode);
        });

        KRATOS_CATCH("")
    }

    std::string Info() const override
    {
        return "BackwardEulerQuasistaticUPwScheme";
    }
}; // Class BackwardEulerQuasistaticUPwScheme

} // namespace Kratos