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
#include "backward_euler_scheme.hpp"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerQuasistaticUPwScheme : public BackwardEulerScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerQuasistaticUPwScheme);

    BackwardEulerQuasistaticUPwScheme()
        : BackwardEulerScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)})
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            // static here means no velocities and accelerations for the displacement/rotation D.O.F.
            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                if (rNode.IsFixed(r_first_order_scalar_variable.first_time_derivative)) continue;

                rNode.FastGetSolutionStepValue(r_first_order_scalar_variable.first_time_derivative) =
                    this->CalculateDerivative(r_first_order_scalar_variable.instance, rNode);
            }
        });

        KRATOS_CATCH("")
    }

    std::string Info() const override { return "BackwardEulerQuasistaticUPwScheme"; }
}; // Class BackwardEulerQuasistaticUPwScheme

} // namespace Kratos