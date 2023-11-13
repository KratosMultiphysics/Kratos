// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Richard Faasse
//

#pragma once

#include "custom_strategies/schemes/newmark_T_scheme.hpp"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerTScheme
    : public NewmarkTScheme<TSparseSpace, TDenseSpace> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerTScheme);

    BackwardEulerTScheme()
        : NewmarkTScheme<TSparseSpace, TDenseSpace>(1.0)
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            const double DeltaTemperature =
                rNode.FastGetSolutionStepValue(TEMPERATURE) -
                rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) = DeltaTemperature / this->GetDeltaTime();
        });

        KRATOS_CATCH("")
    }

}; // Class BackwardEulerTScheme
} // namespace Kratos
