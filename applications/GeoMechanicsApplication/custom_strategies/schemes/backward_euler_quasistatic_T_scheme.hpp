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

// Project includes
#include "custom_strategies/schemes/newmark_quasistatic_T_scheme.hpp"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>

class BackwardEulerQuasistaticTScheme
    : public NewmarkQuasistaticTScheme<TSparseSpace, TDenseSpace> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerQuasistaticTScheme);

    BackwardEulerQuasistaticTScheme()
        : NewmarkQuasistaticTScheme<TSparseSpace, TDenseSpace>(1.0)
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
            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) = DeltaTemperature / mDeltaTime;
        });

        KRATOS_CATCH("")
    }

}; // Class BackwardEulerQuasistaticTScheme
} // namespace Kratos
