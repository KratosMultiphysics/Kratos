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
#include "utilities/parallel_utilities.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class BackwardEulerQuasistaticUPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( BackwardEulerQuasistaticUPwScheme );

    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;

    BackwardEulerQuasistaticUPwScheme() :
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(1.0, 1.0, 1.0)
    {
    }

protected:
    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        mDeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];
        rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT]    = 1.0/mDeltaTime;
        rModelPart.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0/mDeltaTime;

        KRATOS_CATCH("")
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure
        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            // refactor, extract the (a -b)/mDeltaTime that happens 3 times here
            noalias(rNode.FastGetSolutionStepValue(VELOCITY))     = (  rNode.FastGetSolutionStepValue(DISPLACEMENT)
                                                                                   - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1)) / mDeltaTime;

            noalias(rNode.FastGetSolutionStepValue(ACCELERATION)) = (  rNode.FastGetSolutionStepValue(VELOCITY)
                                                                                   - rNode.FastGetSolutionStepValue(VELOCITY,1) ) / mDeltaTime;

            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE)        = (  rNode.FastGetSolutionStepValue(WATER_PRESSURE)
                                                                                  - rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1)) / mDeltaTime;
        });

        KRATOS_CATCH( "" )
    }
}; // Class BackwardEulerQuasistaticUPwScheme

}