//
// Author: Guillermo Casas, gcasas@cimne.upc.edu
//

#include "adams_bashforth_strategy.h"

namespace Kratos {

double AdamsBashforthStrategy::Solve() {
    KRATOS_TRY
    ModelPart& r_model_part = GetModelPart();
    //int step = r_model_part.GetProcessInfo()[FRACTIONAL_STEP]
    if (mFirstStep){
        mFirstStep = false;
        PerformTimeIntegrationOfMotion(1);
    }

    else {
        mFirstStep = true;
        InitializeSolutionStep();
        SearchDEMOperations(r_model_part);
        SearchFEMOperations(r_model_part);
        ForceOperations(r_model_part);
        PerformTimeIntegrationOfMotion(2);
        FinalizeSolutionStep();
    }
    return 0.00;

    KRATOS_CATCH("")
}

} // namespace Kratos
