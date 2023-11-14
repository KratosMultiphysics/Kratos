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
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "generalized_newmark_scheme.hpp"
#include "geo_mechanics_application_variables.h"
#include "geomechanics_time_integration_scheme.hpp"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerTScheme : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerTScheme);

    BackwardEulerTScheme() : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>()
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            const double delta_temperature =
                rNode.FastGetSolutionStepValue(TEMPERATURE) -
                rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) =
                delta_temperature / this->GetDeltaTime();
        });

        KRATOS_CATCH("")
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        this->SetDeltaTime(rModelPart.GetProcessInfo()[DELTA_TIME]);
        rModelPart.GetProcessInfo()[DT_TEMPERATURE_COEFFICIENT] =
            1.0 / (this->GetDeltaTime());

        KRATOS_CATCH("")
    }

    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        for (const auto& r_node : rModelPart.Nodes()) {
            this->CheckSolutionStepsData(r_node, TEMPERATURE);
            this->CheckSolutionStepsData(r_node, DT_TEMPERATURE);
            this->CheckDof(r_node, TEMPERATURE);
        }
    }

}; // Class BackwardEulerTScheme
} // namespace Kratos
