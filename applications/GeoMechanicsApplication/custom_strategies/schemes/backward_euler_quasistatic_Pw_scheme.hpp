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
#include "geomechanics_time_integration_scheme.hpp"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerQuasistaticPwScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerQuasistaticPwScheme);

    BackwardEulerQuasistaticPwScheme()
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>()
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        // Update DtPressure
        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            const double delta_pressure =
                rNode.FastGetSolutionStepValue(WATER_PRESSURE) -
                rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1);
            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) =
                delta_pressure / this->GetDeltaTime();
        });

        KRATOS_CATCH("")
    }

    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        for (const auto& rNode : rModelPart.Nodes()) {
            this->CheckSolutionStepsData(rNode, WATER_PRESSURE);
            this->CheckSolutionStepsData(rNode, DT_WATER_PRESSURE);
            this->CheckDof(rNode, WATER_PRESSURE);
        }
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        this->SetDeltaTime(rModelPart.GetProcessInfo()[DELTA_TIME]);
        rModelPart.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] =
            1.0 / (this->GetDeltaTime());

        KRATOS_CATCH("")
    }

}; // Class BackwardEulerQuasistaticPwScheme

} // namespace Kratos