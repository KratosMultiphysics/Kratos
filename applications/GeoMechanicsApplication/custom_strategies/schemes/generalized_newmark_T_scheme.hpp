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
#include "geo_mechanics_application_variables.h"
#include "geomechanics_time_integration_scheme.hpp"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class GeneralizedNewmarkTScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace> {
public:
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    KRATOS_CLASS_POINTER_DEFINITION(GeneralizedNewmarkTScheme);

    explicit GeneralizedNewmarkTScheme(double theta)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(), mTheta(theta)
    {
    }

    ~GeneralizedNewmarkTScheme() override = default;

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        Scheme<TSparseSpace, TDenseSpace>::Check(rModelPart);
        CheckAllocatedVariables(rModelPart);
        this->CheckBufferSize(rModelPart);

        KRATOS_ERROR_IF(mTheta <= 0) << "Theta has an invalid value\n";

        return 0;

        KRATOS_CATCH("")
    }

    void CheckAllocatedVariables(const ModelPart& rModelPart) const
    {
        for (const auto& r_node : rModelPart.Nodes()) {
            this->CheckSolutionStepsData(r_node, TEMPERATURE);
            this->CheckSolutionStepsData(r_node, DT_TEMPERATURE);
            this->CheckDof(r_node, TEMPERATURE);
        }
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY
        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            const double delta_temperature =
                rNode.FastGetSolutionStepValue(TEMPERATURE) -
                rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            const auto& previous_dt_temperature =
                rNode.FastGetSolutionStepValue(DT_TEMPERATURE, 1);

            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) =
                (delta_temperature - (1.0 - mTheta) * this->GetDeltaTime() * previous_dt_temperature) /
                (mTheta * this->GetDeltaTime());
        });
        KRATOS_CATCH("")
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        this->SetDeltaTime(rModelPart.GetProcessInfo()[DELTA_TIME]);
        rModelPart.GetProcessInfo()[DT_TEMPERATURE_COEFFICIENT] =
            1.0 / (mTheta * this->GetDeltaTime());

        KRATOS_CATCH("")
    }

private:
    double mTheta = 0.0;
}; // Class GeneralizedNewmarkTScheme
} // namespace Kratos