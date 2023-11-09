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
#include "geo_mechanics_scheme.hpp"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticTScheme : public GeoMechanicsScheme<TSparseSpace, TDenseSpace> {
public:

    using BaseType              = Scheme<TSparseSpace,TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    KRATOS_CLASS_POINTER_DEFINITION(NewmarkQuasistaticTScheme);

    explicit NewmarkQuasistaticTScheme(double theta)
        : GeoMechanicsScheme<TSparseSpace, TDenseSpace>(), mTheta(theta)
    {
    }

    ~NewmarkQuasistaticTScheme() override = default;

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
        for (const auto& rNode : rModelPart.Nodes()) {
            KRATOS_ERROR_IF(!rNode.SolutionStepsDataHas(TEMPERATURE))
                << "TEMPERATURE variable is not allocated for node "
                << rNode.Id() << std::endl;

            KRATOS_ERROR_IF(!rNode.SolutionStepsDataHas(DT_TEMPERATURE))
                << "DT_TEMPERATURE variable is not allocated for node "
                << rNode.Id() << std::endl;

            KRATOS_ERROR_IF(!rNode.HasDofFor(TEMPERATURE))
                << "missing TEMPERATURE dof on node " << rNode.Id() << std::endl;
        }
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              TSystemMatrixType& A,
                              TSystemVectorType& Dx,
                              TSystemVectorType& b) override
    {
        this->FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY
        // Update DtPressure
        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            const double delta_temperature =
                rNode.FastGetSolutionStepValue(TEMPERATURE) -
                rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            const auto& previous_dt_temperature =
                rNode.FastGetSolutionStepValue(DT_TEMPERATURE, 1);

            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) =
                1.0 / (mTheta * this->GetDeltaTime()) *
                (delta_temperature - (1.0 - mTheta) * this->GetDeltaTime() * previous_dt_temperature);
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
}; // Class NewmarkQuasistaticTScheme
} // namespace Kratos