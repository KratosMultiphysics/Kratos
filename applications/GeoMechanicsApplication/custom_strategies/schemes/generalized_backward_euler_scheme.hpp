// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "geomechanics_time_integration_scheme.hpp"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class GeneralizedBackwardEulerScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace> {
public:
    GeneralizedBackwardEulerScheme(const Variable<double>& rVariable,
                                   const Variable<double>& rDeltaVariable,
                                   const Variable<double>& rDeltaVariableCoefficient)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(),
          mVariable(rVariable),
          mDeltaVariable(rDeltaVariable),
          mDeltaVariableCoefficient(rDeltaVariableCoefficient)
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            this->UpdateScalarTimeDerivative(rNode, mVariable,
                                             mDeltaVariable);
        });

        KRATOS_CATCH("")
    }

    void UpdateScalarTimeDerivative(Node& rNode,
                                    const Variable<double>& variable,
                                    const Variable<double>& dt_variable) const override
    {
        const double delta_variable = rNode.FastGetSolutionStepValue(variable) -
                                      rNode.FastGetSolutionStepValue(variable, 1);
        rNode.FastGetSolutionStepValue(dt_variable) =
            delta_variable / this->GetDeltaTime();
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        this->SetDeltaTime(rModelPart.GetProcessInfo()[DELTA_TIME]);
        rModelPart.GetProcessInfo()[mDeltaVariableCoefficient] =
            1.0 / (this->GetDeltaTime());

        KRATOS_CATCH("")
    }

    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        for (const auto& rNode : rModelPart.Nodes()) {
            this->CheckSolutionStepsData(rNode, mVariable);
            this->CheckSolutionStepsData(rNode, mDeltaVariable);
            this->CheckDof(rNode, mVariable);
        }
    }

private:
    Variable<double> mVariable;
    Variable<double> mDeltaVariable;
    Variable<double> mDeltaVariableCoefficient;
};

} // namespace Kratos
