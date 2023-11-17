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
class BackwardEulerScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace> {
public:
    BackwardEulerScheme(const Variable<double>& rVariable,
                        const Variable<double>& rDeltaTimeVariable,
                        const Variable<double>& rDeltaTimeVariableCoefficient)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(),
          mVariable(rVariable),
          mDeltaTimeVariable(rDeltaTimeVariable),
          mDeltaTimeVariableCoefficient(rDeltaTimeVariableCoefficient)
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            this->UpdateScalarTimeDerivative(rNode, mVariable, mDeltaTimeVariable);
        });

        KRATOS_CATCH("")
    }

    void UpdateScalarTimeDerivative(Node& rNode,
                                    const Variable<double>& variable,
                                    const Variable<double>& dt_variable) const override
    {
        const auto delta_variable = rNode.FastGetSolutionStepValue(variable, 0) -
                                    rNode.FastGetSolutionStepValue(variable, 1);
        rNode.FastGetSolutionStepValue(dt_variable) =
            delta_variable / this->GetDeltaTime();
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        this->SetDeltaTime(rModelPart.GetProcessInfo()[DELTA_TIME]);
        rModelPart.GetProcessInfo()[mDeltaTimeVariableCoefficient] =
            1.0 / (this->GetDeltaTime());

        KRATOS_CATCH("")
    }

    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        for (const auto& r_node : rModelPart.Nodes()) {
            this->CheckSolutionStepsData(r_node, mVariable);
            this->CheckSolutionStepsData(r_node, mDeltaTimeVariable);
            this->CheckDof(r_node, mVariable);
        }
    }

private:
    Variable<double> mVariable;
    Variable<double> mDeltaTimeVariable;
    Variable<double> mDeltaTimeVariableCoefficient;
};

} // namespace Kratos
