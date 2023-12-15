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

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeneralizedNewmarkScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    GeneralizedNewmarkScheme(double theta,
                             const Variable<double>& rVariable,
                             const Variable<double>& rDeltaTimeVariable,
                             const Variable<double>& rDeltaTimeVariableCoefficient)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(),
          mTheta(theta),
          mVariable(rVariable),
          mDeltaTimeVariable(rDeltaTimeVariable),
          mDeltaTimeVariableCoefficient(rDeltaTimeVariableCoefficient)
    {
        KRATOS_ERROR_IF(this->mTheta <= 0)
            << "Theta must be larger than zero, but got " << this->mTheta << "\n";
    }

protected:
    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        for (const auto& r_node : rModelPart.Nodes())
        {
            this->CheckSolutionStepsData(r_node, mVariable);
            this->CheckSolutionStepsData(r_node, mDeltaTimeVariable);
            this->CheckDof(r_node, mVariable);
        }
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const auto delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        this->SetDeltaTime(delta_time);
        rModelPart.GetProcessInfo()[mDeltaTimeVariableCoefficient] =
            1.0 / (mTheta * delta_time);

        KRATOS_CATCH("")
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(),
                       [this](Node& rNode) {
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
        const auto previous_dt_variable = rNode.FastGetSolutionStepValue(dt_variable, 1);

        rNode.FastGetSolutionStepValue(dt_variable, 0) =
            (delta_variable - (1.0 - mTheta) * this->GetDeltaTime() * previous_dt_variable) /
            (mTheta * this->GetDeltaTime());
    }

private:
    double mTheta = 1.0;
    Variable<double> mVariable;
    Variable<double> mDeltaTimeVariable;
    Variable<double> mDeltaTimeVariableCoefficient;
};

} // namespace Kratos
