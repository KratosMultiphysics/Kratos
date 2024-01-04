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
class BackwardEulerScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    BackwardEulerScheme(const std::vector<FirstOrderVariable>& rFirstOrderVariables,
                        const std::vector<VariableWithTimeDerivatives>& rVariablesWithDerivatives)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(
              rFirstOrderVariables, rVariablesWithDerivatives)
    {
    }

protected:
    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const auto delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        this->SetDeltaTime(delta_time);

        for (const auto& r_first_order_variable : this->GetFirstOrderVariables())
        {
            rModelPart.GetProcessInfo()[r_first_order_variable.delta_time_coefficient] =
                1.0 / delta_time;
        }

        KRATOS_CATCH("")
    }

    void UpdateVectorFirstTimeDerivative(Node& rNode) const override
    {
        for (const auto& r_variable_with_derivative : this->GetSecondOrderVariables())
        {
            SetDerivative(r_variable_with_derivative.first_time_derivative,
                          r_variable_with_derivative.instance, rNode);
        }
    }

    void UpdateVectorSecondTimeDerivative(Node& rNode) const override
    {
        for (const auto& r_variable_with_derivative : this->GetSecondOrderVariables())
        {
            // Make sure that setting the second_time_derivative is done
            // after setting the first_time_derivative.
            SetDerivative(r_variable_with_derivative.second_time_derivative,
                          r_variable_with_derivative.first_time_derivative, rNode);
        }
    }

    template <class T>
    void SetDerivative(const Variable<T>& derivative_variable,
                       const Variable<T>& instance_variable,
                       Node& rNode) const
    {
        rNode.FastGetSolutionStepValue(derivative_variable) =
            (rNode.FastGetSolutionStepValue(instance_variable, 0) -
             rNode.FastGetSolutionStepValue(instance_variable, 1)) /
            this->GetDeltaTime();
    }

    void UpdateScalarTimeDerivative(Node& rNode,
                                    const Variable<double>& variable,
                                    const Variable<double>& dt_variable) const override
    {
        SetDerivative(dt_variable, variable, rNode);
    }
};

} // namespace Kratos
