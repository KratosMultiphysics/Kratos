// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_utilities/sparse_system_utilities.h"
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class IncrementalNewmarkLinearElasticUScheme : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(IncrementalNewmarkLinearElasticUScheme);

    using BaseType              = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    IncrementalNewmarkLinearElasticUScheme(double beta, double gamma)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              {}, {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)}, beta, gamma, std::nullopt)
    {
    }

    void InitializeSolutionStep(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        this->SetTimeFactors(rModelPart);

        // only inititialize solutionstep conditions
        this->BlockForEachActiveCondition(rModelPart, &Condition::InitializeSolutionStep);

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        // only finalize solutionstep conditions
        this->BlockForEachActiveCondition(rModelPart, &Condition::FinalizeSolutionStep);

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        // only finalize non linear iteration conditions
        this->BlockForEachActiveCondition(rModelPart, &Condition::FinalizeNonLinearIteration);

        KRATOS_CATCH("")
    }

    void Predict(ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b) override
    {
        KRATOS_TRY

        // Update (Angular) Acceleration, (Angular) Velocity and DtPressure
        this->PredictVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }

    void Update(ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType&, TSystemVectorType& Dx, TSystemVectorType&) override
    {
        KRATOS_TRY

        // only update derivatives, solution step is updated in the strategy
        TSystemVectorType first_derivative_vector;
        TSystemVectorType second_derivative_vector;

        Geo::SparseSystemUtilities::GetUFirstAndSecondDerivativeVector(
            first_derivative_vector, second_derivative_vector, rDofSet, rModelPart, 1);

        TSystemVectorType delta_first_derivative_vector = TSystemVectorType(rDofSet.size(), 0.0);
        TSparseSpace::UnaliasedAdd(delta_first_derivative_vector,
                                   (this->GetGamma() / (this->GetBeta() * this->GetDeltaTime())), Dx);
        TSparseSpace::UnaliasedAdd(delta_first_derivative_vector,
                                   -(this->GetGamma() / this->GetBeta()), first_derivative_vector);
        TSparseSpace::UnaliasedAdd(delta_first_derivative_vector,
                                   this->GetDeltaTime() * (1 - this->GetGamma() / (2 * this->GetBeta())),
                                   second_derivative_vector);

        // performs: TSystemVectorType delta_second_derivative_vector = Dx * (1 / (mBeta * delta_time * delta_time)) - first_derivative_vector * (1 / (mBeta * delta_time)) - r_second_derivative_vector * (1 / (2 * mBeta));
        TSystemVectorType delta_second_derivative_vector = TSystemVectorType(rDofSet.size(), 0.0);
        TSparseSpace::UnaliasedAdd(delta_second_derivative_vector,
                                   1 / (this->GetBeta() * this->GetDeltaTime() * this->GetDeltaTime()), Dx);
        TSparseSpace::UnaliasedAdd(delta_second_derivative_vector,
                                   -1 / (this->GetBeta() * this->GetDeltaTime()), first_derivative_vector);
        TSparseSpace::UnaliasedAdd(delta_second_derivative_vector, -1 / (2 * this->GetBeta()),
                                   second_derivative_vector);

        // performs: first_derivative_vector += delta_first_derivative_vector;
        TSparseSpace::UnaliasedAdd(first_derivative_vector, 1.0, delta_first_derivative_vector);

        // performs: second_derivative_vector += delta_second_derivative_vector;
        TSparseSpace::UnaliasedAdd(second_derivative_vector, 1.0, delta_second_derivative_vector);

        Geo::SparseSystemUtilities::SetUFirstAndSecondDerivativeVector(
            first_derivative_vector, second_derivative_vector, rModelPart);

        KRATOS_CATCH("")
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        // Empty function because the derivatives are updated outside of the function since it needs more input. And this function is required to be defined
    }

private:
    inline void PredictVariablesDerivatives(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // make sure the timestep size is correct
        this->SetTimeFactors(rModelPart);

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
                if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

                // firstly get current values of first and second derivative
                const array_1d<double, 3> first_derivative_array = rNode.FastGetSolutionStepValue(
                    r_second_order_vector_variable.first_time_derivative, 0);

                const array_1d<double, 3> second_derivative_array = rNode.FastGetSolutionStepValue(
                    r_second_order_vector_variable.second_time_derivative, 0);

                // secondly calculate prediction of first and second derivative
                const array_1d<double, 3> predicted_first_time_derivative =
                    first_derivative_array * (this->GetGamma() / this->GetBeta()) +
                    second_derivative_array *
                        (this->GetDeltaTime() * (this->GetGamma() / (2 * this->GetBeta()) - 1));

                const array_1d<double, 3> predicted_second_time_derivative =
                    first_derivative_array * (1.0 / (this->GetBeta() * this->GetDeltaTime())) +
                    second_derivative_array * (1.0 / (2.0 * this->GetBeta()));

                // update values after calculating both the first and second derivative prediction
                NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                    rNode, r_second_order_vector_variable.first_time_derivative, predicted_first_time_derivative);

                NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                    rNode, r_second_order_vector_variable.second_time_derivative, predicted_second_time_derivative);
            }
        });

        KRATOS_CATCH("")
    }
}; // Class IncrementalNewmarkLinearElasticUScheme

} // namespace Kratos
