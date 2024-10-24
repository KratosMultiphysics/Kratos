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

    void InitializeNonLinIteration(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        // only initialize non linear iteration conditions
        this->BlockForEachActiveCondition(rModelPart, &Condition::InitializeNonLinearIteration);

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
        // no prediction is used in the scheme
    }

    void Update(ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType&, TSystemVectorType& Dx, TSystemVectorType&) override
    {
        KRATOS_TRY

        // only update derivatives, solution step is updated in the strategy.
        // Note that this workflow of only updating the derivatives is specific for this scheme.
        TSystemVectorType first_derivative_vector;
        TSystemVectorType second_derivative_vector;

        Geo::SparseSystemUtilities::GetUFirstAndSecondDerivativeVector(
            first_derivative_vector, second_derivative_vector, rDofSet, rModelPart, 0);

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

}; // Class IncrementalNewmarkLinearElasticUScheme

} // namespace Kratos
