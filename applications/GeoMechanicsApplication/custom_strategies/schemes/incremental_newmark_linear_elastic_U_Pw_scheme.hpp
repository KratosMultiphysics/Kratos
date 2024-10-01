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
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class IncrementalNewmarkLinearElasticUPwScheme : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(IncrementalNewmarkLinearElasticUPwScheme);

    using BaseType              = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    IncrementalNewmarkLinearElasticUPwScheme(double beta, double gamma, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)},
              beta,
              gamma,
              theta)
    {
        int num_threads = ParallelUtilities::GetNumThreads();
        mMassMatrix.resize(num_threads);
        mAccelerationVector.resize(num_threads);
        mDampingMatrix.resize(num_threads);
        mVelocityVector.resize(num_threads);
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
        TSystemVectorType first_derivative_vector = TSystemVectorType(rDofSet.size(), 0.0);
        TSystemVectorType second_derivative_vector = TSystemVectorType(rDofSet.size(), 0.0);


        this->GetFirstAndSecondDerivativeVector(first_derivative_vector,
            second_derivative_vector,
            rModelPart, 1);

        TSystemVectorType delta_first_derivative_vector = TSystemVectorType(rDofSet.size(), 0.0);
        TSparseSpace::UnaliasedAdd(delta_first_derivative_vector, (GetGamma() / (GetBeta() * GetDeltaTime())), Dx);
        TSparseSpace::UnaliasedAdd(delta_first_derivative_vector, -(GetGamma() / GetBeta()), first_derivative_vector);
        TSparseSpace::UnaliasedAdd(delta_first_derivative_vector,
            GetDeltaTime() * (1 - GetGamma() / (2 * GetBeta())), second_derivative_vector);


        // performs: TSystemVectorType delta_second_derivative_vector = Dx * (1 / (mBeta * delta_time * delta_time)) - first_derivative_vector * (1 / (mBeta * delta_time)) - r_second_derivative_vector * (1 / (2 * mBeta));
        TSystemVectorType delta_second_derivative_vector = TSystemVectorType(rDofSet.size(), 0.0);
        TSparseSpace::UnaliasedAdd(delta_second_derivative_vector, 1 / (GetBeta() * GetDeltaTime() * GetDeltaTime()), Dx);
        TSparseSpace::UnaliasedAdd(delta_second_derivative_vector, -1 / (GetBeta() * GetDeltaTime()), first_derivative_vector);
        TSparseSpace::UnaliasedAdd(delta_second_derivative_vector, -1 / (2 * GetBeta()), second_derivative_vector);

        // performs: first_derivative_vector += delta_first_derivative_vector;
        TSparseSpace::UnaliasedAdd(first_derivative_vector, 1.0, delta_first_derivative_vector);

        // performs: second_derivative_vector += delta_second_derivative_vector;
        TSparseSpace::UnaliasedAdd(second_derivative_vector, 1.0, delta_second_derivative_vector);

        this->SetFirstAndSecondDerivativeVector(first_derivative_vector, second_derivative_vector, rModelPart);



        //UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(Condition&                     rCurrentCondition,
                                      LocalSystemMatrixType&         LHS_Contribution,
                                      LocalSystemVectorType&         RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(Element&                       rCurrentElement,
                                      LocalSystemMatrixType&         LHS_Contribution,
                                      LocalSystemVectorType&         RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Element&                       rCurrentElement,
                                  LocalSystemVectorType&         RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Condition&                     rCurrentCondition,
                                  LocalSystemVectorType&         rRHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIds,
                                  const ProcessInfo&             rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentCondition, rRHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(rEquationIds, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Condition&                     rCurrentCondition,
                                  LocalSystemMatrixType&         LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Element&                       rCurrentElement,
                                  LocalSystemMatrixType&         LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        // Empty function because the derivatives are updated outside of the function since it needs more input. And this function is required to be defined 
    }

protected:


    void UpdateVectorFirstTimeDerivative(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

            const array_1d<double, 3> updated_first_derivative =
                rNode.FastGetSolutionStepValue(r_second_order_vector_variable.first_time_derivative, 1) +
                (1.0 - GetGamma()) * this->GetDeltaTime() *
                rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 1) +
                GetGamma() * this->GetDeltaTime() *
                rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 0);

            NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                rNode, r_second_order_vector_variable.first_time_derivative, updated_first_derivative);
        }
    }

    void UpdateVectorSecondTimeDerivative(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

            const array_1d<double, 3> updated_second_time_derivative =
                ((rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 0) -
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 1)) -
                    this->GetDeltaTime() * rNode.FastGetSolutionStepValue(
                        r_second_order_vector_variable.first_time_derivative, 1) -
                    (0.5 - GetBeta()) * this->GetDeltaTime() * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 1)) /
                (GetBeta() * this->GetDeltaTime() * this->GetDeltaTime());

            NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                rNode, r_second_order_vector_variable.second_time_derivative, updated_second_time_derivative);
        }
    }

    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo&     CurrentProcessInfo)
    {
        KRATOS_TRY

        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) +=
                (1.0 / (this->GetBeta() * this->GetDeltaTime() * this->GetDeltaTime())) * M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += (this->GetGamma() / (this->GetBeta() * this->GetDeltaTime())) * C;

        KRATOS_CATCH("")
    }

    void AddDynamicsToRHS(Condition&             rCurrentCondition,
                          LocalSystemVectorType& RHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo&     CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        // adding inertia contribution
        if (M.size1() != 0) {
            rCurrentCondition.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        // adding damping contribution
        if (C.size1() != 0) {
            rCurrentCondition.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }

        KRATOS_CATCH("")
    }

    void AddDynamicsToRHS(Element&               rCurrentElement,
                          LocalSystemVectorType& RHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo&     CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        // adding inertia contribution
        if (M.size1() != 0) {
            rCurrentElement.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        // adding damping contribution
        if (C.size1() != 0) {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }

        KRATOS_CATCH("")
    }

    inline void PredictVariablesDerivatives(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // make sure the timestep size is correct
        this->SetTimeFactors(rModelPart);

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            
            for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
                if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

                // firstly get current values of first and second derivative
                const array_1d<double, 3> first_derivative_array =
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.first_time_derivative, 0);

                const array_1d<double, 3> second_derivative_array =
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 0);

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


    void GetFirstAndSecondDerivativeVector(TSystemVectorType& rFirstDerivativeVector,
        TSystemVectorType& rSecondDerivativeVector,
        ModelPart& rModelPart,
        IndexType          i)
    {
        block_for_each(rModelPart.Nodes(),
            [&rFirstDerivativeVector, &rSecondDerivativeVector, i, this](Node& rNode) {
                if (rNode.IsActive()) {
                    GetDerivativesForVariable(DISPLACEMENT_X, rNode, rFirstDerivativeVector,
                        rSecondDerivativeVector, i);
                    GetDerivativesForVariable(DISPLACEMENT_Y, rNode, rFirstDerivativeVector,
                        rSecondDerivativeVector, i);

                    const std::vector<const Variable<double>*> optional_variables = {
                        &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z };

                    for (const auto p_variable : optional_variables) {
                        GetDerivativesForOptionalVariable(*p_variable, rNode, rFirstDerivativeVector,
                            rSecondDerivativeVector, i);
                    }
                }
            });
    }

    void SetFirstAndSecondDerivativeVector(TSystemVectorType& rFirstDerivativeVector,
        TSystemVectorType& rSecondDerivativeVector,
        ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [&rFirstDerivativeVector, &rSecondDerivativeVector, this](Node& rNode) {
            if (rNode.IsActive()) {
                SetDerivativesForVariable(DISPLACEMENT_X, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
                SetDerivativesForVariable(DISPLACEMENT_Y, rNode, rFirstDerivativeVector, rSecondDerivativeVector);

                const std::vector<const Variable<double>*> optional_variables = {
                    &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z };

                for (const auto p_variable : optional_variables) {
                    SetDerivativesForOptionalVariable(*p_variable, rNode, rFirstDerivativeVector,
                        rSecondDerivativeVector);
                }
            }
            });
    }

    void GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
        const Node& rNode,
        TSystemVectorType& rFirstDerivativeVector,
        TSystemVectorType& rSecondDerivativeVector,
        IndexType               i) const
    {
        if (rNode.HasDofFor(rVariable)) {
            GetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector, i);
        }
    }

    void SetDerivativesForOptionalVariable(const Variable<double>& rVariable,
        Node& rNode,
        const TSystemVectorType& rFirstDerivativeVector,
        const TSystemVectorType& rSecondDerivativeVector)
    {
        if (rNode.HasDofFor(rVariable)) {
            SetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
        }
    }

    void GetDerivativesForVariable(const Variable<double>& rVariable,
        const Node& rNode,
        TSystemVectorType& rFirstDerivativeVector,
        TSystemVectorType& rSecondDerivativeVector,
        IndexType               i) const
    {
        const auto& r_first_derivative = rVariable.GetTimeDerivative();
        const auto& r_second_derivative = r_first_derivative.GetTimeDerivative();

        const auto equation_id = rNode.GetDof(rVariable).EquationId();
        rFirstDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(r_first_derivative, i);
        rSecondDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(r_second_derivative, i);
    }

    void SetDerivativesForVariable(const Variable<double>& rVariable,
        Node& rNode,
        const TSystemVectorType& rFirstDerivativeVector,
        const TSystemVectorType& rSecondDerivativeVector)
    {
        const auto& r_first_derivative = rVariable.GetTimeDerivative();
        const auto& r_second_derivative = r_first_derivative.GetTimeDerivative();

        const auto equation_id = rNode.GetDof(rVariable).EquationId();
        rNode.FastGetSolutionStepValue(r_first_derivative) = rFirstDerivativeVector[equation_id];
        rNode.FastGetSolutionStepValue(r_second_derivative) = rSecondDerivativeVector[equation_id];
    }


private:
    std::vector<Matrix> mMassMatrix;
    std::vector<Vector> mAccelerationVector;
    std::vector<Matrix> mDampingMatrix;
    std::vector<Vector> mVelocityVector;

}; // Class NewmarkDynamicUPwScheme

} // namespace Kratos
