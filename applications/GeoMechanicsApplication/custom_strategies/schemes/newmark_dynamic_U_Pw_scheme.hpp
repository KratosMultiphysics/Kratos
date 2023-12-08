// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class NewmarkDynamicUPwScheme
    : public NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkDynamicUPwScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>::mBeta;
    using NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>::mGamma;

    NewmarkDynamicUPwScheme(double beta, double gamma, double theta)
        : NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>(beta, gamma, theta)
    {
        // Allocate auxiliary memory
        int num_threads = ParallelUtilities::GetNumThreads();
        mMassMatrix.resize(num_threads);
        mAccelerationVector.resize(num_threads);
        mDampingMatrix.resize(num_threads);
        mVelocityVector.resize(num_threads);
    }

    void Predict(ModelPart& rModelPart,
                 DofsArrayType& rDofSet,
                 TSystemMatrixType& A,
                 TSystemVectorType& Dx,
                 TSystemVectorType& b) override
    {
        KRATOS_TRY

        PredictVariables(rModelPart);
        // Update (Angular) Acceleration, (Angular) Velocity and DtPressure
        this->UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }
    void PredictVariables(const ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(),
                       [this](Node& rNode) { PredictVariablesForNode(rNode); });
    }

    void PredictVariablesForNode(Node& rNode)
    {
        for (const auto& variable_derivative : this->GetVariableDerivatives())
        {
            if (!rNode.SolutionStepsDataHas(variable_derivative.instance))
                continue;
            PredictVariableForNode(rNode, variable_derivative);
        }
    }

    void PredictVariableForNode(Node& rNode, const VariableWithTimeDerivatives& rVariableWithDerivatives)
    {
        std::vector<std::string> components = {"X", "Y"};
        if (rNode.HasDofFor(this->GetComponentFromVectorVariable(
                rVariableWithDerivatives.instance, "Z")))
            components.emplace_back("Z");

        for (const auto& component : components)
        {
            const auto& instance_component = this->GetComponentFromVectorVariable(
                rVariableWithDerivatives.instance, component);
            const auto& first_time_derivative_component = this->GetComponentFromVectorVariable(
                rVariableWithDerivatives.first_time_derivative, component);
            const auto& second_time_derivative_component = this->GetComponentFromVectorVariable(
                rVariableWithDerivatives.second_time_derivative, component);

            const double previous_variable =
                rNode.FastGetSolutionStepValue(instance_component, 1);
            const double current_first_time_derivative =
                rNode.FastGetSolutionStepValue(first_time_derivative_component, 0);
            const double previous_first_time_derivative =
                rNode.FastGetSolutionStepValue(first_time_derivative_component, 1);
            const double current_second_time_derivative =
                rNode.FastGetSolutionStepValue(second_time_derivative_component, 0);
            const double previous_second_time_derivative =
                rNode.FastGetSolutionStepValue(second_time_derivative_component, 1);

            if (rNode.IsFixed(second_time_derivative_component))
            {
                rNode.FastGetSolutionStepValue(instance_component) =
                    previous_variable + this->GetDeltaTime() * previous_first_time_derivative +
                    this->GetDeltaTime() * this->GetDeltaTime() *
                        ((0.5 - mBeta) * previous_second_time_derivative +
                         mBeta * current_second_time_derivative);
            }
            else if (rNode.IsFixed(first_time_derivative_component))
            {
                rNode.FastGetSolutionStepValue(instance_component) =
                    previous_variable +
                    this->GetDeltaTime() *
                        ((mBeta / mGamma) * (current_first_time_derivative - previous_first_time_derivative) +
                         previous_first_time_derivative);
            }
            else if (!rNode.IsFixed(instance_component))
            {
                rNode.FastGetSolutionStepValue(instance_component) =
                    previous_variable + this->GetDeltaTime() * previous_first_time_derivative +
                    0.5 * this->GetDeltaTime() * this->GetDeltaTime() * previous_second_time_derivative;
            }
        }
    }

    void CalculateSystemContributions(Condition& rCurrentCondition,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(Element& rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution,
                                             CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Element& rCurrentElement,
                                  LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
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

    void CalculateRHSContribution(Condition& rCurrentCondition,
                                  LocalSystemVectorType& rRHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIds,
                                  const ProcessInfo& rCurrentProcessInfo) override
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

    void CalculateLHSContribution(Condition& rCurrentCondition,
                                  LocalSystemMatrixType& LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Element& rCurrentElement,
                                  LocalSystemMatrixType& LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

protected:
    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) +=
                (1.0 / (mBeta * this->GetDeltaTime() * this->GetDeltaTime())) * M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += (mGamma / (mBeta * this->GetDeltaTime())) * C;

        KRATOS_CATCH("")
    }

    void AddDynamicsToRHS(Condition& rCurrentCondition,
                          LocalSystemVectorType& RHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        // adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        // adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentCondition.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }

        KRATOS_CATCH("")
    }

    void AddDynamicsToRHS(Element& rCurrentElement,
                          LocalSystemVectorType& RHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        // adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        // adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }

        KRATOS_CATCH("")
    }

private:
    std::vector<Matrix> mMassMatrix;
    std::vector<Vector> mAccelerationVector;
    std::vector<Matrix> mDampingMatrix;
    std::vector<Vector> mVelocityVector;
}; // Class NewmarkDynamicUPwScheme

} // namespace Kratos