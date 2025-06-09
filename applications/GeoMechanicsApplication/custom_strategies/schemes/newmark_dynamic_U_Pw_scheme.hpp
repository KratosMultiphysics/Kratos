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
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class NewmarkDynamicUPwScheme : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkDynamicUPwScheme);

    using BaseType              = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    NewmarkDynamicUPwScheme(double beta, double gamma, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)},
              beta,
              gamma,
              theta)
    {
        // Allocate auxiliary memory
        int num_threads = ParallelUtilities::GetNumThreads();
        mMassMatrix.resize(num_threads);
        mAccelerationVector.resize(num_threads);
        mDampingMatrix.resize(num_threads);
        mVelocityVector.resize(num_threads);
    }

    void Predict(ModelPart& rModelPart, DofsArrayType&, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        PredictVariables(rModelPart);
        // Update (Angular) Acceleration, (Angular) Velocity and DtPressure
        this->UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }

    void PredictVariables(const ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [this](Node& rNode) { PredictVariablesForNode(rNode); });
    }

    void PredictVariablesForNode(Node& rNode)
    {
        for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;
            PredictVariableForNode(rNode, r_second_order_vector_variable);
        }
    }

    void PredictVariableForNode(Node& rNode, const SecondOrderVectorVariable& rSecondOrderVariables)
    {
        const std::vector<std::string> components = {"X", "Y", "Z"};

        for (const auto& component : components) {
            const auto& instance_component = VariablesUtilities::GetComponentFromVectorVariable(
                rSecondOrderVariables.instance.Name(), component);

            if (!rNode.HasDofFor(instance_component)) continue;

            const auto& first_time_derivative_component = VariablesUtilities::GetComponentFromVectorVariable(
                rSecondOrderVariables.first_time_derivative.Name(), component);
            const auto& second_time_derivative_component = VariablesUtilities::GetComponentFromVectorVariable(
                rSecondOrderVariables.second_time_derivative.Name(), component);

            const double previous_variable = rNode.FastGetSolutionStepValue(instance_component, 1);
            const double current_first_time_derivative =
                rNode.FastGetSolutionStepValue(first_time_derivative_component, 0);
            const double previous_first_time_derivative =
                rNode.FastGetSolutionStepValue(first_time_derivative_component, 1);
            const double current_second_time_derivative =
                rNode.FastGetSolutionStepValue(second_time_derivative_component, 0);
            const double previous_second_time_derivative =
                rNode.FastGetSolutionStepValue(second_time_derivative_component, 1);
            if (rNode.IsFixed(second_time_derivative_component)) {
                rNode.FastGetSolutionStepValue(instance_component) =
                    previous_variable + this->GetDeltaTime() * previous_first_time_derivative +
                    this->GetDeltaTime() * this->GetDeltaTime() *
                        ((0.5 - this->GetBeta()) * previous_second_time_derivative +
                         this->GetBeta() * current_second_time_derivative);
            } else if (rNode.IsFixed(first_time_derivative_component)) {
                rNode.FastGetSolutionStepValue(instance_component) =
                    previous_variable +
                    this->GetDeltaTime() * ((this->GetBeta() / this->GetGamma()) *
                                                (current_first_time_derivative - previous_first_time_derivative) +
                                            previous_first_time_derivative);
            } else if (!rNode.IsFixed(instance_component)) {
                rNode.FastGetSolutionStepValue(instance_component) =
                    previous_variable + this->GetDeltaTime() * previous_first_time_derivative +
                    0.5 * this->GetDeltaTime() * this->GetDeltaTime() * previous_second_time_derivative;
            }
        }
    }

    void CalculateSystemContributions(Condition&                     rCurrentCondition,
                                      LocalSystemMatrixType&         rLHS_Contribution,
                                      LocalSystemVectorType&         rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationIds,
                                      const ProcessInfo&             rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToLHS(rLHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentCondition, rRHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(rEquationIds, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(Element&                       rCurrentElement,
                                      LocalSystemMatrixType&         rLHS_Contribution,
                                      LocalSystemVectorType&         rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationIds,
                                      const ProcessInfo&             rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToLHS(rLHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentElement, rRHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(rEquationIds, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Element&                       rCurrentElement,
                                  LocalSystemVectorType&         rRHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIds,
                                  const ProcessInfo&             rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentElement, rRHS_Contribution, mMassMatrix[thread],
                               mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(rEquationIds, rCurrentProcessInfo);

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
                                  LocalSystemMatrixType&         rLHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIds,
                                  const ProcessInfo&             rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToLHS(rLHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(rEquationIds, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Element&                       rCurrentElement,
                                  LocalSystemMatrixType&         rLHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIds,
                                  const ProcessInfo&             rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToLHS(rLHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(rEquationIds, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

protected:
    void AddDynamicsToLHS(LocalSystemMatrixType& rLHS_Contribution,
                          LocalSystemMatrixType& rM,
                          LocalSystemMatrixType& rC,
                          const ProcessInfo&     rCurrentProcessInfo)
    {
        KRATOS_TRY

        // adding mass contribution
        if (rM.size1() != 0)
            noalias(rLHS_Contribution) +=
                (1.0 / (this->GetBeta() * this->GetDeltaTime() * this->GetDeltaTime())) * rM;

        // adding damping contribution
        if (rC.size1() != 0)
            noalias(rLHS_Contribution) += (this->GetGamma() / (this->GetBeta() * this->GetDeltaTime())) * rC;

        KRATOS_CATCH("")
    }

    void AddDynamicsToRHS(Condition&             rCurrentCondition,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& rM,
                          LocalSystemMatrixType& rC,
                          const ProcessInfo&     rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        // adding inertia contribution
        if (rM.size1() != 0) {
            rCurrentCondition.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, mAccelerationVector[thread]);
        }

        // adding damping contribution
        if (rC.size1() != 0) {
            rCurrentCondition.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(rRHS_Contribution) -= prod(rC, mVelocityVector[thread]);
        }

        KRATOS_CATCH("")
    }

    void AddDynamicsToRHS(Element&               rCurrentElement,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& rM,
                          LocalSystemMatrixType& rC,
                          const ProcessInfo&     rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        // adding inertia contribution
        if (rM.size1() != 0) {
            rCurrentElement.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, mAccelerationVector[thread]);
        }

        // adding damping contribution
        if (rC.size1() != 0) {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(rRHS_Contribution) -= prod(rC, mVelocityVector[thread]);
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