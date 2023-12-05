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
        int NumThreads = ParallelUtilities::GetNumThreads();
        mMassMatrix.resize(NumThreads);
        mAccelerationVector.resize(NumThreads);
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    void Predict(ModelPart& rModelPart,
                 DofsArrayType& rDofSet,
                 TSystemMatrixType& A,
                 TSystemVectorType& Dx,
                 TSystemVectorType& b) override
    {
        KRATOS_TRY

        // Predict Displacements on free nodes and update Acceleration, Velocity and DtPressure
        block_for_each(
            rModelPart.Nodes(),
            [&](Node& rNode)
            {
                // refactor such that all mathematical formulas appear once only.
                if (rNode.IsFixed(ACCELERATION_X))
                {
                    const double& PreviousDisplacement =
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_X, 1);
                    const double& PreviousVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_X, 1);
                    const double& PreviousAcceleration =
                        rNode.FastGetSolutionStepValue(ACCELERATION_X, 1);
                    const double& CurrentAcceleration =
                        rNode.FastGetSolutionStepValue(ACCELERATION_X);

                    rNode.FastGetSolutionStepValue(DISPLACEMENT_X) =
                        PreviousDisplacement + this->GetDeltaTime() * PreviousVelocity +
                        this->GetDeltaTime() * this->GetDeltaTime() *
                            ((0.5 - mBeta) * PreviousAcceleration + mBeta * CurrentAcceleration);
                }
                else if (rNode.IsFixed(VELOCITY_X))
                {
                    const double& PreviousDisplacement =
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_X, 1);
                    const double& PreviousVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_X, 1);
                    const double& CurrentVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_X);
                    rNode.FastGetSolutionStepValue(DISPLACEMENT_X) =
                        PreviousDisplacement +
                        this->GetDeltaTime() *
                            ((mBeta / mGamma) * (CurrentVelocity - PreviousVelocity) +
                             PreviousVelocity);
                }
                else if (!rNode.IsFixed(DISPLACEMENT_X))
                {
                    const double& PreviousDisplacement =
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_X, 1);
                    const double& PreviousVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_X, 1);
                    const double& PreviousAcceleration =
                        rNode.FastGetSolutionStepValue(ACCELERATION_X, 1);

                    rNode.FastGetSolutionStepValue(DISPLACEMENT_X) =
                        PreviousDisplacement + this->GetDeltaTime() * PreviousVelocity +
                        0.5 * this->GetDeltaTime() * this->GetDeltaTime() * PreviousAcceleration;
                }

                if (rNode.IsFixed(ACCELERATION_Y))
                {
                    const double& PreviousDisplacement =
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 1);
                    const double& PreviousVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_Y, 1);
                    const double& PreviousAcceleration =
                        rNode.FastGetSolutionStepValue(ACCELERATION_Y, 1);
                    const double& CurrentAcceleration =
                        rNode.FastGetSolutionStepValue(ACCELERATION_Y);

                    rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) =
                        PreviousDisplacement + this->GetDeltaTime() * PreviousVelocity +
                        this->GetDeltaTime() * this->GetDeltaTime() *
                            ((0.5 - mBeta) * PreviousAcceleration + mBeta * CurrentAcceleration);
                }
                else if (rNode.IsFixed(VELOCITY_Y))
                {
                    const double& PreviousDisplacement =
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 1);
                    const double& PreviousVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_Y, 1);
                    const double& CurrentVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_Y);
                    rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) =
                        PreviousDisplacement +
                        this->GetDeltaTime() *
                            ((mBeta / mGamma) * (CurrentVelocity - PreviousVelocity) +
                             PreviousVelocity);
                }
                else if (!rNode.IsFixed(DISPLACEMENT_Y))
                {
                    const double& PreviousDisplacement =
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 1);
                    const double& PreviousVelocity =
                        rNode.FastGetSolutionStepValue(VELOCITY_Y, 1);
                    const double& PreviousAcceleration =
                        rNode.FastGetSolutionStepValue(ACCELERATION_Y, 1);

                    rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) =
                        PreviousDisplacement + this->GetDeltaTime() * PreviousVelocity +
                        0.5 * this->GetDeltaTime() * this->GetDeltaTime() * PreviousAcceleration;
                }

                // For 3D cases
                if (rNode.HasDofFor(DISPLACEMENT_Z))
                {
                    if (rNode.IsFixed(ACCELERATION_Z))
                    {
                        const double& PreviousDisplacement =
                            rNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 1);
                        const double& PreviousVelocity =
                            rNode.FastGetSolutionStepValue(VELOCITY_Z, 1);
                        const double& PreviousAcceleration =
                            rNode.FastGetSolutionStepValue(ACCELERATION_Z, 1);
                        const double& CurrentAcceleration =
                            rNode.FastGetSolutionStepValue(ACCELERATION_Z);

                        rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) =
                            PreviousDisplacement + this->GetDeltaTime() * PreviousVelocity +
                            this->GetDeltaTime() * this->GetDeltaTime() *
                                ((0.5 - mBeta) * PreviousAcceleration + mBeta * CurrentAcceleration);
                    }
                    else if (rNode.IsFixed(VELOCITY_Z))
                    {
                        const double& PreviousDisplacement =
                            rNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 1);
                        const double& PreviousVelocity =
                            rNode.FastGetSolutionStepValue(VELOCITY_Z, 1);
                        const double& CurrentVelocity =
                            rNode.FastGetSolutionStepValue(VELOCITY_Z);
                        rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) =
                            PreviousDisplacement +
                            this->GetDeltaTime() *
                                ((mBeta / mGamma) * (CurrentVelocity - PreviousVelocity) +
                                 PreviousVelocity);
                    }
                    else if (!rNode.IsFixed(DISPLACEMENT_Z))
                    {
                        const double& PreviousDisplacement =
                            rNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 1);
                        const double& PreviousVelocity =
                            rNode.FastGetSolutionStepValue(VELOCITY_Z, 1);
                        const double& PreviousAcceleration =
                            rNode.FastGetSolutionStepValue(ACCELERATION_Z, 1);

                        rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) =
                            PreviousDisplacement + this->GetDeltaTime() * PreviousVelocity +
                            0.5 * this->GetDeltaTime() * this->GetDeltaTime() * PreviousAcceleration;
                    }
                }

                noalias(rNode.FastGetSolutionStepValue(ACCELERATION)) =
                    (1.0 / (mBeta * this->GetDeltaTime() * this->GetDeltaTime())) *
                    ((rNode.FastGetSolutionStepValue(DISPLACEMENT) -
                      rNode.FastGetSolutionStepValue(DISPLACEMENT, 1)) -
                     this->GetDeltaTime() * rNode.FastGetSolutionStepValue(VELOCITY, 1) -
                     (0.5 - mBeta) * this->GetDeltaTime() * this->GetDeltaTime() *
                         rNode.FastGetSolutionStepValue(ACCELERATION, 1));

                noalias(rNode.FastGetSolutionStepValue(VELOCITY)) =
                    rNode.FastGetSolutionStepValue(VELOCITY, 1) +
                    (1.0 - mGamma) * this->GetDeltaTime() *
                        rNode.FastGetSolutionStepValue(ACCELERATION, 1) +
                    mGamma * this->GetDeltaTime() *
                        rNode.FastGetSolutionStepValue(ACCELERATION);

                this->UpdateScalarTimeDerivative(rNode, WATER_PRESSURE, DT_WATER_PRESSURE);
            });

        KRATOS_CATCH("")
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