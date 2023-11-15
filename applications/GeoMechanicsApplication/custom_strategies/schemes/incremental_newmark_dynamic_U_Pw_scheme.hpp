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
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class IncrementalNewmarkDynamicUPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(IncrementalNewmarkDynamicUPwScheme);

    using BaseType              = Scheme<TSparseSpace,TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mBeta;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mGamma;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mTheta;

    IncrementalNewmarkDynamicUPwScheme(double beta, double gamma, double theta)
        : NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(beta, gamma, theta)
    {
        //Allocate auxiliary memory
        int NumThreads = ParallelUtilities::GetNumThreads();
        mMassMatrix.resize(NumThreads);
        mAccelerationVector.resize(NumThreads);
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
        mPreviousRhs.resize(NumThreads);
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        // initialize element as it is required to initialize solution step and non linear iteration
        BaseType::InitializeElements(rModelPart);

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        // Initialize solutionstep and non-linear iteration here as LHS is required to be constant throughout the computation
        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive) rElement.InitializeSolutionStep(rCurrentProcessInfo);
            });

        // initialize non linear iteration, here as LHS is requried to be constant 
        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive) rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            });


        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

            SetTimeFactors(rModelPart);

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive) rCondition.InitializeSolutionStep(rCurrentProcessInfo);
            });

        KRATOS_CATCH("")
    }

    void InitializeNonLinIteration(ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

            const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();


        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive) rCondition.InitializeNonLinearIteration(rCurrentProcessInfo);
            });

        KRATOS_CATCH("")
    }


    void FinalizeNonLinIteration(ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

            const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        // finalize non linear interation for element in order to calculate stresses within element
        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive) rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive) rCondition.FinalizeNonLinearIteration(rCurrentProcessInfo);
            });

        KRATOS_CATCH("")
    }

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
		// no prediction required
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        // here dynamics are only added to lhs, dynamics are added to the right hand side within the builder_and_solver
        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        // no rhs contribution from elements
        this->CalculateLHSContribution(
            rCurrentElement,
            LHS_Contribution,
            EquationId, CurrentProcessInfo);

        if (RHS_Contribution.size() != EquationId.size())
            RHS_Contribution.resize(EquationId.size(), false);
        
        TSparseSpace::SetToZero(RHS_Contribution);



        // int thread = OpenMPUtils::ThisThread();
        //
        // rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        //
        // TSparseSpace::SetToZero(RHS_Contribution);
        //
        // rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);
        //
        // rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);
        //
        // this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);
        //
        // //this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);
        //
        // rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);


        KRATOS_CATCH( "" )
    }

    void CalculateRHSContribution(
        Element &rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
       
            KRATOS_ERROR << "Build_and_solver should not request rhs contribution from IncrementalNewmarkDynamicUPwScheme" << std::endl;;
    
        // int thread = OpenMPUtils::ThisThread();
        //
        //
        // rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);
        //
        // rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);
        //
        // this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);
        //
        // rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIds,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        //KRATOS_TRY
        KRATOS_ERROR << "Build_and_solver should not request rhs contribution from IncrementalNewmarkDynamicUPwScheme" << std::endl;;
        // int thread = OpenMPUtils::ThisThread();
        //
        // rCurrentCondition.CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
        //
        // rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);
        //
        // rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);
        //
        // this->AddDynamicsToRHS(rCurrentCondition, rRHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(rEquationIds, rCurrentProcessInfo);

        //KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
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

    void CalculateLHSContribution(
        Element &rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

protected:
    /// Member Variables
    std::vector< Matrix > mMassMatrix;
    std::vector< Vector > mAccelerationVector;
    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;
    std::vector< Vector > mPreviousRhs;

    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) += (1.0/(mBeta*mDeltaTime*mDeltaTime))*M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += (mGamma/(mBeta*mDeltaTime))*C;

        KRATOS_CATCH( "" )
    }

    // void AddDynamicsToRHS(Condition& rCurrentCondition,
    //     LocalSystemVectorType& RHS_Contribution,
    //     LocalSystemMatrixType& M,
    //     LocalSystemMatrixType& C,
    //     const ProcessInfo& CurrentProcessInfo)
    // {
    //     KRATOS_TRY
    //
    //         int thread = OpenMPUtils::ThisThread();
    //
    //     if (M.size1() != 0 || C.size1() != 0)
    //     {
    //         rCurrentCondition.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
    //         rCurrentCondition.GetFirstDerivativesVector(mVelocityVector[thread], 0);
    //
    //         //adding inertia contribution
    //         if (M.size1() != 0)
    //         {
    //
    //             auto m_part = mVelocityVector[thread] * (1 / (mBeta * mDeltaTime)) + mAccelerationVector[thread] * (1 / (2 * mBeta));
    //
    //             noalias(RHS_Contribution) += prod(M, m_part);
    //         }
    //
    //         //adding damping contribution
    //         if (C.size1() != 0)
    //         {
    //             auto c_part = mVelocityVector[thread] * (mGamma / mBeta) + mAccelerationVector[thread] * (mDeltaTime * (mGamma / (2 * mBeta) - 1));
    //
    //             noalias(RHS_Contribution) += prod(C, c_part);
    //         }
    //     }
    //
    //     KRATOS_CATCH("")
    // }
    //
    // void AddDynamicsToRHS(Element &rCurrentElement,
    //                       LocalSystemVectorType& RHS_Contribution,
    //                       LocalSystemMatrixType& M,
    //                       LocalSystemMatrixType& C,
    //                       const ProcessInfo& CurrentProcessInfo)
    // {
    //     KRATOS_TRY
    //
    //     int thread = OpenMPUtils::ThisThread();
    //
    //     if (M.size1() != 0 || C.size1() != 0)
    //     {
    //         rCurrentElement.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
    //         rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);
    //
    //         //adding inertia contribution
    //         if (M.size1() != 0)
    //         {
    //
    //             auto m_part = mVelocityVector[thread] * (1 / (mBeta * mDeltaTime)) + mAccelerationVector[thread] * (1 / (2 * mBeta));
    //
    //             noalias(RHS_Contribution) += prod(M, mAccelerationVector[thread]);
    //         }
    //
    //         //adding damping contribution
    //         if (C.size1() != 0)
    //         {
    //             auto c_part = mVelocityVector[thread] * (mGamma / mBeta) + mAccelerationVector[thread] * (mDeltaTime * (mGamma / (2 * mBeta) - 1));
    //
    //             noalias(RHS_Contribution) += prod(C, mVelocityVector[thread]);
    //         }
    //     }
    //
    //     KRATOS_CATCH( "" )
    // }


    void Update(ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

            int NumThreads = ParallelUtilities::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k + 1];

            //Update Displacement and Pressure (DOFs)
            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof) {
                if (itDof->IsFree()) itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx, itDof->EquationId());
            }
        }

       // this->UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }


    void FinalizeSolutionStep(ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        NewmarkQuasistaticUPwScheme::FinalizeSolutionStep(rModelPart, A, Dx, b);

        this->UpdateVariablesDerivatives(rModelPart);
    }


    virtual inline void UpdateVariablesDerivatives(ModelPart& rModelPart)
    {
        KRATOS_TRY

            //Update Acceleration, Velocity and DtPressure
            block_for_each(rModelPart.Nodes(), [&](Node& rNode) {

	            const auto delta_u = rNode.FastGetSolutionStepValue(DISPLACEMENT) - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);

	            const auto delta_a = (delta_u * (1 / (mBeta * mDeltaTime * mDeltaTime))
	                - rNode.FastGetSolutionStepValue(VELOCITY) * (1 / (mBeta * mDeltaTime))
	                - rNode.FastGetSolutionStepValue(ACCELERATION) * (1 / (2 * mBeta)));

	            const auto delta_v = delta_u * (mGamma / (mBeta * mDeltaTime))
	                - rNode.FastGetSolutionStepValue(VELOCITY) * (mGamma / mBeta)
	                + rNode.FastGetSolutionStepValue(ACCELERATION) * (mDeltaTime * (1 - mGamma / (2 * mBeta)));


	            rNode.FastGetSolutionStepValue(ACCELERATION) += delta_a;
	            rNode.FastGetSolutionStepValue(VELOCITY) += delta_v;

                });

        KRATOS_CATCH("")
    }


}; // Class IncrementalNewmarkDynamicUPwScheme

}