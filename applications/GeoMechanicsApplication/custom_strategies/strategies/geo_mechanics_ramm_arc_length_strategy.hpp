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

// Project includes
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoMechanicsRammArcLengthStrategy
    : public GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicsRammArcLengthStrategy);

    using BaseType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using GrandMotherType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using MotherType = GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;
    using SparseSpaceType          = TSparseSpace;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;
    using TSystemVectorPointerType = typename BaseType::TSystemVectorPointerType;
    using GrandMotherType::mCalculateReactionsFlag;
    using GrandMotherType::mInitializeWasPerformed;
    using GrandMotherType::mMaxIterationNumber;
    using GrandMotherType::mpA; // Tangent matrix
    using GrandMotherType::mpb; // Residual vector of iteration i
    using GrandMotherType::mpBuilderAndSolver;
    using GrandMotherType::mpConvergenceCriteria;
    using GrandMotherType::mpDx; // Delta x of iteration i
    using GrandMotherType::mpScheme;
    using GrandMotherType::mReformDofSetAtEachStep;
    using MotherType::mSubModelPartList;
    using MotherType::mVariableNames;

    GeoMechanicsRammArcLengthStrategy(ModelPart&                      model_part,
                                      typename TSchemeType::Pointer   pScheme,
                                      typename TLinearSolver::Pointer pNewLinearSolver,
                                      typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                      typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                      Parameters&                             rParameters,
                                      int                                     MaxIterations = 30,
                                      bool CalculateReactions                               = false,
                                      bool ReformDofSetAtEachStep                           = false,
                                      bool MoveMeshFlag                                     = false)
        : GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              model_part,
              pScheme,
              pNewLinearSolver,
              pNewConvergenceCriteria,
              pNewBuilderAndSolver,
              rParameters,
              MaxIterations,
              CalculateReactions,
              ReformDofSetAtEachStep,
              MoveMeshFlag)
    {
        mDesiredIterations = rParameters["desired_iterations"].GetInt();
        mMaxRadiusFactor   = rParameters["max_radius_factor"].GetDouble();
        mMinRadiusFactor   = rParameters["min_radius_factor"].GetDouble();
    }

    void Initialize() override
    {
        KRATOS_TRY

        if (!mInitializeWasPerformed) {
            MotherType::Initialize();

            if (!mInitializeArcLengthWasPerformed) {
                // set up the system
                if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false) {
                    // setting up the list of the DOFs to be solved
                    mpBuilderAndSolver->SetUpDofSet(mpScheme, BaseType::GetModelPart());

                    // shaping correctly the system
                    mpBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
                }

                // Compute initial radius (mRadius_0)
                mpBuilderAndSolver->ResizeAndInitializeVectors(mpScheme, mpA, mpDx, mpb,
                                                               BaseType::GetModelPart());
                TSystemMatrixType& mA  = *mpA;
                TSystemVectorType& mDx = *mpDx;
                TSystemVectorType& mb  = *mpb;
                TSparseSpace::SetToZero(mA);
                TSparseSpace::SetToZero(mDx);
                TSparseSpace::SetToZero(mb);

                mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

                mRadius_0 = TSparseSpace::TwoNorm(mDx);
                mRadius   = mRadius_0;

                // Compute vector of reference external force (mf)
                this->InitializeSystemVector(mpf);
                TSystemVectorType& mf = *mpf;
                TSparseSpace::SetToZero(mf);

                mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mf);

                // Initialize the loading factor Lambda
                mLambda     = 0.0;
                mLambda_old = 1.0;

                // Initialize Norm of solution
                mNormxEquilibrium = 0.0;

                mInitializeArcLengthWasPerformed = true;

                KRATOS_INFO("Ramm's Arc Length Strategy") << "Strategy Initialized" << std::endl;
            }
        }

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        GrandMotherType::InitializeSolutionStep();

        this->SaveInitializeSystemVector(mpf);
        this->InitializeSystemVector(mpDxf);
        this->InitializeSystemVector(mpDxb);
        this->InitializeSystemVector(mpDxPred);
        this->InitializeSystemVector(mpDxStep);

        KRATOS_CATCH("")
    }

    bool SolveSolutionStep() override
    {
        // ********** Prediction phase **********
        KRATOS_INFO("Ramm's Arc Length Strategy")
            << "ARC-LENGTH RADIUS: " << mRadius / mRadius_0 << " X initial radius" << std::endl;

        // Initialize variables
        DofsArrayType&     rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA      = *mpA;
        TSystemVectorType& mDx     = *mpDx;
        TSystemVectorType& mb      = *mpb;
        TSystemVectorType& mf      = *mpf;
        TSystemVectorType& mDxb    = *mpDxb;
        TSystemVectorType& mDxf    = *mpDxf;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;

        // initializing the parameters of the iteration loop
        double       NormDx;
        unsigned int iteration_number                                  = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        bool is_converged =
            mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDxf);

        // Note: This is not so efficient, but I want to solve mA*mDxf=mf without losing mf
        this->BuildWithDirichlet(mA, mDxf, mb);
        noalias(mb) = mf;
        mpBuilderAndSolver->SystemSolve(mA, mDxf, mb);

        // update results
        double DLambda = mRadius / TSparseSpace::TwoNorm(mDxf);
        mDLambdaStep   = DLambda;
        mLambda += DLambda;
        noalias(mDxPred) = DLambda * mDxf;
        noalias(mDxStep) = mDxPred;
        this->Update(rDofSet, mA, mDxPred, mb);

        // move the mesh if needed
        if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();

        // ********** Correction phase (iteration cycle) **********
        if (is_converged) {
            mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDxf, mb);
            if (mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(mb);
                mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
            }
            is_converged =
                mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDxf, mb);
        }

        while (!is_converged && iteration_number++ < mMaxIterationNumber) {
            // setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxf);

            // Note: This is not so efficient, but I want to solve mA*mDxf=mf without losing mf
            this->BuildWithDirichlet(mA, mDxf, mb);
            noalias(mb) = mf;
            mpBuilderAndSolver->SystemSolve(mA, mDxf, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxb);

            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDxb, mb);

            DLambda = -TSparseSpace::Dot(mDxPred, mDxb) / TSparseSpace::Dot(mDxPred, mDxf);

            noalias(mDx) = mDxb + DLambda * mDxf;

            // Check solution before update
            if (mNormxEquilibrium > 1.0e-10) {
                NormDx = TSparseSpace::TwoNorm(mDx);

                if ((NormDx / mNormxEquilibrium) > 1.0e3 ||
                    (std::abs(DLambda) / std::abs(mLambda - mDLambdaStep)) > 1.0e3) {
                    is_converged = false;
                    break;
                }
            }

            // update results
            mDLambdaStep += DLambda;
            mLambda += DLambda;
            noalias(mDxStep) += mDx;
            this->Update(rDofSet, mA, mDx, mb);

            // move the mesh if needed
            if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();

            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            // *** Check Convergence ***
            if (is_converged) {
                if (mpConvergenceCriteria->GetActualizeRHSflag()) {
                    TSparseSpace::SetToZero(mb);
                    mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
                }
                is_converged =
                    mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
            }

        } // While

        // Check iteration_number
        if (iteration_number >= mMaxIterationNumber) {
            is_converged = true;
            // plots a warning if the maximum number of iterations is exceeded
            if (BaseType::GetModelPart().GetCommunicator().MyPID() == 0) {
                this->MaxIterationsExceeded();
            }
        }

        // calculate reactions if required
        if (mCalculateReactionsFlag) {
            mpBuilderAndSolver->CalculateReactions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }

        BaseType::GetModelPart().GetProcessInfo()[IS_CONVERGED] = is_converged;

        return is_converged;
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        unsigned int iteration_number = BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER];

        // Update the radius
        mRadius = mRadius * sqrt(double(mDesiredIterations) / double(iteration_number));

        DofsArrayType&     rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA      = *mpA;
        TSystemVectorType& mDx     = *mpDx;
        TSystemVectorType& mb      = *mpb;

        if (BaseType::GetModelPart().GetProcessInfo()[IS_CONVERGED]) {
            // Modify the radius to advance faster when convergence is achieved
            if (mRadius > mMaxRadiusFactor * mRadius_0) mRadius = mMaxRadiusFactor * mRadius_0;
            else if (mRadius < mMinRadiusFactor * mRadius_0) mRadius = mMinRadiusFactor * mRadius_0;

            // Update Norm of x
            mNormxEquilibrium = this->CalculateReferenceDofsNorm(rDofSet);
        } else {
            std::cout << "************ NO CONVERGENCE: restoring equilibrium path ************" << std::endl;

            TSystemVectorType& mDxStep = *mpDxStep;

            // update results
            mLambda -= mDLambdaStep;
            noalias(mDx) = -mDxStep;
            this->Update(rDofSet, mA, mDx, mb);

            // move the mesh if needed
            if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
        }

        BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_LAMBDA]        = mLambda;
        BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_RADIUS_FACTOR] = mRadius / mRadius_0;

        mpScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        mpBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        // Cleaning memory after the solution
        mpScheme->Clean();

        if (mReformDofSetAtEachStep) // deallocate the systemvectors
        {
            this->ClearStep();
        }

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        KRATOS_TRY

        SparseSpaceType::Clear(mpf);
        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);

        TSystemVectorType& mf      = *mpf;
        TSystemVectorType& mDxf    = *mpDxf;
        TSystemVectorType& mDxb    = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;

        SparseSpaceType::Resize(mf, 0);
        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        GrandMotherType::Clear();

        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // bool IsConverged() override
    // {
    //     KRATOS_TRY
    //     KRATOS_INFO("Arc-Length:IsConverged()") << std::endl;

    //     bool IsConverged = true;

    //     // Note: Initialize() needs to be called beforehand

    //     this->InitializeSolutionStep();

    //     this->Predict();

    //     // Solve the problem with constant load
    //     IsConverged = this->CheckConvergence();

    //     this->FinalizeSolutionStep();

    //     return IsConverged;

    //     KRATOS_CATCH("")
    // }

    virtual void UpdateLoads()
    {
        KRATOS_TRY

        mLambda = BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_LAMBDA];
        mRadius = (BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_RADIUS_FACTOR]) * mRadius_0;

        // Update External Loads
        this->UpdateExternalLoads();

        KRATOS_CATCH("")
    }

protected:
    /// Member Variables
    TSystemVectorPointerType mpf;      /// Vector of reference external forces
    TSystemVectorPointerType mpDxf;    /// Delta x of A*Dxf=f
    TSystemVectorPointerType mpDxb;    /// Delta x of A*Dxb=b
    TSystemVectorPointerType mpDxPred; /// Delta x of prediction phase
    TSystemVectorPointerType mpDxStep; /// Delta x of the current step

    unsigned int mDesiredIterations; /// This is used to calculate the radius of the next step

    bool mInitializeArcLengthWasPerformed = false;

    double mMaxRadiusFactor,
        mMinRadiusFactor;        /// Used to limit the radius of the arc length strategy
    double mRadius, mRadius_0;   /// Radius of the arc length strategy
    double mLambda, mLambda_old; /// Loading factor
    double mNormxEquilibrium;    /// Norm of the solution vector in equilibrium
    double mDLambdaStep;         /// Delta lambda of the current step

    int Check() override
    {
        KRATOS_TRY

        return MotherType::Check();

        KRATOS_CATCH("")
    }

    void InitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (pv == nullptr) {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), false);
    }

    void SaveInitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (pv == nullptr) {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), true);
    }

    void BuildWithDirichlet(TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
    {
        KRATOS_TRY

        mpBuilderAndSolver->Build(mpScheme, BaseType::GetModelPart(), mA, mb);
        mpBuilderAndSolver->ApplyDirichletConditions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        KRATOS_CATCH("")
    }

    virtual void Update(DofsArrayType& rDofSet, TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
    {
        KRATOS_TRY

        // Update scheme
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Update External Loads
        this->UpdateExternalLoads();

        KRATOS_CATCH("")
    }

    void ClearStep()
    {
        KRATOS_TRY

        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);

        TSystemVectorType& mDxf    = *mpDxf;
        TSystemVectorType& mDxb    = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;

        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        GrandMotherType::Clear();

        KRATOS_CATCH("")
    }

    void UpdateExternalLoads()
    {
        // Update External Loads
        for (unsigned int i = 0; i < mVariableNames.size(); i++) {
            ModelPart&         rSubModelPart = *(mSubModelPartList[i]);
            const std::string& VariableName  = mVariableNames[i];

            if (KratosComponents<Variable<double>>::Has(VariableName)) {
                const Variable<double>& var = KratosComponents<Variable<double>>::Get(VariableName);

#pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(), NodesBegin, NodesEnd);

                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode) {
                        double& rvalue = itNode->FastGetSolutionStepValue(var);
                        rvalue *= (mLambda / mLambda_old);
                    }
                }
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(VariableName)) {
                typedef Variable<double> component_type;
                const component_type&    varx =
                    KratosComponents<component_type>::Get(VariableName + std::string("_X"));
                const component_type& vary =
                    KratosComponents<component_type>::Get(VariableName + std::string("_Y"));
                const component_type& varz =
                    KratosComponents<component_type>::Get(VariableName + std::string("_Z"));

#pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(), NodesBegin, NodesEnd);

                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode) {
                        double& rvaluex = itNode->FastGetSolutionStepValue(varx);
                        rvaluex *= mLambda / mLambda_old;
                        double& rvaluey = itNode->FastGetSolutionStepValue(vary);
                        rvaluey *= mLambda / mLambda_old;
                        double& rvaluez = itNode->FastGetSolutionStepValue(varz);
                        rvaluez *= mLambda / mLambda_old;
                    }
                }
            } else {
                KRATOS_ERROR
                    << "One variable of the applied loads has a non supported type. Variable: " << VariableName
                    << std::endl;
            }
        }

        // Save the applied Lambda factor
        mLambda_old = mLambda;
    }
}; // Class GeoMechanicsRammArcLengthStrategy

} // namespace Kratos