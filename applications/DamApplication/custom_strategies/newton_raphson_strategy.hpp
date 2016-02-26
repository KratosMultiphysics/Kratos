//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined(KRATOS_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_NEWTON_RAPHSON_STRATEGY

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"

#include "dam_application.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class NewtonRaphsonStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(NewtonRaphsonStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef TSparseSpace SparseSpaceType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Constructor
    NewtonRaphsonStrategy(ModelPart& model_part,
                          typename TSchemeType::Pointer pScheme,
                          typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                          double Dofs_Relative_Tolerance,
                          double Residual_Relative_Tolerance,
                          int MaxIterations,
                          bool CalculateReactions,
                          bool ReformDofSetAtEachStep,
                          bool MoveMeshFlag) : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

        //set flags to default values
        mMaxIterationNumber = MaxIterations;
        mDofs_Relative_Tolerance = Dofs_Relative_Tolerance;
        mResidual_Relative_Tolerance = Residual_Relative_Tolerance;
        mCalculateReactionsFlag = CalculateReactions;
        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //saving the scheme
        mpScheme = pScheme;

        //setting up the default builder and solver
        mpBuilderAndSolver = pNewBuilderAndSolver;

        //tells to the builder and solver if the reactions have to be Calculated or not
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        //tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        KRATOS_CATCH( "" )
    }

    //------------------------------------------------------------------------------------

    //Destructor
    virtual ~NewtonRaphsonStrategy() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetEchoLevel(int Level)
    {
        BaseType::mEchoLevel = Level;
        mpBuilderAndSolver->SetEchoLevel(Level);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check()
    {
        KRATOS_TRY

        BaseType::Check();

        mpScheme->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize()
    {
        KRATOS_TRY

        std::cout << "Initializing Newton-Raphson Strategy" << std::endl;

        //Initialize Scheme
        if(mpScheme->SchemeIsInitialized() == false)
            mpScheme->Initialize(BaseType::GetModelPart());

        //Initialize Elements
        if(mpScheme->ElementsAreInitialized() == false)
            mpScheme->InitializeElements(BaseType::GetModelPart());

        //Initialize The Conditions
        if (mpScheme->ConditionsAreInitialized() == false)
            mpScheme->InitializeConditions(BaseType::GetModelPart());

        //set up the system
        if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false)
        {
            //setting up the list of the DOFs to be solved
            mpBuilderAndSolver->SetUpDofSet(mpScheme, BaseType::GetModelPart());

            //shaping correctly the system
            mpBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
        }

        //Initialize system matrix and vectors
        mpBuilderAndSolver->ResizeAndInitializeVectors(mpA, mpDx, mpb, BaseType::GetModelPart().Elements(),
                                                        BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
        InitializeSystemVector(mpTotalx);
        InitializeSystemVector(mpb0);

        //Initialize vectors of initial residual and total x
        TSystemVectorType& mb0 = *mpb0;
        TSparseSpace::SetToZero(mb0);
        mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb0);

        if(TSparseSpace::TwoNorm(mb0) < 1.0e-20)
            KRATOS_THROW_ERROR( std::logic_error, "Norm of the initial residual < 1.0e-20. One must initialize the strategy with some external load applied", "" )
            
        TSystemVectorType& mTotalx = *mpTotalx;
        TSparseSpace::SetToZero(mTotalx);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double Solve()
    {
        KRATOS_TRY

        // ********** Initialize **********

        //initialize solution step
        InitializeSolutionStep();
        
        //Initialize variables
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mDx = *mpDx;

        TSystemVectorType& mb0 = *mpb0;
        TSystemVectorType& mTotalx = *mpTotalx;
        
        // ********** Prediction **********

        //updates the database with a prediction of the solution
        Predict();

        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mDx);
        TSparseSpace::SetToZero(mb);

        mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        mTotalx += mDx;

        //update results
        DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true)
            BaseType::MoveMesh();

        // ********** Correction (iterations loop) **********

        //initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 0;
        double dofs_ratio = 1000.0;
        double residual_ratio = 1000.0;

        while (((dofs_ratio > mDofs_Relative_Tolerance) || (residual_ratio > mResidual_Relative_Tolerance)) && (iteration_number < mMaxIterationNumber))
        {
            //setting the number of iteration
            iteration_number += 1;

            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDx);

            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

            mTotalx += mDx;

            //update results
            mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true)
                BaseType::MoveMesh();

            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            // *** Convergence criterion ***

            dofs_ratio = TSparseSpace::TwoNorm(mDx)/TSparseSpace::TwoNorm(mTotalx);
            std::cout << "    ITERATION " << iteration_number << std::endl;
            std::cout << "        Dofs Ratio = " << dofs_ratio << std::endl;

            TSparseSpace::SetToZero(mb);
            mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
            residual_ratio = TSparseSpace::TwoNorm(mb)/TSparseSpace::TwoNorm(mb0);
            std::cout << "        Residual Ratio = " << residual_ratio << std::endl;
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number == mMaxIterationNumber && (dofs_ratio > mDofs_Relative_Tolerance || residual_ratio > mResidual_Relative_Tolerance))
            std::cout << "************ WARNING: Maximum number of iterations exceeded ************" << std::endl;
        else if (dofs_ratio < mDofs_Relative_Tolerance && residual_ratio < mResidual_Relative_Tolerance)
            std::cout << "Convergence is achieved" << std::endl;

        // ********** Finalize **********

        //calculate reactions if required
        if (mCalculateReactionsFlag == true)
            mpBuilderAndSolver->CalculateReactions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        FinalizeSolutionStep();

        //Deallocate the system vectors and matrix
        if (mReformDofSetAtEachStep == true)
            ClearStep();

        return 0.0;

        KRATOS_CATCH( "" )
    } // Solve()

    //------------------------------------------------------------------------------------

    void InitializeSolutionStep()
    {
        KRATOS_TRY
            
        if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false)
        {
            //setting up the list of the DOFs to be solved
            mpBuilderAndSolver->SetUpDofSet(mpScheme, BaseType::GetModelPart());

            //shaping correctly the system
            mpBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());

            //setting up the system matrix and vectors
            mpBuilderAndSolver->ResizeAndInitializeVectors(mpA, mpDx, mpb, BaseType::GetModelPart().Elements(),
                                                            BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
            InitializeSystemVector(mpb0);
        }

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        mpScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        KRATOS_CATCH( "" )
    }
    
    //------------------------------------------------------------------------------------

    void Predict()
    {
        KRATOS_TRY

        DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mDx = *mpDx;

        mpScheme->Predict(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        //move the mesh if needed
        if (this->MoveMeshFlag() == true)
            BaseType::MoveMesh();

        KRATOS_CATCH( "" )
    }

    //------------------------------------------------------------------------------------

    void FinalizeSolutionStep()
    {
        KRATOS_TRY
        
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        mpScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        mpBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        mpScheme->Clean();

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Clear()
    {
        KRATOS_TRY

        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Clear(mpDx);
        SparseSpaceType::Clear(mpb0);
        SparseSpaceType::Clear(mpTotalx);

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb0 = *mpb0;
        TSystemVectorType& mTotalx = *mpTotalx;

        SparseSpaceType::Resize(mA, 0, 0);
        SparseSpaceType::Resize(mb, 0);
        SparseSpaceType::Resize(mDx, 0);
        SparseSpaceType::Resize(mb0, 0);
        SparseSpaceType::Resize(mTotalx, 0);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->Clear();

        mpScheme->Clear();

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    TSystemMatrixPointerType mpA;

    TSystemVectorPointerType mpb;
    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb0;
    TSystemVectorPointerType mpTotalx;

    typename TSchemeType::Pointer mpScheme;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    double mDofs_Relative_Tolerance, mResidual_Relative_Tolerance;

    unsigned int mMaxIterationNumber;

    bool mReformDofSetAtEachStep, mCalculateReactionsFlag;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (pv == NULL)
        {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), false);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ClearStep()
    {
        KRATOS_TRY

        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Clear(mpDx);

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mDx = *mpDx;

        SparseSpaceType::Resize(mA, 0, 0);
        SparseSpaceType::Resize(mb, 0);
        SparseSpaceType::Resize(mDx, 0);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->Clear();

        mpScheme->Clear();

        KRATOS_CATCH("");
    }

}; // Class NewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_NEWTON_RAPHSON_STRATEGY  defined
