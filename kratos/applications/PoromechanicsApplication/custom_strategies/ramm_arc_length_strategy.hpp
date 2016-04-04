//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_RAMM_ARC_LENGTH_STRATEGY)
#define KRATOS_RAMM_ARC_LENGTH_STRATEGY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class RammArcLengthStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(RammArcLengthStrategy);

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
    RammArcLengthStrategy(ModelPart& model_part,
                          typename TSchemeType::Pointer pScheme,
                          typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                          double Dofs_Relative_Tolerance,
                          double Residual_Relative_Tolerance,
                          int MaxIterations,
                          int DesiredIterations,
                          double MaxRadius,
                          double MinRadius,
                          bool CalculateReactions,
                          bool ReformDofSetAtEachStep,
                          bool MoveMeshFlag) : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

        //saving the scheme
        mpScheme = pScheme;

        //setting up the default builder and solver
        mpBuilderAndSolver = pNewBuilderAndSolver;

        //set flags to default values
        mDofs_Relative_Tolerance = Dofs_Relative_Tolerance;
        mResidual_Relative_Tolerance = Residual_Relative_Tolerance;

        mMaxIterationNumber = MaxIterations;
        mDesiredIterations = DesiredIterations;

        mMaxRadiusFactor = MaxRadius;
        mMinRadiusFactor = MinRadius;

        mCalculateReactionsFlag = CalculateReactions;
        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //tells to the builder and solver if the reactions have to be Calculated or not
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        //tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        //Initialize ProcessInfo variables
        model_part.GetProcessInfo()[NO_CONVERGENCE] = 0;
        
        KRATOS_CATCH( "" )
    }

    //------------------------------------------------------------------------------------

    //Destructor
    virtual ~RammArcLengthStrategy() {}

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

        std::cout << "Initializing Ramm's Arc Length Strategy" << std::endl;

        //Initialize Scheme
        if(mpScheme->SchemeIsInitialized() == false)
            mpScheme->Initialize(BaseType::GetModelPart());

        //Initialize Elements
        if(mpScheme->ElementsAreInitialized() == false)
            mpScheme->InitializeElements(BaseType::GetModelPart());

        //Initialize Conditions
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
        this->InitializeSystemVector(mpTotalx);
        this->InitializeSystemVector(mpf);
        this->InitializeSystemVector(mpDxf);
        this->InitializeSystemVector(mpDxb);
        this->InitializeSystemVector(mpDxPred);
        this->InitializeSystemVector(mpDxStep);

        //Initialize vector of reference external force, vector of total x, and initial radius
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mTotalx = *mpTotalx;
        TSystemVectorType& mf = *mpf;
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mTotalx);
        TSparseSpace::SetToZero(mf);

        mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mTotalx, mf);

        double Normf = TSparseSpace::TwoNorm(mf);
        if(Normf < 1.0e-20)
            KRATOS_THROW_ERROR( std::logic_error, "Norm of the initial residual < 1.0e-20. One must initialize the strategy with some external load applied", "" )
            
        mRadius_0 = TSparseSpace::TwoNorm(mTotalx);
        mRadius = mRadius_0;

        TSparseSpace::SetToZero(mTotalx);

        //Initialize the loading factor Lambda
        mLambda = 0.0;
        mLambda_old = 1.0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double Solve()
    {
        KRATOS_TRY

        // ********** Initialize **********

        //initialize solution step
        this->InitializeSolutionStep();

        std::cout << "ARC-LENGTH RADIUS: " << mRadius << " (" << mRadius/mRadius_0 << " x initial radius)" << std::endl;
        
        //Initialize variables
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mDxb = *mpDxb;

        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mf = *mpf;
        TSystemVectorType& mDxf = *mpDxf;

        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;
        TSystemVectorType& mTotalx = *mpTotalx;
        
        // ********** Prediction **********

        //updates the database with a prediction of the solution
        this->Predict();

        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDxf);

        this->BuildWithDirichlet(mA, mDx, mb);
        
        mpBuilderAndSolver->SystemSolve(mA, mDxf, mf);

        double DLambda = mRadius/TSparseSpace::TwoNorm(mDxf);
        double DLambdaStep = DLambda;
        mLambda += DLambda;
        
        mDxPred = DLambda*mDxf;
        mDxStep = mDxPred;
        mTotalx += mDxPred;

        //update results
        this->Update(mA, mDxPred, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true)
            BaseType::MoveMesh();

        // ********** Correction (iterations loop) **********

        //initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 0;
        double dofs_ratio = 1000.0;
        double residual_ratio = 1000.0;
        bool possible_convergence = true;
        double NormDx, Normx;
        const double Normf = TSparseSpace::TwoNorm(mf);

        while ((dofs_ratio > mDofs_Relative_Tolerance || residual_ratio > mResidual_Relative_Tolerance) && iteration_number < mMaxIterationNumber && possible_convergence == true)
        {
            //setting the number of iteration
            iteration_number += 1;

            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxb);

            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDxb, mb);

            TSparseSpace::SetToZero(mDxf);

            mpBuilderAndSolver->SystemSolve(mA, mDxf, mf);

            DLambda = -TSparseSpace::Dot(mDxPred, mDxb)/TSparseSpace::Dot(mDxPred, mDxf);
            
            mDx = mDxb + DLambda*mDxf;

            //Check solution
            NormDx = TSparseSpace::TwoNorm(mDx);
            Normx = TSparseSpace::TwoNorm(mTotalx-mDxStep);
            if(Normx > 1e-10)
            {
                if( (NormDx/Normx) > 5e2 || (fabs(DLambda)/fabs(mLambda-DLambdaStep)) > 5e2 || (TSparseSpace::TwoNorm(mb)/Normf) > 5e2 )
                {
                    possible_convergence = false;
                    break;
                }
            }
            
            DLambdaStep += DLambda;
            mLambda += DLambda;
            
            mDxStep += mDx;
            mTotalx += mDx;

            //update results
            this->Update(mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true)
                BaseType::MoveMesh();

            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            // *** Convergence criterion ***

            dofs_ratio = NormDx/TSparseSpace::TwoNorm(mTotalx);
            std::cout << "    ITERATION " << iteration_number << std::endl;
            std::cout << "        Dofs Ratio = " << dofs_ratio << std::endl;

            TSparseSpace::SetToZero(mb);
            mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
            residual_ratio = TSparseSpace::TwoNorm(mb)/Normf;
            std::cout << "        Residual Ratio = " << residual_ratio << std::endl;
        } //while

        if (dofs_ratio < mDofs_Relative_Tolerance && residual_ratio < mResidual_Relative_Tolerance)
            std::cout << "Convergence is achieved" << std::endl;
        else if (iteration_number == mMaxIterationNumber && (dofs_ratio > mDofs_Relative_Tolerance || residual_ratio > mResidual_Relative_Tolerance))
        {
            std::cout << "************ WARNING: Maximum number of iterations exceeded ************" << std::endl;
            possible_convergence = false;
        }

        // ********** Finalize **********

        //calculate reactions if required
        if (mCalculateReactionsFlag == true && possible_convergence == true)
            mpBuilderAndSolver->CalculateReactions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        this->FinalizeSolutionStep(iteration_number, possible_convergence, DLambdaStep);

        //Deallocate the system vectors and matrix
        if (mReformDofSetAtEachStep == true)
            this->ClearStep();

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
            this->InitializeSystemVector(mpDxf);
            this->InitializeSystemVector(mpDxb);
            this->InitializeSystemVector(mpDxPred);
            this->InitializeSystemVector(mpDxStep);
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

    void FinalizeSolutionStep(unsigned int iteration_number, bool possible_convergence, double DLambdaStep)
    {
        KRATOS_TRY

        mRadius = mRadius*sqrt(double(mDesiredIterations)/double(iteration_number));

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        if (possible_convergence == true)
        {
            if (mRadius > mMaxRadiusFactor*mRadius_0)
                mRadius = mMaxRadiusFactor*mRadius_0;
            else if(mRadius < mMinRadiusFactor*mRadius_0)
                mRadius = mMinRadiusFactor*mRadius_0;

            BaseType::GetModelPart().GetProcessInfo()[NO_CONVERGENCE] = 0;
        }
        else
        {
            std::cout << "************ IMPOSSIBLE CONVERGENCE: redressing equilibrium path ************" << std::endl;

            TSystemVectorType& mDxStep = *mpDxStep;
            TSystemVectorType& mTotalx = *mpTotalx;

            mLambda -= DLambdaStep;
            mTotalx -= mDxStep;

            mDx = -mDxStep;

            //update results
            this->Update(mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true)
                BaseType::MoveMesh();

            BaseType::GetModelPart().GetProcessInfo()[NO_CONVERGENCE] = 1;
        }

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
        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);
        SparseSpaceType::Clear(mpTotalx);
        SparseSpaceType::Clear(mpf);

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;
        TSystemVectorType& mTotalx = *mpTotalx;
        TSystemVectorType& mf = *mpf;

        SparseSpaceType::Resize(mA, 0, 0);
        SparseSpaceType::Resize(mb, 0);
        SparseSpaceType::Resize(mDx, 0);
        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);
        SparseSpaceType::Resize(mTotalx, 0);
        SparseSpaceType::Resize(mf, 0);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->Clear();

        mpScheme->Clear();

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    TSystemMatrixPointerType mpA; //Tangent matrix

    TSystemVectorPointerType mpb; //Residual vector of iteration i
    TSystemVectorPointerType mpDx; //Delta x of iteration i
    TSystemVectorPointerType mpTotalx; //Total accumulated x

    typename TSchemeType::Pointer mpScheme;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    double mDofs_Relative_Tolerance, mResidual_Relative_Tolerance;

    unsigned int mMaxIterationNumber;

    bool mReformDofSetAtEachStep, mCalculateReactionsFlag;

    // Arc Length Variables
    unsigned int mDesiredIterations; //This is used to calculate the radius of the next step
    double mMaxRadiusFactor, mMinRadiusFactor; //Used to limit the radius of the arc length strategy
    double mRadius, mRadius_0; //Radius of the arc length strategy
    double mLambda, mLambda_old; // Loading factor
    TSystemVectorPointerType mpf; //Vector of reference external forces
    TSystemVectorPointerType mpDxf; //Delta x of A*Dxf=f
    TSystemVectorPointerType mpDxb; //Delta x of A*Dxb=b
    TSystemVectorPointerType mpDxPred; //Delta x of prediction phase

    // Recovery Variables
    TSystemVectorPointerType mpDxStep; //Delta x of the current step

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

    void BuildWithDirichlet(TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
    {
        KRATOS_TRY

        mpBuilderAndSolver->Build(mpScheme, BaseType::GetModelPart(), mA, mb);
        mpBuilderAndSolver->ApplyDirichletConditions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Update(TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
    {
        KRATOS_TRY

        DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        
        //Update External Conditions
        
        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(5)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(5)->NodesEnd(); i++)
        {
            double& PointLoadX = (i)->FastGetSolutionStepValue(FORCE_X);
            PointLoadX *= (mLambda/mLambda_old);
        }

        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(6)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(6)->NodesEnd(); i++)
        {
            double& PointLoadY = (i)->FastGetSolutionStepValue(FORCE_Y);
            PointLoadY *= (mLambda/mLambda_old);
        }

        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(7)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(7)->NodesEnd(); i++)
        {
            double& PointLoadZ = (i)->FastGetSolutionStepValue(FORCE_Z);
            PointLoadZ *= (mLambda/mLambda_old);
        }
        
        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(8)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(8)->NodesEnd(); i++)
        {
            double& LineLoadX = (i)->FastGetSolutionStepValue(FACE_LOAD_X);
            LineLoadX *= (mLambda/mLambda_old);
        }

        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(9)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(9)->NodesEnd(); i++)
        {
            double& LineLoadY = (i)->FastGetSolutionStepValue(FACE_LOAD_Y);
            LineLoadY *= (mLambda/mLambda_old);
        }

        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(10)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(10)->NodesEnd(); i++)
        {
            double& LineLoadZ = (i)->FastGetSolutionStepValue(FACE_LOAD_Z);
            LineLoadZ *= (mLambda/mLambda_old);
        }

        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(11)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(11)->NodesEnd(); i++)
        {
            double& NormalContactStress = (i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
            NormalContactStress *= (mLambda/mLambda_old);
        }

        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(12)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(12)->NodesEnd(); i++)
        {
            double& TangentialContactStress = (i)->FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS);
            TangentialContactStress *= (mLambda/mLambda_old);
        }

        for(ModelPart::NodeIterator i = BaseType::GetModelPart().pGetMesh(13)->NodesBegin(); i < BaseType::GetModelPart().pGetMesh(13)->NodesEnd(); i++)
        {
            double& FluidFlux = (i)->FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
            FluidFlux *= (mLambda/mLambda_old);
        }
        
        mLambda_old = mLambda;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ClearStep()
    {
        KRATOS_TRY

        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Clear(mpDx);
        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;

        SparseSpaceType::Resize(mA, 0, 0);
        SparseSpaceType::Resize(mb, 0);
        SparseSpaceType::Resize(mDx, 0);
        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->Clear();

        mpScheme->Clear();

        KRATOS_CATCH("");
    }

}; // Class RammArcLengthStrategy

} // namespace Kratos

#endif // KRATOS_RAMM_ARC_LENGTH_STRATEGY  defined
