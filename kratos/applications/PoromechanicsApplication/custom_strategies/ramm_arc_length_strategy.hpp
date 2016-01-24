//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_RAMM_ARC_LENGTH_STRATEGY)
#define KRATOS_RAMM_ARC_LENGTH_STRATEGY

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"

#include "poromechanics_application.h"

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
        InitializeSystemVector(mpf);
        InitializeSystemVector(mpDxf);
        InitializeSystemVector(mpDxb);
        InitializeSystemVector(mpDxPred);
        InitializeSystemVector(mpDxStep);

        //Initialize vector of reference external force, vector of total x, and initial radius
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mTotalx = *mpTotalx;
        TSystemVectorType& mf = *mpf;
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mTotalx);
        TSparseSpace::SetToZero(mf);

        mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mTotalx, mf);

        mRadius_0 = TSparseSpace::TwoNorm(mTotalx);
        mRadius = mRadius_0;

        TSparseSpace::SetToZero(mTotalx);

        //Initialize the loading factor Lambda
        mLambda = 0.0;
        mLambda_old = 1;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double Solve()
    {
        KRATOS_TRY

        // ********** Initialize **********

        //initialize solution step
        InitializeSolutionStep();

        std::cout << "RADIUS: " << mRadius << " (" << mRadius/mRadius_0 << " x initial radius)" << std::endl;
        
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
        Predict();

        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDxf);

        BuildWithDirichlet(mA, mDx, mb);
        
        mpBuilderAndSolver->SystemSolve(mA, mDxf, mf);

        double DLambda = mRadius/TSparseSpace::TwoNorm(mDxf);
        double DLambdaStep = DLambda;
        mLambda += DLambda;
        
        mDxPred = DLambda*mDxf;
        mDxStep = mDxPred;
        mTotalx += mDxPred;

        //update results
        Update(mA, mDxPred, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true)
            BaseType::MoveMesh();

        // ********** Correction (iterations loop) **********

        //initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 0;
        double dofs_ratio = 1000, residual_ratio = 1000;
        bool possible_convergence = true;
        double NormDx, Normx;
        double Normf = TSparseSpace::TwoNorm(mf);

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
            if( isnan(NormDx) || isinf(NormDx) || isnan(DLambda) || isinf(DLambda) )
            {
                possible_convergence = false;
                break;
            }
            else if(Normx > 1e-10)
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
            Update(mA, mDx, mb);

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

        if (iteration_number == mMaxIterationNumber && (dofs_ratio > mDofs_Relative_Tolerance || residual_ratio > mResidual_Relative_Tolerance))
        {
            std::cout << "************ WARNING: Maximum number of iterations exceeded ************" << std::endl;
            possible_convergence = false;
        }
        else if (dofs_ratio < mDofs_Relative_Tolerance && residual_ratio < mResidual_Relative_Tolerance)
            std::cout << "Convergence is achieved" << std::endl;

        // ********** Finalize **********

        //calculate reactions if required
        if (mCalculateReactionsFlag == true && possible_convergence == true)
            mpBuilderAndSolver->CalculateReactions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        FinalizeSolutionStep(iteration_number, possible_convergence, DLambdaStep);

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
            InitializeSystemVector(mpDxf);
            InitializeSystemVector(mpDxb);
            InitializeSystemVector(mpDxPred);
            InitializeSystemVector(mpDxStep);
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
            Update(mA, mDx, mb);

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

        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin(); i != BaseType::GetModelPart().NodesEnd(); ++i)
        {
            array_1d<double,3>& PointLoad = (i)->FastGetSolutionStepValue(POINT_LOAD);
            PointLoad = (mLambda/mLambda_old)*PointLoad;

            array_1d<double,3>& LineLoad = (i)->FastGetSolutionStepValue(LINE_LOAD);
            LineLoad = (mLambda/mLambda_old)*LineLoad;

            array_1d<double,3>& SurfaceLoad = (i)->FastGetSolutionStepValue(SURFACE_LOAD);
            SurfaceLoad = (mLambda/mLambda_old)*SurfaceLoad;

            double& NormalContactStress = (i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
            NormalContactStress = (mLambda/mLambda_old)*NormalContactStress;

            double& TangentialContactStress = (i)->FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS);
            TangentialContactStress = (mLambda/mLambda_old)*TangentialContactStress;

            double& FluidFlux = (i)->FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
            FluidFlux = (mLambda/mLambda_old)*FluidFlux;
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
