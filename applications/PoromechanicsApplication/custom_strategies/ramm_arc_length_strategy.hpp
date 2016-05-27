//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_RAMM_ARC_LENGTH_STRATEGY)
#define KRATOS_RAMM_ARC_LENGTH_STRATEGY

// Project includes
#include "custom_strategies/newton_raphson_strategy.hpp"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class RammArcLengthStrategy : public NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
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
    typedef NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    using MotherType::mpA; //Tangent matrix
    using MotherType::mpb; //Residual vector of iteration i
    using MotherType::mpDx; //Delta x of iteration i
    using MotherType::mpTotalx; //Total accumulated x
    using MotherType::mpScheme;
    using MotherType::mpBuilderAndSolver;
    using MotherType::mDofs_Relative_Tolerance;
    using MotherType::mResidual_Relative_Tolerance;
    using MotherType::mMaxIterationNumber;
    using MotherType::mReformDofSetAtEachStep;
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mNormb0;

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
                          bool MoveMeshFlag) : NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                                                                    pNewBuilderAndSolver,Dofs_Relative_Tolerance,Residual_Relative_Tolerance,
                                                                    MaxIterations,CalculateReactions,ReformDofSetAtEachStep,MoveMeshFlag)
    {
        KRATOS_TRY

        //set flags to default values
        mDesiredIterations = DesiredIterations;
        mMaxRadiusFactor = MaxRadius;
        mMinRadiusFactor = MinRadius;
        
        KRATOS_CATCH( "" )
    }

    //------------------------------------------------------------------------------------

    //Destructor
    virtual ~RammArcLengthStrategy() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check()
    {
        KRATOS_TRY
        
        int err = MotherType::Check();
        if(err != 0) return err;
                
        //check for variables keys (verify that the variables are correctly initialized)
        if(FORCE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"FORCE has Key zero! (check if the application is correctly registered", "" )
        if(FACE_LOAD.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"FACE_LOAD has Key zero! (check if the application is correctly registered", "" )
        if(NORMAL_CONTACT_STRESS.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"NORMAL_CONTACT_STRESS has Key zero! (check if the application is correctly registered", "" )
        if(TANGENTIAL_CONTACT_STRESS.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"TANGENTIAL_CONTACT_STRESS has Key zero! (check if the application is correctly registered", "" )
        if(NORMAL_FLUID_FLUX.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"NORMAL_FLUID_FLUX has Key zero! (check if the application is correctly registered", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=BaseType::GetModelPart().NodesBegin(); it!=BaseType::GetModelPart().NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(FORCE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "FORCE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(FACE_LOAD) == false)
                KRATOS_THROW_ERROR( std::logic_error, "FACE_LOAD variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(NORMAL_CONTACT_STRESS) == false)
                KRATOS_THROW_ERROR( std::logic_error, "NORMAL_CONTACT_STRESS variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(TANGENTIAL_CONTACT_STRESS) == false)
                KRATOS_THROW_ERROR( std::logic_error, "TANGENTIAL_CONTACT_STRESS variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(NORMAL_FLUID_FLUX) == false)
                KRATOS_THROW_ERROR( std::logic_error, "NORMAL_FLUID_FLUX variable is not allocated for node ", it->Id() )
        }

        return err;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize()
    {
        KRATOS_TRY

        std::cout << "Initializing Ramm's Arc Length Strategy" << std::endl;

        //Initialize ProcessInfo variables
        BaseType::GetModelPart().GetProcessInfo()[NO_CONVERGENCE] = 0;
        
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

        mNormb0 = TSparseSpace::TwoNorm(mf);
        if(mNormb0 < 1.0e-20)
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
                if( (NormDx/Normx) > 5e2 || (fabs(DLambda)/fabs(mLambda-DLambdaStep)) > 5e2 || (TSparseSpace::TwoNorm(mb)/mNormb0) > 5e2 )
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
            residual_ratio = TSparseSpace::TwoNorm(mb)/mNormb0;
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

        MotherType::FinalizeSolutionStep();
        
        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Clear()
    {
        KRATOS_TRY
        
        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);
        SparseSpaceType::Clear(mpf);

        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;
        TSystemVectorType& mf = *mpf;

        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);
        SparseSpaceType::Resize(mf, 0);
        
        MotherType::Clear();
        
        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    unsigned int mDesiredIterations; //This is used to calculate the radius of the next step
    double mMaxRadiusFactor, mMinRadiusFactor; //Used to limit the radius of the arc length strategy
    double mRadius, mRadius_0; //Radius of the arc length strategy
    double mLambda, mLambda_old; // Loading factor
    TSystemVectorPointerType mpf; //Vector of reference external forces
    TSystemVectorPointerType mpDxf; //Delta x of A*Dxf=f
    TSystemVectorPointerType mpDxb; //Delta x of A*Dxb=b
    TSystemVectorPointerType mpDxPred; //Delta x of prediction phase
    TSystemVectorPointerType mpDxStep; //Delta x of the current step

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

        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);

        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;

        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        MotherType::ClearStep();

        KRATOS_CATCH("");
    }

}; // Class RammArcLengthStrategy

} // namespace Kratos

#endif // KRATOS_RAMM_ARC_LENGTH_STRATEGY  defined
