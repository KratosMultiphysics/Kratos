//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_POROMECHANICS_RAMM_ARC_LENGTH_NONLOCAL_STRATEGY)
#define KRATOS_POROMECHANICS_RAMM_ARC_LENGTH_NONLOCAL_STRATEGY

// Project includes
#include "custom_strategies/strategies/poromechanics_ramm_arc_length_strategy.hpp"
#include "custom_utilities/nonlocal_damage_2D_utilities.hpp"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class PoromechanicsRammArcLengthNonlocalStrategy : public PoromechanicsRammArcLengthStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsRammArcLengthNonlocalStrategy);
    
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> Grandx2MotherType;
    typedef PoromechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> GrandMotherType;
    typedef PoromechanicsRammArcLengthStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef TSparseSpace SparseSpaceType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    using Grandx2MotherType::mpScheme;
    using Grandx2MotherType::mpBuilderAndSolver;
    using Grandx2MotherType::mpConvergenceCriteria;
    using Grandx2MotherType::mpA; //Tangent matrix
    using Grandx2MotherType::mpb; //Residual vector of iteration i
    using Grandx2MotherType::mpDx; //Delta x of iteration i
    using Grandx2MotherType::mCalculateReactionsFlag;
    using Grandx2MotherType::mSolutionStepIsInitialized;
    using Grandx2MotherType::mMaxIterationNumber;
    using MotherType::mpf;
    using MotherType::mpDxf;
    using MotherType::mpDxb;
    using MotherType::mpDxPred;
    using MotherType::mpDxStep;
    using MotherType::mRadius;
    using MotherType::mRadius_0;
    using MotherType::mLambda;
    using MotherType::mNormxEquilibrium;
    using MotherType::mDLambdaStep;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoromechanicsRammArcLengthNonlocalStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : PoromechanicsRammArcLengthStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, rParameters, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {            
            if(model_part.GetProcessInfo()[DOMAIN_SIZE]==2)
            {
                mNonlocalDamageUtility = NonlocalDamage2DUtilities(rParameters);
            }
            else
            {
                KRATOS_THROW_ERROR( std::invalid_argument,"NONLOCAL DAMAGE IS NOT AVAILABLE FOR 3D CASES YET", "" )
            }
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~PoromechanicsRammArcLengthNonlocalStrategy() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep()
    {
        KRATOS_TRY

        if (mSolutionStepIsInitialized == false)
		{
            MotherType::InitializeSolutionStep();
            
            mNonlocalDamageUtility.SearchGaussPointsNeighbours(BaseType::GetModelPart());
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	bool SolveSolutionStep()
	{
        // ********** Prediction phase **********
                
        std::cout << "ARC-LENGTH RADIUS: " << mRadius/mRadius_0 << " X initial radius" << std::endl;
        
        // Initialize variables
		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mf = *mpf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;
        
        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        mNonlocalDamageUtility.CalculateNonlocalEquivalentStrain(BaseType::GetModelPart().GetProcessInfo());
        
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDxf);
        
        // Note: This is not so efficient, but I want to solve mA*mDxf=mf without losing mf
        this->BuildWithDirichlet(mA, mDxf, mb);
        noalias(mb) = mf;
        mpBuilderAndSolver->SystemSolve(mA, mDxf, mb);
        
        //update results
        double DLambda = mRadius/TSparseSpace::TwoNorm(mDxf);
        mDLambdaStep = DLambda;
        mLambda += DLambda;
        noalias(mDxPred) = DLambda*mDxf;
        noalias(mDxStep) = mDxPred;
        this->Update(rDofSet, mA, mDxPred, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
        
        // ********** Correction phase (iteration cicle) **********

        //initializing the parameters of the iteration loop
        bool is_converged = false;
        mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDxf, mb);
        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            TSparseSpace::SetToZero(mb);
            mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
        }
        is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDxf, mb);

        unsigned int iteration_number = 0;
        double NormDx;
        
        while (is_converged == false && iteration_number < mMaxIterationNumber)
        {
            //setting the number of iteration
            iteration_number += 1;
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            
            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            mNonlocalDamageUtility.CalculateNonlocalEquivalentStrain(BaseType::GetModelPart().GetProcessInfo());
            
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
            
            DLambda = -TSparseSpace::Dot(mDxPred, mDxb)/TSparseSpace::Dot(mDxPred, mDxf);
            
            noalias(mDx) = mDxb + DLambda*mDxf;
            
            //Check solution before update
            if( mNormxEquilibrium > 1.0e-10 )
            {
                NormDx = TSparseSpace::TwoNorm(mDx);
                
                if( (NormDx/mNormxEquilibrium) > 1.0e3 || (std::abs(DLambda)/std::abs(mLambda-mDLambdaStep)) > 1.0e3 )
                {
                    is_converged = false;
                    break;
                }
            }
            
            //update results
            mDLambdaStep += DLambda;
            mLambda += DLambda;
            noalias(mDxStep) += mDx;
            this->Update(rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
            
            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            mNonlocalDamageUtility.CalculateNonlocalEquivalentStrain(BaseType::GetModelPart().GetProcessInfo());
            
            // *** Check Convergence ***
            
            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(mb);
                mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
            }
            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        }//While
        
        // Check iteration_number 
        if (iteration_number >= mMaxIterationNumber)
        {
            is_converged = true;
            //plots a warning if the maximum number of iterations is exceeded
            if(BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                this->MaxIterationsExceeded();
            }
        }
        
        //calculate reactions if required
        if (mCalculateReactionsFlag == true)
        {
            mpBuilderAndSolver->CalculateReactions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        
        BaseType::GetModelPart().GetProcessInfo()[IS_CONVERGED] = is_converged;
        
		return is_converged;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeSolutionStep()
    {
        KRATOS_TRY
        
        MotherType::FinalizeSolutionStep();
        
        mNonlocalDamageUtility.Clear();

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
        
    NonlocalDamage2DUtilities mNonlocalDamageUtility; //TODO: this should be a general class
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; // Class PoromechanicsRammArcLengthNonlocalStrategy

} // namespace Kratos

#endif // KRATOS_POROMECHANICS_RAMM_ARC_LENGTH_NONLOCAL_STRATEGY  defined