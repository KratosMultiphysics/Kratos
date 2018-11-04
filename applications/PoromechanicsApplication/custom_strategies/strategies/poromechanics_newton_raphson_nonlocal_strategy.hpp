
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_POROMECHANICS_NEWTON_RAPHSON_NONLOCAL_STRATEGY)
#define KRATOS_POROMECHANICS_NEWTON_RAPHSON_NONLOCAL_STRATEGY

// Project includes
#include "custom_strategies/strategies/poromechanics_newton_raphson_strategy.hpp"
#include "custom_utilities/nonlocal_damage_utilities.hpp"
#include "custom_utilities/nonlocal_damage_2D_utilities.hpp"
#include "custom_utilities/nonlocal_damage_3D_utilities.hpp"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class PoromechanicsNewtonRaphsonNonlocalStrategy : public PoromechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsNewtonRaphsonNonlocalStrategy);
    
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> GrandMotherType;
    typedef PoromechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    using GrandMotherType::mpScheme;
    using GrandMotherType::mpBuilderAndSolver;
    using GrandMotherType::mpConvergenceCriteria;
    using GrandMotherType::mpA; //Tangent matrix
    using GrandMotherType::mpb; //Residual vector of iteration i
    using GrandMotherType::mpDx; //Delta x of iteration i
    using GrandMotherType::mCalculateReactionsFlag;
    using GrandMotherType::mSolutionStepIsInitialized;
    using GrandMotherType::mMaxIterationNumber;
    using GrandMotherType::mInitializeWasPerformed;
    using MotherType::mpParameters;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoromechanicsNewtonRaphsonNonlocalStrategy(
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
        ) : PoromechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, rParameters, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag) 
        {
            mNonlocalDamageIsInitialized = false;
            mSearchNeighboursAtEachStep = rParameters["search_neighbours_step"].GetBool();
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~PoromechanicsNewtonRaphsonNonlocalStrategy() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize() override
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false)
		{
            MotherType::Initialize();
            
            if(mNonlocalDamageIsInitialized == false)
            {
                if(BaseType::GetModelPart().GetProcessInfo()[DOMAIN_SIZE]==2)
                {
                    mpNonlocalDamageUtility = new NonlocalDamage2DUtilities();
                }
                else
                {
                    mpNonlocalDamageUtility = new NonlocalDamage3DUtilities();
                }
                mpNonlocalDamageUtility->SearchGaussPointsNeighbours(mpParameters,BaseType::GetModelPart());
                
                mNonlocalDamageIsInitialized = true;
            }
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mSolutionStepIsInitialized == false)
		{
            GrandMotherType::InitializeSolutionStep();
            
            if(mNonlocalDamageIsInitialized == false)
            {
                if(BaseType::GetModelPart().GetProcessInfo()[DOMAIN_SIZE]==2)
                {
                    mpNonlocalDamageUtility = new NonlocalDamage2DUtilities();
                }
                else
                {
                    mpNonlocalDamageUtility = new NonlocalDamage3DUtilities();
                }
                mpNonlocalDamageUtility->SearchGaussPointsNeighbours(mpParameters,BaseType::GetModelPart());
                
                mNonlocalDamageIsInitialized = true;
            }
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	bool SolveSolutionStep() override
	{
        // ********** Prediction phase **********
        
        // Initialize variables
		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain(mpParameters,BaseType::GetModelPart().GetProcessInfo());
        
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDx);

        mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
        
        //update results
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
        
        // ********** Correction phase (iteration cicle) **********

        //initializing the parameters of the iteration loop
        
        bool is_converged = false;
        mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            TSparseSpace::SetToZero(mb);
            mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
        }
        is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        
        unsigned int iteration_number = 0;
        
        while (is_converged == false && iteration_number < mMaxIterationNumber)
        {
            //setting the number of iteration
            iteration_number += 1;
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            
            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain(mpParameters,BaseType::GetModelPart().GetProcessInfo());
            
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDx);
            
            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
                        
            //update results
            mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
            
            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain(mpParameters,BaseType::GetModelPart().GetProcessInfo());
            
            // *** Check Convergence ***
            
            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(mb);
                mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
            }
            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        }//While
        
        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            this->MaxIterationsExceeded();
        }
        
        //calculate reactions if required
        if (mCalculateReactionsFlag == true)
        {
            mpBuilderAndSolver->CalculateReactions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
                
		return is_converged;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY
        
        GrandMotherType::FinalizeSolutionStep();
        
        if(mSearchNeighboursAtEachStep == true)
        {
            delete mpNonlocalDamageUtility;
            mNonlocalDamageIsInitialized = false;
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Clear() override
    {
        KRATOS_TRY

        GrandMotherType::Clear();

        if(mSearchNeighboursAtEachStep == false)
        {
            delete mpNonlocalDamageUtility;
            mNonlocalDamageIsInitialized = false;
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    NonlocalDamageUtilities* mpNonlocalDamageUtility;
    bool mNonlocalDamageIsInitialized;
    bool mSearchNeighboursAtEachStep;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool CheckConvergence() override
    {
        // ********** Prediction phase **********
        
        // Initialize variables
		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;
        
        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain(mpParameters,BaseType::GetModelPart().GetProcessInfo());
                
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDx);
        
        mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
        
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
        
        mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        
        unsigned int iteration_number = 0;
        bool is_converged = false;
        double dofs_ratio = 1000.0;
        double ReferenceDofsNorm;
        double NormDx;
        
        // ********** Correction phase (iteration cicle) **********
        
        while (is_converged == false && iteration_number < mMaxIterationNumber)
        {
            //setting the number of iteration
            iteration_number += 1;
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            
            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain(mpParameters,BaseType::GetModelPart().GetProcessInfo());
            
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDx);
            
            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
            
            mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
            
            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            mpNonlocalDamageUtility->CalculateNonlocalEquivalentStrain(mpParameters,BaseType::GetModelPart().GetProcessInfo());
            
            NormDx = TSparseSpace::TwoNorm(mDx);
            ReferenceDofsNorm = this->CalculateReferenceDofsNorm(rDofSet);
            dofs_ratio = NormDx/ReferenceDofsNorm;
            KRATOS_INFO("Newton Raphson Nonlocal Strategy") << "TEST ITERATION: " << iteration_number << std::endl;
            KRATOS_INFO("Newton Raphson Nonlocal Strategy") << "    Dofs Ratio = " << dofs_ratio << std::endl;
            
            if(dofs_ratio <= 1.0e-3)
                is_converged = true;
        }
        
        return is_converged;
    }

}; // Class PoromechanicsNewtonRaphsonNonlocalStrategy

} // namespace Kratos

#endif // KRATOS_POROMECHANICS_NEWTON_RAPHSON_NONLOCAL_STRATEGY  defined