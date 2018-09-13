//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $

#if !defined(KRATOS_ARC_LENGTH_STRATEGY )
#define  KRATOS_ARC_LENGTH_STRATEGY

/* External includes */
#include "boost/smart_ptr.hpp"
#include <iomanip>

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "custom_strategies/builder_and_solvers/static_general_builder_and_solver.h"


#define ALEN_HSP 1
#define ALEN_CYL 2

#define ALEN_TYPE ALEN_CYL

#define ALEN_CYL_ALPHA2 0.1

namespace Kratos
{

/// Short class definition.
/**   Detail class definition.
 */
template<class TSparseSpace,
         class TDenseSpace, 
         class TLinearSolver
         >
class ArcLengthStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(ArcLengthStrategy);

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
	
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef StaticGeneralBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > TDefaultBuilderAndSolverType;

public:

    ArcLengthStrategy(ModelPart& model_part,
                      typename TSchemeType::Pointer pScheme,
                      typename TLinearSolver::Pointer pNewLinearSolver,
                      typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
					  double arcLength, 
					  double psi_u, 
					  double psi_f, 
					  double u_ref,
                      int MaxIterations = 30,
                      bool CalculateReactions = false,
                      bool ReformDofSetAtEachStep = false,
                      bool MoveMeshFlag = false
					  )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
		, mArcLength2(arcLength*arcLength)
		, mArcLength_0(arcLength)
		, mDeltaTime_0(-1.0)
		, mPsiU2(psi_u*psi_u)
		, mPsiF2(psi_f*psi_f)
		, mURef2(u_ref*u_ref)
		, mDeltaUHat(NULL)
		, mDeltaUBar(NULL)
		, mDeltaU(NULL)
		, mDeltaUStep(NULL)
		, mPHat(NULL)
		, mDeltaLambdaStep(0.0)
		, mCurrentLambda(0.0)
		, mSignLastDeltaLambdaStep(0.0)
		, mArcLengthInitialized(false)
	{        
        mIsConverged = false;

        mKeepSystemConstantDuringIterations = false;

        //set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        //saving the scheme
        mpScheme = pScheme;

        //saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        //setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new TDefaultBuilderAndSolverType(mpLinearSolver));

        //set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        //tells to the builder and solver if the reactions have to be Calculated or not
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        //tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        //set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        //by default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);
    }

    virtual ~ArcLengthStrategy()
    {
    }

public:
	
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    }
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    }

    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    }
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    }

    void SetInitializePerformedFlag(bool InitializePerformedFlag)
    {
		mInitializeWasPerformed = InitializePerformedFlag;
    }
    bool GetInitializePerformedFlag()
    {
		return mInitializeWasPerformed;
    }

    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        mReformDofSetAtEachStep = flag;
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

	void SetKeepSystemConstantDuringIterations(bool value)
    {
        mKeepSystemConstantDuringIterations = value;
    }
    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }
	
    void SetMaxIterationNumber(unsigned int MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }
    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    void SetEchoLevel(int Level)
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

	TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;
        return mA;
    }
	
public:

    void Predict()
    {
        /*if (mInitializeWasPerformed == false) Initialize();
        if (mSolutionStepIsInitialized == false) InitializeSolutionStep();
        DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        mpScheme->Predict(BaseType::GetModelPart(), rDofSet, *mpA, *mpDx, *mpb);
        if (this->MoveMeshFlag()) BaseType::MoveMesh();*/
    }

	void PrepareFirstSolve()
	{
		/* 
		===============================================================================
		Solution step initialization
		===============================================================================
		*/

		mIsConverged = false;

		ModelPart& model = BaseType::GetModelPart();

		mCurrentLambda = model.GetProcessInfo()[LAMBDA];

		if (!mInitializeWasPerformed) this->Initialize();

		int finalize_level   = model.GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL];
		bool verbose_allowed = model.GetCommunicator().MyPID() == 0;
		bool verbose         = this->GetEchoLevel() > 0 && verbose_allowed;

		model.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 0;

		if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true)
		{
			mpBuilderAndSolver->SetUpDofSet(mpScheme, model);
			mpBuilderAndSolver->SetUpSystem(model);
		}

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();

		//if (verbose) PromptCurrentTime();

		this->Predict(); 

		// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		//if (!mSolutionStepIsInitialized) this->InitializeSolutionStep();
		// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if (!mSolutionStepIsInitialized) 
		{
			this->InitializeSolutionStep();
			this->ArcLength_Initialize();
			this->ArcLength_InitializeSolutionStep();
		}
		else
		{
			this->ArcLength_Initialize();
		}
		// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	}
	
	double Solve()
	{
		/* 
		===============================================================================
		Solution step initialization
		===============================================================================
		*/

		mIsConverged = false;

		ModelPart& model = BaseType::GetModelPart();

		if (!mInitializeWasPerformed) this->Initialize();

		int finalize_level   = model.GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL];
		bool verbose_allowed = model.GetCommunicator().MyPID() == 0;
		bool verbose         = this->GetEchoLevel() > 0 && verbose_allowed;

		model.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 0;

		if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true)
		{
			mpBuilderAndSolver->SetUpDofSet(mpScheme, model);
			mpBuilderAndSolver->SetUpSystem(model);
		}

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();

		if (verbose) PromptCurrentTime();

		this->Predict(); 

		// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		//if (!mSolutionStepIsInitialized) this->InitializeSolutionStep();
		// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if (!mSolutionStepIsInitialized) 
		{
			this->InitializeSolutionStep();
			this->ArcLength_Initialize();
			this->ArcLength_InitializeSolutionStep();
		}
		else
		{
			this->ArcLength_Initialize();
		}
		// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		TSystemMatrixType& mA	= *mpA;
		TSystemVectorType& mDx	= *mpDx;
		TSystemVectorType& mb	= *mpb;

		if(finalize_level == RVE_STRATEGY_FINALIZE_STEP_LEVEL___FINALIZE_ONLY)
		{
			// just call the finalize method, but first updated the RHS.
			TSparseSpace::SetToZero(mb);
			mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
			model.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 1;
			this->FinalizeSolutionStep(mA, mDx, mb);
			return 0.0;
		}
		else if(finalize_level == RVE_STRATEGY_FINALIZE_STEP_LEVEL___ABORT_ONLY)
		{
			// abort the previous solution step.
			// the trial values will not be committed
			this->AbortSolutionStep();
			return 0.0;
		}

		if(TSparseSpace::Size(mDx) == 0)
		{
			/* 
			===============================================================================
			If all DOFs are set, just update the RHS and move on
			===============================================================================
			*/

			mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
			if(this->CheckIntegrationErrors() != 0.0)
			{
				if(verbose) this->PromptIntegrationErrors();
				mIsConverged = false;
				model.GetProcessInfo()[NL_ITERATION_NUMBER] = mMaxIterationNumber;
			}
			else
			{
				if (verbose) std::cout << "No Free DOFs - Finalizing solution step" << std::endl;
				mIsConverged = true;
			}
		}
		else
		{
			/* 
			===============================================================================
			Begin the N-R equilibrium iteration loop
			===============================================================================
			*/
			unsigned int iteration_number = 0;
			model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

			// for initial unbalance
			TSystemVectorType& PHat = *mPHat;
			//mpScheme->InitializeNonLinIteration(model, mA, mDx, mb);
			//mIsConverged = mpConvergenceCriteria->PreCriteria(model, rDofSet, mA, mDx, mb);
			//TSparseSpace::SetToZero(mDx);
			//TSparseSpace::SetToZero(mb);
			//mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
			//noalias(mb) += mCurrentLambda * PHat; // Note: we cannot increase loads here in c++ (we do it in python)
			//mIsConverged = mpConvergenceCriteria->PostCriteria(model, rDofSet, mA, mDx, mb); // to initialize the criterion

			while(true)
			{
				/* 
				===============================================================================
				Initialize the current nonlinear iteration
				===============================================================================
				*/
				iteration_number++;
				model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
				mpScheme->InitializeNonLinIteration(model, mA, mDx, mb);
				mIsConverged = mpConvergenceCriteria->PreCriteria(model, rDofSet, mA, mDx, mb);
				/* 
				===============================================================================
				Build system matrix
				===============================================================================
				*/
				TSparseSpace::SetToZero(mA);
				//mpBuilderAndSolver->BuildLHS(mpScheme, model, mA); // LHS
				TSparseSpace::SetToZero(mb);
				mpBuilderAndSolver->Build(mpScheme, model, mA, mb);
				noalias(mb) += mCurrentLambda * PHat; // Note: we cannot increase loads here in c++ (we do it in python)

				/* 
				===============================================================================
				Abort iteration process in case of integration errors,
				otherwise solve the system of equations
				===============================================================================
				*/
				if(this->CheckIntegrationErrors() != 0.0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				else
				{
					TSparseSpace::SetToZero(mDx);
					mpBuilderAndSolver->SystemSolve(mA, mDx, mb); // solve for DU
				}
				/* 
				===============================================================================
				Update with arc length
				===============================================================================
				*/
				rDofSet = mpBuilderAndSolver->GetDofSet();
				bool update_res = this->ArcLength_Update();

				if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
				if(!update_res) {
					if(verbose) this->PromptArcLengthUpdateErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				/* 
				===============================================================================
				Build updated residual vector
				===============================================================================
				*/
				TSparseSpace::SetToZero(mb);
				mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
				noalias(mb) += mCurrentLambda * PHat; // Note: we cannot increase loads here in c++ (we do it in python)
				/* 
				===============================================================================
				Abort iteration process in case of integration errors
				===============================================================================
				*/
				if(this->CheckIntegrationErrors() != 0.0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				/* 
				===============================================================================
				Finalize nonlinear iteration
				===============================================================================
				*/
				mpScheme->FinalizeNonLinIteration(model, mA, mDx, mb);
				/* 
				===============================================================================
				Check convergence criteria.
				Also check for divergence and if so, abort solution step.
				===============================================================================
				*/
				model.GetProcessInfo()[ITERATION_CONVERGENCE_FLAG] = 0;
				mIsConverged = mpConvergenceCriteria->PostCriteria(model, rDofSet, mA, mDx, mb);
				if(model.GetProcessInfo()[ITERATION_CONVERGENCE_FLAG] != 0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				if(mIsConverged || (iteration_number > mMaxIterationNumber) ) break;
			}

			// refesh RHS (we DON'T do this in BuilderAndSolver to check for integration errors!!!!)
			TSparseSpace::SetToZero(mb);
			mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
			noalias(mb) += mCurrentLambda * PHat; // Note: we cannot increase loads here in c++ (we do it in python)

			if(this->CheckIntegrationErrors() != 0.0)
			{
				if(verbose) this->PromptIntegrationErrors();
				mIsConverged = false;
				iteration_number = mMaxIterationNumber;
				model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
			}
		}

		/* 
		===============================================================================
		Perform the finalization of this solution step
		based on the STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL
		===============================================================================
		*/

		model.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 1;

		if(finalize_level == RVE_STRATEGY_FINALIZE_STEP_LEVEL___ALWAYS_FINALIZE)
		{
			if(verbose && (mIsConverged == false)) PromptMaxIterationsExceeded();
			if(verbose) PromptSeparator();
			this->FinalizeSolutionStep(mA, mDx, mb);
			if(verbose) PromptSeparator();
		}
		else if(finalize_level == RVE_STRATEGY_FINALIZE_STEP_LEVEL___FINALIZE_OR_ABORT)
		{
			if(mIsConverged)
			{
				if(verbose) PromptSeparator();
				this->FinalizeSolutionStep(mA, mDx, mb);
				if(verbose) PromptSeparator();
			}
			else
			{
				if(verbose) PromptMaxIterationsExceeded();
				this->AbortSolutionStep();
			}
		}
		else if(finalize_level == RVE_STRATEGY_FINALIZE_STEP_LEVEL___ALWAYS_ABORT)
		{
			if(verbose && (mIsConverged == false)) PromptMaxIterationsExceeded();
			this->AbortSolutionStep();
		}
		else if(finalize_level == RVE_STRATEGY_FINALIZE_STEP_LEVEL___DO_NOTHING)
		{
			if(verbose && (mIsConverged == false)) PromptMaxIterationsExceeded();
		}

		return 0.00;
	} 

    bool IsConverged()
    {
		return mIsConverged;
    }

    void CalculateOutputData()
    {
        mpScheme->CalculateOutputData(BaseType::GetModelPart(), mpBuilderAndSolver->GetDofSet(), *mpA, *mpDx, *mpb);
    }

    void Clear()
    {
        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Resize(*mpA, 0, 0);

        SparseSpaceType::Clear(mpDx);
        SparseSpaceType::Resize(*mpDx, 0);

        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Resize(*mpb, 0);

        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->Clear();

        mpScheme->Clear();
    }

protected:

    void Initialize()
    {
		ModelPart & model = BaseType::GetModelPart();

        if(!mpScheme->SchemeIsInitialized())
			mpScheme->Initialize(model);

        if(!mpScheme->ElementsAreInitialized())
			mpScheme->InitializeElements(model);

        if(!mpConvergenceCriteria->mConvergenceCriteriaIsInitialized)
			mpConvergenceCriteria->Initialize(model);

		BaseType::mStiffnessMatrixIsBuilt = false;

        mInitializeWasPerformed = true;

		mIsConverged = false;
    }

    void InitializeSolutionStep()
    {
		ModelPart & model = BaseType::GetModelPart();

        mpBuilderAndSolver->ResizeAndInitializeVectors(mpScheme, 
					mpA, mpDx, mpb, 
					model);

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;
		
        mpBuilderAndSolver->InitializeSolutionStep(model, mA, mDx, mb);
        mpScheme->InitializeSolutionStep(model, mA, mDx, mb);

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		mpConvergenceCriteria->InitializeSolutionStep(model, rDofSet, mA, mDx, mb);
		
        mSolutionStepIsInitialized = true;

		mIsConverged = false;

		(*mDeltaUHat              ) = (*mConverged_DeltaUHat            );
		(*mDeltaUBar              ) = (*mConverged_DeltaUBar            );
		(*mDeltaU                 ) = (*mConverged_DeltaU               );
		(*mDeltaUStep             ) = (*mConverged_DeltaUStep           );
		(*mPHat                   ) = (*mConverged_PHat                 );
		mDeltaLambdaStep            = mConverged_DeltaLambdaStep         ;
		mCurrentLambda              = mConverged_CurrentLambda           ;
		mSignLastDeltaLambdaStep    = mConverged_SignLastDeltaLambdaStep ;
    }

	void FinalizeSolutionStep(TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
	{
		ModelPart & model = BaseType::GetModelPart();

		if (mCalculateReactionsFlag) {
            mpBuilderAndSolver->CalculateReactions(mpScheme, model, mA, mDx, mb);
		}

        mpScheme->FinalizeSolutionStep(model, mA, mDx, mb);
		mpScheme->Clean();
		
        mpBuilderAndSolver->FinalizeSolutionStep(model, mA, mDx, mb);

		if(mReformDofSetAtEachStep) 
			this->Clear();

		mSolutionStepIsInitialized = false;

		(*mConverged_DeltaUHat            ) = (*mDeltaUHat              );
		(*mConverged_DeltaUBar            ) = (*mDeltaUBar              );
		(*mConverged_DeltaU               ) = (*mDeltaU                 );
		(*mConverged_DeltaUStep           ) = (*mDeltaUStep             );
		(*mConverged_PHat                 ) = (*mPHat                   );
		mConverged_DeltaLambdaStep          = mDeltaLambdaStep           ;
		mConverged_CurrentLambda            = mCurrentLambda             ;
		mConverged_SignLastDeltaLambdaStep  = mSignLastDeltaLambdaStep   ;
	}
	
	void AbortSolutionStep()
	{
		ModelPart & model = BaseType::GetModelPart();

		if(BaseType::MoveMeshFlag())
		{
			for (ModelPart::NodeIterator i = model.NodesBegin(); i != model.NodesEnd(); ++i)
			{
				(i)->Coordinates() = (i)->GetInitialPosition() + (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
			}
		}

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
		{
			i_dof->GetSolutionStepValue() = i_dof->GetSolutionStepValue(1);
		}

		mSolutionStepIsInitialized = false;

		//mCurrentLambda = model.GetProcessInfo()[LAMBDA];
		(*mDeltaUHat              ) = (*mConverged_DeltaUHat              );
		(*mDeltaUBar              ) = (*mConverged_DeltaUBar              );
		(*mDeltaU                 ) = (*mConverged_DeltaU                 );
		(*mDeltaUStep             ) = (*mConverged_DeltaUStep             );
		(*mPHat                   ) = (*mConverged_PHat                   );
		mDeltaLambdaStep            = mConverged_DeltaLambdaStep        ;
		mCurrentLambda              = mConverged_CurrentLambda          ;
		mSignLastDeltaLambdaStep    = mConverged_SignLastDeltaLambdaStep;
	}
	
	/* 
	 * ==========================================================================
	 *                                utilities
	 * ==========================================================================
	 */

	double CheckIntegrationErrors()
	{
		ModelPart & model = BaseType::GetModelPart();
		std::vector<double> error_codes;
		for(ModelPart::ElementIterator elem_iter = model.ElementsBegin(); elem_iter != model.ElementsEnd(); elem_iter++)
		{
			Element& elem = *elem_iter;
			elem.GetValueOnIntegrationPoints(CONSTITUTIVE_INTEGRATION_ERROR_CODE, error_codes, model.GetProcessInfo());
			for(size_t i = 0; i < error_codes.size(); i++)
				if(error_codes[i] != 0.0)
					return error_codes[i];
		}
		return 0.0;
	}

	inline void PromptSeparator()
	{
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

	inline void PromptCurrentTime()
	{
		ModelPart& mp = BaseType::GetModelPart();
		const ProcessInfo& pinfo = mp.GetProcessInfo();
		std::stringstream ss;
		ss << std::endl;
		ss << "------------------------------------------------------------------------" << std::endl;
		ss << " TIME:  " << std::scientific << pinfo[TIME] << "  - DELTA TIME: " << pinfo[DELTA_TIME] << std::endl;
		ss << "------------------------------------------------------------------------" << std::endl;
		ss << std::endl;
		std::cout << ss.str();
	}

	inline void PromptMaxIterationsExceeded()
    {
		std::stringstream ss;
		ss << "------------------------------------------------------------------------" << std::endl;
		ss << " ERROR: max number of iterations exceeded" << std::endl;
		ss << "------------------------------------------------------------------------" << std::endl;
		std::cout << ss.str();
    }
	
	inline void PromptIntegrationErrors()
	{
		std::cout << "The equlibrium iteration process will be aborted because" 
				  << std::endl 
				  << " some elements reported integration errors" 
				  << std::endl;
	}

	inline void PromptArcLengthUpdateErrors()
	{
		std::cout << "The equlibrium iteration process will be aborted because" 
				  << std::endl 
				  << " the arc-length update procedure didn't converge" 
				  << std::endl;
	}

	int Check()
    {
        KRATOS_TRY

        BaseType::Check();

		ModelPart & model = BaseType::GetModelPart();

        mpBuilderAndSolver->Check(model);
        mpScheme->Check(model);
        mpConvergenceCriteria->Check(model);

        KRATOS_CATCH("")
		
		return 0;
    }

	/* 
	 * ==========================================================================
	 *                           arclength utilities
	 * ==========================================================================
	 */

#define ARCLEN_INIT_VEC(X,size) \
	if(X == NULL) { \
		TSystemVectorPointerType dummy = TSystemVectorPointerType(new TSystemVectorType(0)); \
		X.swap(dummy); \
	} \
	if (X->size() != size) X->resize(size, false); \
	X->clear();

	void ArcLength_Initialize()
	{
		// quick return
		if(mArcLengthInitialized) return;

		// initialize and resize auxiliary vectors
		ModelPart & model = BaseType::GetModelPart();
		unsigned int size = mpBuilderAndSolver->GetEquationSystemSize();
		
		ARCLEN_INIT_VEC(mDeltaUHat,  size);
		ARCLEN_INIT_VEC(mDeltaUBar,  size);
		ARCLEN_INIT_VEC(mDeltaU,     size);
		ARCLEN_INIT_VEC(mDeltaUStep, size);
		ARCLEN_INIT_VEC(mPHat,       size);

		ARCLEN_INIT_VEC(mConverged_DeltaUHat,  size);
		ARCLEN_INIT_VEC(mConverged_DeltaUBar,  size);
		ARCLEN_INIT_VEC(mConverged_DeltaU,     size);
		ARCLEN_INIT_VEC(mConverged_DeltaUStep, size);
		ARCLEN_INIT_VEC(mConverged_PHat,       size);

		mConverged_DeltaLambdaStep		   = 0.0;
		mConverged_CurrentLambda		   = 0.0;
		mConverged_SignLastDeltaLambdaStep = 1.0;

		// determination of PHat
		/*
		phat is the reference external force vector that we need to control with arc length.
		NOTE: it is required that before running the first call to this->Solve()
		      all external forces (to be controlled) should be set entirely (as if LAMBDA = 1)
			  in this way we can get PHat from the RHS assuming that at the end of the last
			  stage the residual is zero( or almost... numerically converged)
		*/
		TSystemVectorType& PHat = *mPHat;

		PHat.clear();
		mpBuilderAndSolver->BuildRHS(mpScheme, model, PHat); // now we have: RHS = Fext_n - R_n + dF_n+1 = dF_n+1
		mCurrentLambda = 0.0; // we assume that in python the LAMBDA is set to 1 on stage initialization!

		this->mConverged_PHat = this->mPHat;

		// check for a reference load
		bool haveLoad = false;
		for (unsigned int i=0; i<size; i++) {
			if ( PHat(i) != 0.0 ) {
				haveLoad = true;
				break;
			}
		}
		if (!haveLoad) {
			std::cout << "ArcLength Strategy:\nError: Zero external load vector to be piloted by the arc length\n";
			exit(-1);
		}

		// mark as initialized
		mArcLengthInitialized = true;
	}

	void ArcLength_InitializeSolutionStep()
	{
		ModelPart & model = BaseType::GetModelPart();

		double dtime = model.GetProcessInfo()[DELTA_TIME];
		if(mDeltaTime_0 == -1.0)
			mDeltaTime_0 = dtime;
		double new_arc_len = mArcLength_0 * dtime / mDeltaTime_0;
		if(mDeltaTime_0 == 0.0) {
			new_arc_len = mArcLength_0;
			mDeltaTime_0 = -1.0;
		}
		else {
			std::cout << "NEW ARC LENGTH ==================================================\n";
			KRATOS_WATCH(mDeltaTime_0);
			KRATOS_WATCH(new_arc_len);
			mArcLength2 = new_arc_len*new_arc_len;
		}

		if (mDeltaLambdaStep < 0.0)
			mSignLastDeltaLambdaStep = -1.0;
		else
			mSignLastDeltaLambdaStep = +1.0;
		//mSignLastDeltaLambdaStep = 1.0;

		// references
		TSystemMatrixType& mA     = *mpA;
		TSystemVectorType& dUHat  = *mDeltaUHat;
		TSystemVectorType& PHat   = *mPHat;
		TSystemVectorType& dU     = *mDeltaU;
		TSystemVectorType& dUStep = *mDeltaUStep;

		// determine mDeltaUHat
		mpBuilderAndSolver->BuildLHS(mpScheme, model, mA); // LHS in mA
		mpBuilderAndSolver->SystemSolve(mA, dUHat, PHat); // solve for dUHat using mA and PHat
		
		// determine delta lambda(1) == dlambda
		double dLambda = 0.0;
#if (ALEN_TYPE == ALEN_HSP)
		dLambda = std::sqrt(mArcLength2/( (mPsiU2/mURef2*inner_prod(dUHat,dUHat) ) + mPsiF2*inner_prod(PHat,PHat) ));
		dLambda *= mSignLastDeltaLambdaStep; // base sign of load change on what was happening last step
#elif (ALEN_TYPE == ALEN_CYL)
		dLambda = std::sqrt(mArcLength2/(inner_prod(dUHat,dUHat)+ALEN_CYL_ALPHA2));
		dLambda *= mSignLastDeltaLambdaStep; // base sign of load change on what was happening last step
#endif

		mDeltaLambdaStep = dLambda;
		mCurrentLambda += dLambda;

		// determine delta U(1) == dU
		noalias(dU) = dLambda*dUHat;
		dUStep = dU;
		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		mpScheme->Update(model, rDofSet, mA, dU, PHat);

		(*mpDx) = dU;
		if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
	}

	bool ArcLength_Update()
	{
		ModelPart & model = BaseType::GetModelPart();
		TSystemMatrixType& mA = *mpA;
		TSystemVectorType& mb = *mpb;
		TSystemVectorType& dX = *mpDx;
		TSystemVectorType& dUBar = *mDeltaUBar;
		TSystemVectorType& dUHat = *mDeltaUHat;
		TSystemVectorType& PHat = *mPHat;
		TSystemVectorType& dUStep = *mDeltaUStep;
		TSystemVectorType& dU = *mDeltaU;

		// copy
		dUBar = dX; // have to do this as the SOE is gonna change
		mpBuilderAndSolver->SystemSolve(mA, dUHat, PHat); // solve for dUHat using mA and PHat

#if (ALEN_TYPE == ALEN_HSP)
		
		// determine the coefficients of our quadratic equation
		double a1 = mPsiU2/mURef2*(inner_prod(dUHat,dUHat)) + mPsiF2 * inner_prod(PHat,PHat);
		double a2 = 2.0 *(mPsiU2/mURef2*(inner_prod(dUHat,dUBar)+inner_prod(dUHat,dUStep)) + mPsiF2*mDeltaLambdaStep * inner_prod(PHat,PHat));
		Vector aux(dUStep+dUBar);
		double a3 = mPsiU2/mURef2 * inner_prod(aux,aux) - mArcLength2 + (mDeltaLambdaStep*mDeltaLambdaStep)*mPsiF2 * inner_prod(PHat,PHat) ;

		// 935
		/*ModelPart::NodeType::Pointer cnode = model.pGetNode(935);
		unsigned int eid = cnode->GetDof(DISPLACEMENT_X).EquationId();
		double a1 = mPsiU2/mURef2*(dUHat[eid]*dUHat[eid]) + mPsiF2 * (PHat[eid]*PHat[eid]);
		double a2 = 2.0 *(mPsiU2/mURef2*((dUHat[eid]*dUBar[eid])+(dUHat[eid]*dUStep[eid])) + mPsiF2*mDeltaLambdaStep * (PHat[eid]*PHat[eid]));
		Vector aux(dUStep+dUBar);
		double a3 = mPsiU2/mURef2 * (aux[eid]*aux[eid]) - mArcLength2 + (mDeltaLambdaStep*mDeltaLambdaStep)*mPsiF2 * (PHat[eid]*PHat[eid]) ;*/

		// check for a solution to quadratic
		double b24ac = a2*a2 - a1*a3;
		if (b24ac < 0) {
			std::cout << "Arc Length Strategy: Error - Imaginary roots due to multiple instalibities\n";
			return false;
		}
		double dLambda;
		if (a1 == 0.0) {
			std::cout << "Arc Length Strategy: Error - zero denominator";
			dLambda = -a3/(2.0*a2);
		}
		else
		{
			// determine the roots of the quadratic
			double sqrtb24ac = sqrt(b24ac);
			double dlambda1 = (-a2 + sqrtb24ac)/a1;
			double dlambda2 = (-a2 - sqrtb24ac)/a1;

			// 1
			double val = inner_prod(dUHat,dUStep);
			double costheta1 = inner_prod(dUStep,dUStep) + inner_prod(dUBar,dUStep);
			double costheta2 = costheta1 + dlambda2*val;
			costheta1 += dlambda1*val;

			// 935
			/*double val = (dUHat[eid]*dUStep[eid]);
			double costheta1 = (dUStep[eid]*dUStep[eid]) + (dUBar[eid]*dUStep[eid]);
			double costheta2 = costheta1 + dlambda2*val;
			costheta1 += dlambda1*val;*/

			// choose dLambda based on angle between incremental displacement before
			// and after this step -- want positive
			if (costheta1 > costheta2)
				dLambda = dlambda1;
			else
				dLambda = dlambda2;
		}

#elif (ALEN_TYPE == ALEN_CYL)
		
		// determine the coeeficients of our quadratic equation
		double a = ALEN_CYL_ALPHA2 + inner_prod(dUHat,dUHat);
		double b = 2.0 * ( ALEN_CYL_ALPHA2*mDeltaLambdaStep + inner_prod(dUHat,dUBar) + inner_prod(dUStep,dUHat) );
		double c = 2.0 * inner_prod(dUStep,dUBar) + inner_prod(dUBar,dUBar);

		// check for a solution to quadratic
		double b24ac = b*b - 4.0*a*c;
		if (b24ac < 0.0) {
			std::cout << "ArcLength - imaginary roots due to multiple instability";
			std::cout << " directions - initial load increment was too large\n";
			std::cout << "a: " << a << " b: " << b << " c: " << c << " b24ac: " << b24ac << std::endl;
			return false;
		}
		double a2 = 2.0*a;
		if (a2 == 0.0) {
			std::cout << "ArcLength - zero denominator";
			std::cout << " alpha was set to 0.0 and zero reference load\n";
			return false;
		}

		// determine the roots of the quadratic
		double sqrtb24ac = std::sqrt(b24ac);
		double dlambda1 = (-b + sqrtb24ac)/a2;
		double dlambda2 = (-b - sqrtb24ac)/a2;
		
		double val = inner_prod(dUHat, dUStep);
		double theta1 = inner_prod(dUStep, dUStep) + inner_prod(dUBar, dUStep);
		//    double theta2 = theta1 + dlambda2*val;
		theta1 += dlambda1*val;
		
		// choose dLambda based on angle between incremental displacement before
		// and after this step -- want positive
		double dLambda;
		if (theta1 > 0.0)
			dLambda = dlambda1;
		else
			dLambda = dlambda2;

		// ALTERNATIVE TEST! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*ModelPart::NodeType& inode = model.GetNode(935);
		unsigned int idy = inode.GetDof(DISPLACEMENT_X).EquationId();
		array_1d<double,3> iu_n  = inode.FastGetSolutionStepValue(DISPLACEMENT,0);
		double uy1 = iu_n(1) + dUBar(idy) + dlambda1*dUHat(idy);
		double uy2 = iu_n(1) + dUBar(idy) + dlambda2*dUHat(idy);
		if(uy1 < uy2)
			dLambda = dlambda1;
		else
			dLambda = dlambda2;*/
		// ALTERNATIVE TEST! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#endif

		

		// determine delta U(i)
		dU = dUBar;
		dU += dLambda*dUHat;

		// update dU and dlambda
		dUStep += dU;
		mDeltaLambdaStep += dLambda;
		mCurrentLambda += dLambda;

		// update the model
		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		mpScheme->Update(model, rDofSet, mA, dU, mb);

		dX = dU;

		return true;
	}

	void LineSearchCalculation(double& rCurrentAlpha,double& rPreviousAlpha)
	{
		//return;
		ModelPart& r_model_part = BaseType::GetModelPart();
		TSystemMatrixType& A  = *mpA;
		TSystemVectorType& Dx = *mpDx;
		TSystemVectorType& b  = *mpb;
		TSystemVectorType& PHat = *mPHat;

		//LineSearchUtility<TSparseSpace, TDenseSpace, TLinearSolver> LineSearch;

		/*double ComputedAlpha = LineSearch.ExecuteLineSearch(
			mpBuilderAndSolver, 
			mpScheme, this->GetModelPart(), 
			mA, mDx, mb, 
			rCurrentAlpha, rPreviousAlpha); */

		double current_alpha = rCurrentAlpha;
		double previous_alpha = rPreviousAlpha;

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		TSystemVectorType  ReferenceDx = Dx;

		double R0= inner_prod(Dx,b);
		Dx = ReferenceDx *1.0;
		//mpScheme->Update(r_model_part,rDofSet,A,Dx,b);
		this->ArcLength_Update(); if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();

		mpBuilderAndSolver->BuildRHS(mpScheme, r_model_part, b);
		noalias(b) += mCurrentLambda * PHat; // Note: we cannot increase loads here in c++ (we do it in python)

		double R1= inner_prod(ReferenceDx,b);
		Dx *= (-1.0);
		//mpScheme->Update(r_model_part,rDofSet,A,Dx,b);
		this->ArcLength_Update(); if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
		b.clear();

		if(R0*R1<0){
			double R2 = R1;
			if(fabs(R1)<fabs(R0))
				R2=R0;
			double R0start = R0;
			double alpha = 0.0;
			double nabla = 0.0;
			double delta = 1.0;
			double CurrentAlpha  = 1.0; 
			int iterations=0;
			int max_iterations = 10;
			while(fabs(R2/R0start)>0.3 && iterations<max_iterations && (R1*R0)<0.0 && fabs(R1)>1.0e-7 && fabs(R0)>1.0e-7) {
				alpha = 0.5*(nabla+delta);
				CurrentAlpha  = alpha;
				Dx = ReferenceDx * CurrentAlpha;
				//mpScheme->Update(r_model_part,rDofSet,A,Dx,b);
				this->ArcLength_Update(); if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();

				mpBuilderAndSolver->BuildRHS(mpScheme, r_model_part, b);
				noalias(b) += mCurrentLambda * PHat; // Note: we cannot increase loads here in c++ (we do it in python)

				R2 = inner_prod(ReferenceDx,b);
				// ** Restore Current Displacement, Velocity, Acceleration
				Dx *= (-1.0);
				//mpScheme->Update(r_model_part,rDofSet,A,Dx,b);
				this->ArcLength_Update(); if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
				b.clear();
				if(R2*R1<0.0){
					//slope_k-1 = slope_k
					nabla = alpha;
					R0 = R2;
				}
				else if(R2*R0<0.0){
					//slope_k-1 = slope_k
					delta = alpha;
					R1 = R2;
				}
				else{
					break;
				}  
				iterations++;
			}
			std::cout<<" [ LINE SEARCH: (Iterations: "<<iterations<<", alpha: "<<CurrentAlpha<<") ] "<<std::endl;
			previous_alpha = current_alpha;
			current_alpha  = CurrentAlpha;
		}
		if(current_alpha>1 || current_alpha<=0)
			current_alpha=1;

		//Restore Current Displacement, Velocity, Acceleration
		Dx      = ReferenceDx;
		double ComputedAlpha = current_alpha;

		rPreviousAlpha = rCurrentAlpha;
		rCurrentAlpha  = ComputedAlpha;

		Dx *= rCurrentAlpha;

		rPreviousAlpha = rCurrentAlpha;
		rCurrentAlpha  = 1.0;
	}

protected:

    ArcLengthStrategy(const ArcLengthStrategy& Other)
    {
    }

protected:

    typename TSchemeType::Pointer                mpScheme;
    typename TLinearSolver::Pointer              mpLinearSolver;
    typename TBuilderAndSolverType::Pointer      mpBuilderAndSolver;
    typename TConvergenceCriteriaType::Pointer   mpConvergenceCriteria;

    TSystemVectorPointerType  mpDx;
    TSystemVectorPointerType  mpb;
    TSystemMatrixPointerType  mpA;

    bool            mReformDofSetAtEachStep;
	bool            mKeepSystemConstantDuringIterations;
	unsigned int    mMaxIterationNumber;
    bool            mCalculateReactionsFlag;
    bool            mInitializeWasPerformed;
	bool            mSolutionStepIsInitialized;
	bool            mIsConverged;
	
	bool            mAllowNonConvergence;

	double                   mArcLength2;
	double                   mArcLength_0;
	double                   mDeltaTime_0;
	double                   mPsiU2;
	double                   mPsiF2;
	double                   mURef2;
	TSystemVectorPointerType mDeltaUHat;
	TSystemVectorPointerType mDeltaUBar;
	TSystemVectorPointerType mDeltaU;
	TSystemVectorPointerType mDeltaUStep;
	TSystemVectorPointerType mPHat;
	double                   mDeltaLambdaStep;
	double                   mCurrentLambda;
	double                   mSignLastDeltaLambdaStep;
	bool                     mArcLengthInitialized;

	TSystemVectorPointerType mConverged_DeltaUHat;
	TSystemVectorPointerType mConverged_DeltaUBar;
	TSystemVectorPointerType mConverged_DeltaU;
	TSystemVectorPointerType mConverged_DeltaUStep;
	TSystemVectorPointerType mConverged_PHat;
	double                   mConverged_DeltaLambdaStep;
	double                   mConverged_CurrentLambda;
	double                   mConverged_SignLastDeltaLambdaStep;

}; /* Class ArcLengthStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_ARC_LENGTH_STRATEGY  defined */
