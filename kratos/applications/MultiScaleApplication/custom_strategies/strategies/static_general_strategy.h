//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $

#if !defined(KRATOS_STATIC_GENERAL_STRATEGY )
#define  KRATOS_STATIC_GENERAL_STRATEGY

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

namespace Kratos
{

/// Short class definition.
/**   Detail class definition.
 */
template<class TSparseSpace,
         class TDenseSpace, 
         class TLinearSolver
         >
class StaticGeneralStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(StaticGeneralStrategy);

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

    StaticGeneralStrategy(ModelPart& model_part,
                          typename TSchemeType::Pointer pScheme,
                          typename TLinearSolver::Pointer pNewLinearSolver,
                          typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                          int MaxIterations = 30,
                          bool CalculateReactions = false,
                          bool ReformDofSetAtEachStep = false,
                          bool MoveMeshFlag = false
						  )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
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

    StaticGeneralStrategy(ModelPart& model_part,
                          typename TSchemeType::Pointer pScheme,
                          typename TLinearSolver::Pointer pNewLinearSolver,
                          typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                          typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                          int MaxIterations = 30,
                          bool CalculateReactions = false,
                          bool ReformDofSetAtEachStep = false,
                          bool MoveMeshFlag = false
						  )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
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
        mpBuilderAndSolver = pNewBuilderAndSolver;

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

    virtual ~StaticGeneralStrategy()
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
        if (mInitializeWasPerformed == false) Initialize();
        if (mSolutionStepIsInitialized == false) InitializeSolutionStep();
        DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        mpScheme->Predict(BaseType::GetModelPart(), rDofSet, *mpA, *mpDx, *mpb);
        if (this->MoveMeshFlag()) BaseType::MoveMesh();
    }

	double Solve()
	{
		/* ----------------------------------------------------------
		* Perfomed operations:
		* - strategy initialization
		* - dofset/system initialization (if required)
		* - prediction
		* - initialization of the solution step
		*-----------------------------------------------------------*/

		mIsConverged = false;

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		ModelPart& model = BaseType::GetModelPart();

		if (!mInitializeWasPerformed) this->Initialize();

		model.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 0;

		// SUPPORT - 4
		if(model.GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] == 4)
		{
			model.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 1;
			this->FinalizeSolutionStep(*mpA, *mpDx, *mpb);
			return 0.0;
		}
		// SUPPORT - 4 - END

		if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true)
		{
			mpBuilderAndSolver->SetUpDofSet(mpScheme, model);
			mpBuilderAndSolver->SetUpSystem(model);
		}

		if (this->GetEchoLevel() != 0 && model.GetCommunicator().MyPID() == 0 ) PromptCurrentTime();

		this->Predict(); 

		if (!mSolutionStepIsInitialized) this->InitializeSolutionStep();

		TSystemMatrixType& mA	= *mpA;
		TSystemVectorType& mDx	= *mpDx;
		TSystemVectorType& mb	= *mpb;

		if(TSparseSpace::Size(mDx) == 0)
		{
			if (this->GetEchoLevel() != 0 && model.GetCommunicator().MyPID() == 0 )
				std::cout << "No Free DOFs - Finalizing solution step" << std::endl;
			mIsConverged = true;
		}
		else
		{
			/* ----------------------------------------------------------
			* Perfomed operations:
			* - Main Newton Raphson Iteration cicle
			*-----------------------------------------------------------*/

			unsigned int iteration_number = 0;
			model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

			while (true)
			{
				/* ----------------------------------------------------------
				* Perfomed operations:
				* - Initialize the nonlinear iteration
				*-----------------------------------------------------------*/

				iteration_number++;
				model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

				mpScheme->InitializeNonLinIteration(model, mA, mDx, mb);
				mIsConverged = mpConvergenceCriteria->PreCriteria(model, rDofSet, mA, mDx, mb);

				/* ----------------------------------------------------------
				* Perfomed operations:
				* - Build the system and solve
				*   0 -> build LHS just once
				*   1 -> build LHS at the beginning of each solution step
				*   2 -> build LHS at each iteration
				*-----------------------------------------------------------*/

				bool rebuildLHS = true;
				if(BaseType::mRebuildLevel == 0) 
					rebuildLHS = !BaseType::mStiffnessMatrixIsBuilt;
				else if(BaseType::mRebuildLevel == 1) 
					rebuildLHS = iteration_number == 1;

				TSparseSpace::SetToZero(mDx);
				TSparseSpace::SetToZero(mb);
				if(rebuildLHS)
				{
					TSparseSpace::SetToZero(mA);
					mpBuilderAndSolver->BuildAndSolve(mpScheme, model, mA, mDx, mb);
				}
				else
				{
					mpBuilderAndSolver->BuildRHSAndSolve(mpScheme, model, mA, mDx, mb);
				}
				
				/* ----------------------------------------------------------
				* Perfomed operations:
				* - Finalize the nonlinear iteration
				*-----------------------------------------------------------*/

				rDofSet = mpBuilderAndSolver->GetDofSet();
				mpScheme->Update(model, rDofSet, mA, mDx, mb);

				if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();

				mpScheme->FinalizeNonLinIteration(model, mA, mDx, mb);
				
				/* ----------------------------------------------------------
				* Perfomed operations:
				* - Convergence test
				*-----------------------------------------------------------*/

				/*if(iteration_number == 1)
				{
					if(this->GetEchoLevel() > 0 && model.GetCommunicator().MyPID() == 0)
						std::cout << "[INITIAL TEST]\n";
					mIsConverged = mpConvergenceCriteria->PostCriteria(model, rDofSet, mA, mDx, mb);
					TSparseSpace::SetToZero(mb);
					mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
				}*/
				mIsConverged = mpConvergenceCriteria->PostCriteria(model, rDofSet, mA, mDx, mb);

				if(mIsConverged || (iteration_number > mMaxIterationNumber) ) break;
			}
		}
		

		/* ----------------------------------------------------------
		* Perfomed operations:
		* - Perform the finalization of this solution step
		*   based on the STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL
		*   - 0 -> Always finalize this solution step (DEFAULT)
		*   - 1 -> Finalize the solution step 
		*          only if convergence is achieved, otherwise
		*          the solution step is aborted and we go back
		*          to the previous (converged) state
		*   - 2 -> Always abort this solution step
		*   - 3 -> Never finalize (nor abort) the solution step
		*          This is useful when the solution step has been
		*          used for small perturbations
		*-----------------------------------------------------------*/

		model.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 1;

		int finalize_level = model.GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL];
		bool verbose_allowed = model.GetCommunicator().MyPID() == 0;
		bool verbose = this->GetEchoLevel() > 0 && verbose_allowed;
		if(finalize_level == 0)
		{
			if(verbose && (mIsConverged == false)) PromptMaxIterationsExceeded();
			if(verbose) PromptSeparator();
			this->FinalizeSolutionStep(mA, mDx, mb);
			if(verbose) PromptSeparator();
		}
		else if(finalize_level == 1)
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
		else if(finalize_level == 2)
		{
			if(verbose && (mIsConverged == false)) PromptMaxIterationsExceeded();
			this->AbortSolutionStep();
		}
		else // finalize_level = 3
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

        mpBuilderAndSolver->ResizeAndInitializeVectors(
					mpA, mpDx, mpb, 
					model.Elements(), 
					model.Conditions(), 
					model.GetProcessInfo());

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        mpBuilderAndSolver->InitializeSolutionStep(model, mA, mDx, mb);
        mpScheme->InitializeSolutionStep(model, mA, mDx, mb);

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		mpConvergenceCriteria->InitializeSolutionStep(model, rDofSet, mA, mDx, mb);

        mSolutionStepIsInitialized = true;

		mIsConverged = false;
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
	}
	
	inline void PromptSeparator()
	{
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

	inline void PromptCurrentTime()
	{
		std::stringstream ss;
		ss << std::endl;
		ss << "------------------------------------------------------------------------" << std::endl;
		ss << " TIME:  " << std::scientific << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
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

    StaticGeneralStrategy(const StaticGeneralStrategy& Other)
    {
    }

}; /* Class StaticGeneralStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_STATIC_GENERAL_STRATEGY  defined */
