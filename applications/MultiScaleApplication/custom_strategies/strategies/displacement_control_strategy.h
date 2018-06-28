//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $

#if !defined(KRATOS_DISPLACEMENT_CONTROL_STRATEGY )
#define  KRATOS_DISPLACEMENT_CONTROL_STRATEGY

/* External includes */
#include "boost/smart_ptr.hpp"
#include <iomanip>
#include <boost/math/special_functions/fpclassify.hpp>

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "custom_strategies/builder_and_solvers/static_general_builder_and_solver.h"

#define DISPCTRL_INIT_VEC(X, newsize) \
	if(X == NULL) { \
		TSystemVectorPointerType dummy = TSystemVectorPointerType(new TSystemVectorType(0)); \
		X.swap(dummy); \
	} \
	if (X->size() != newsize) X->resize(newsize, false); \
	X->clear();

namespace Kratos
{

/// Short class definition.
/**   Detail class definition.
 */
template<class TSparseSpace,
         class TDenseSpace, 
         class TLinearSolver
         >
class DisplacementControlStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(DisplacementControlStrategy);

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

    DisplacementControlStrategy(ModelPart& model_part,
                          typename TSchemeType::Pointer pScheme,
                          typename TLinearSolver::Pointer pNewLinearSolver,
                          typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                          int MaxIterations = 30,
                          bool CalculateReactions = false,
                          bool ReformDofSetAtEachStep = false,
                          bool MoveMeshFlag = false,
						  unsigned int ControlNodeID = 0,
						  unsigned int ControlDofID = 0,
						  double TargetDisplacement = 0
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
		
		mFactor = 1.0;
		mTotalFactor = 1.0;
		mControlNodeID = ControlNodeID;
		mControlDofID = ControlDofID;
		mTargetDisplacement = TargetDisplacement;
		m_U_q_initial = 0.0;
		m_dU_bar_q = 0.0;

		mFirstStepDone = false;
    }

    DisplacementControlStrategy(ModelPart& model_part,
                          typename TSchemeType::Pointer pScheme,
                          typename TLinearSolver::Pointer pNewLinearSolver,
                          typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                          typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                          int MaxIterations = 30,
                          bool CalculateReactions = false,
                          bool ReformDofSetAtEachStep = false,
                          bool MoveMeshFlag = false,
						  unsigned int ControlNodeID = 0,
						  unsigned int ControlDofID = 0,
						  double TargetDisplacement = 0
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
		
		mFactor = 1.0;
		mTotalFactor = 1.0;
		mControlNodeID = ControlNodeID;
		mControlDofID = ControlDofID;
		mTargetDisplacement = TargetDisplacement;
		m_U_q_initial = 0.0;
		m_dU_bar_q = 0.0;

		mFirstStepDone = false;
    }

    virtual ~DisplacementControlStrategy()
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
	
	void SetLoadFactors(double factor, double totalFactor)
	{
		mFactor = factor;
		mTotalFactor = totalFactor;
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

	void PrepareFirstSolve()
	{
		mIsConverged = false;

		ModelPart& model = BaseType::GetModelPart();

		if (!mInitializeWasPerformed) this->Initialize();

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

		if (!mSolutionStepIsInitialized) this->InitializeSolutionStep();

		TSystemVectorType& mDx	= *mpDx;
		unsigned int sys_size = TSparseSpace::Size(mDx);
		
		if(sys_size == 0)
		{
			std::cout << "Displacement Control Strategy:\nError: Zero external load vector to be piloted by the strategy\n";
			exit(-1);
		}
		else
		{
			// init params
			mLambda = 0.0;
			mCycle = 0;
			
			// check controlled node
			ModelPart::NodeType::Pointer cnode = model.pGetNode(mControlNodeID);
			bool has_control_node = false;
			if(cnode) {
				if(cnode->SolutionStepsDataHas(DISPLACEMENT)) {
					if(mControlDofID == 0) {
						if(cnode->HasDofFor(DISPLACEMENT_X)) 
							if(cnode->pGetDof(DISPLACEMENT_X)->IsFree())
								has_control_node = true;
					}
					else if(mControlDofID == 1) {
						if(cnode->HasDofFor(DISPLACEMENT_Y)) 
							if(cnode->pGetDof(DISPLACEMENT_Y)->IsFree())
								has_control_node = true;
					}
					else if(mControlDofID == 2) {
						if(cnode->HasDofFor(DISPLACEMENT_Z)) 
							if(cnode->pGetDof(DISPLACEMENT_Z)->IsFree())
								has_control_node = true;
					}
				}
			}
			if(!has_control_node) {
				std::cout << "Displacement Control Strategy:\nError: Control Node is not existent\n" << std::endl;
				exit(-1);
			}
			
			// get initial displacement on controlled node
			m_U_q_initial = 0.0;
			if(mControlDofID == 0)
				m_U_q_initial = cnode->FastGetSolutionStepValue(DISPLACEMENT_X);
			else if(mControlDofID == 1) 
				m_U_q_initial = cnode->FastGetSolutionStepValue(DISPLACEMENT_Y);
			else 
				m_U_q_initial = cnode->FastGetSolutionStepValue(DISPLACEMENT_Z);
			KRATOS_WATCH(m_U_q_initial);
			
			// initialize arc length vectors
			DISPCTRL_INIT_VEC(mpFhat, sys_size);
			// determine Fhat
			// assuming that this is the first stage, or that the previous stage has converged
			// the RHS vector will contain only the external forces and 0 internal forces (balanced without external forces)
			// Here we further assume that LAMBDA is set to 1
			TSystemVectorType& mFhat = *mpFhat;
			mpBuilderAndSolver->BuildRHS(mpScheme, model, mFhat);
			// check for a reference load
			bool haveLoad = false;
			for (unsigned int i=0; i<sys_size; i++) {
				if ( mFhat(i) != 0.0 ) {
					haveLoad = true;
					break;
				}
			}
			if (!haveLoad) {
				std::cout << "Displacement Control Strategy:\nError: Zero external load vector to be piloted by the strategy\n";
				exit(-1);
			}
		}
		
		// do not finalize this step. this is a pre-calculation to obtain Fhat!
		//mSolutionStepIsInitialized = false;
		this->FinalizeSolutionStep(*mpA, *mpDx, *mpb);
		mFirstStepDone = false;
	}
	
	double Solve()
	{
		mCycle++; // increment cycle counter
		
		/* 
		===============================================================================
		Solution step initialization
		===============================================================================
		*/
		
		mIsConverged = false;

		ModelPart& model = BaseType::GetModelPart();

		if (!mInitializeWasPerformed) this->Initialize();

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

		if (!mSolutionStepIsInitialized) this->InitializeSolutionStep();

		TSystemMatrixType& mA	         = *mpA;
		TSystemVectorType& mDx	         = *mpDx;
		TSystemVectorType& mb	         = *mpb;
		TSystemVectorType& mFhat         = *mpFhat; // reference external load to be piloted by the strategy
		
		unsigned int sys_size = TSparseSpace::Size(mDx);
		
		// control node
		ModelPart::NodeType::Pointer cnode = model.pGetNode(mControlNodeID);
		// control equation id
		unsigned int cequid = 0;
		if(mControlDofID == 0)
			cequid = cnode->pGetDof(DISPLACEMENT_X)->EquationId();
		else if(mControlDofID == 1) 
			cequid = cnode->pGetDof(DISPLACEMENT_Y)->EquationId();
		else 
			cequid = cnode->pGetDof(DISPLACEMENT_Z)->EquationId();
		
		// calculate initial imposed target displacement increment
		if(mCycle == 1)
		{
			double time_beging = model.GetProcessInfo()[TIME];
			double time_end = model.GetProcessInfo()[END_TIME];
			double delta_time = model.GetProcessInfo()[DELTA_TIME];
			double duration = time_end-time_beging+delta_time;
			m_dU_bar_q = (mTargetDisplacement - m_U_q_initial) * delta_time / duration;
			KRATOS_WATCH(mTargetDisplacement);
			KRATOS_WATCH(m_dU_bar_q);
			KRATOS_WATCH(delta_time);
			KRATOS_WATCH(duration);
		}
		
		// lambda
		double current_lambda = mLambda;
		
		if(sys_size == 0)
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
			TSystemVectorType dFext (sys_size, 0.0);
			TSystemVectorType Uhat (sys_size, 0.0);
			TSystemVectorType Uaux (sys_size, 0.0);
			
			/* 
			===============================================================================
			Begin the N-R equilibrium iteration loop
			===============================================================================
			*/

			unsigned int iteration_number = 0;
			model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

			// calculate the current stiffness matrix
			TSparseSpace::SetToZero(mb);
			if(mReformDofSetAtEachStep || !mFirstStepDone)
			{
				std::cout << "FORMING PREVIOUS STIFFNESS\n";
				// we don't have a previous converged stiffness
				TSparseSpace::SetToZero(mA);
				mpBuilderAndSolver->Build(mpScheme, model, mA, mb); // <<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK INT ERROR
			}// use previous converged stiffness
			
			// solve for the tangential displacement due to the reference (total) load
			mpBuilderAndSolver->SystemSolve(mA, Uhat, mFhat);
			double U_hat_pilot_q = Uhat(cequid);
			KRATOS_WATCH(U_hat_pilot_q);
			
			// update the load scale factor
			double delta_lambda = m_dU_bar_q*mTotalFactor / U_hat_pilot_q;
			current_lambda += delta_lambda;
			KRATOS_WATCH(mTotalFactor);
			KRATOS_WATCH(delta_lambda);
			KRATOS_WATCH(current_lambda);
			
			// calculate the incremental external force vector
			noalias(dFext) = delta_lambda * mFhat;
			
			// Solve for the displacement increment, and update
			TSparseSpace::SetToZero(mDx);
			mpBuilderAndSolver->SystemSolve(mA, mDx, dFext);
			rDofSet = mpBuilderAndSolver->GetDofSet();
			mpScheme->Update(model, rDofSet, mA, mDx, mb);
			if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
			
			// Calculate the internal force vector, and the residual
			TSparseSpace::SetToZero(mb);
			mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
			noalias(mb) += current_lambda * mFhat;
			
			while (true)
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
				
				// update the stiffness matrix
				TSparseSpace::SetToZero(mb);
				// ----------------------------------------------------------
				// REBUILD LHS ALWAYS
				/*TSparseSpace::SetToZero(mA);
				mpBuilderAndSolver->Build(mpScheme, model, mA, mb);*/
				// REBUILD LHS ONLY ONCE PER STEP
				if(iteration_number == 1)
				{
					TSparseSpace::SetToZero(mA);
					mpBuilderAndSolver->Build(mpScheme, model, mA, mb);
				}
				else
				{
					mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
				}
				// ----------------------------------------------------------
				noalias(mb) += current_lambda * mFhat;
				if(this->CheckIntegrationErrors() != 0.0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				
				// update the load ratio corrector
				mpBuilderAndSolver->SystemSolve(mA, Uaux, mb);
				// check solution validity
				if(!this->CheckSolutionValidity(Uaux))
				{
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				double U_tilde = Uaux(cequid);
				mpBuilderAndSolver->SystemSolve(mA, Uaux, mFhat);
				// check solution validity
				if(!this->CheckSolutionValidity(Uaux))
				{
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				double U_hat = Uaux(cequid);
				double d_lambda = U_tilde / U_hat;
				current_lambda -= d_lambda;
				
				// modify the old residual vector
				noalias(mb) -= d_lambda * mFhat;
				
				// solve for the displacement correction and update
				mpBuilderAndSolver->SystemSolve(mA, mDx, mb);
				// check solution validity
				if(!this->CheckSolutionValidity(mDx))
				{
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				mpScheme->Update(model, rDofSet, mA, mDx, mb);
				if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
				mpScheme->FinalizeNonLinIteration(model, mA, mDx, mb);
				
				// compute the new internal nodal forces and residual vector
				TSparseSpace::SetToZero(mb);
				mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
				noalias(mb) += current_lambda * mFhat;
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

				if(mIsConverged && iteration_number == 1)
					mIsConverged = false;

				if(mIsConverged || (iteration_number > mMaxIterationNumber) ) break;
			}

			// refesh RHS (we DON'T do this in BuilderAndSolver to check for integration errors!!!!)
			TSparseSpace::SetToZero(mb);
			mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
			noalias(mb) += current_lambda * mFhat;
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

		/* 
		===============================================================================
		save data
		===============================================================================
		*/
		
		if(mIsConverged)
		{
			mLambda = current_lambda;
		}
		
		model.GetProcessInfo()[LAMBDA_OUTPUT] = mLambda;

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

		mFirstStepDone = true;
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

	bool CheckSolutionValidity(TSystemVectorType& mDx)
	{
		bool is_valid = true;
		unsigned int eqsize = TSparseSpace::Size(mDx);
		for(unsigned int i = 0; i < eqsize; i++) {
			double idx = mDx[i];
			if(!( (boost::math::isfinite)(idx) )) {
				std::cout << "found a non finite value in the solution vector (nan or inifinite)\n";
				is_valid = false;
				break;
			}
		}
		return is_valid;
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

protected:

    DisplacementControlStrategy(const DisplacementControlStrategy& Other)
    {
    }

protected:

    typename TSchemeType::Pointer                mpScheme;
    typename TLinearSolver::Pointer              mpLinearSolver;
    typename TBuilderAndSolverType::Pointer      mpBuilderAndSolver;
    typename TConvergenceCriteriaType::Pointer   mpConvergenceCriteria;

    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb;
    TSystemMatrixPointerType mpA;
	
	TSystemVectorPointerType mpFhat;
	double mLambda;
	double mFactor;
	double mTotalFactor;
	unsigned int mCycle;

	unsigned int mControlNodeID;
	unsigned int mControlDofID; // 0=X, 1=Y, 2=Z
	double mTargetDisplacement;
	double m_U_q_initial;
	double m_dU_bar_q;
	
    bool            mReformDofSetAtEachStep;
	bool            mKeepSystemConstantDuringIterations;
	unsigned int    mMaxIterationNumber;
    bool            mCalculateReactionsFlag;
    bool            mInitializeWasPerformed;
	bool            mSolutionStepIsInitialized;
	bool            mIsConverged;
	
	bool            mAllowNonConvergence;

	bool            mFirstStepDone;

}; /* Class DisplacementControlStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_DISPLACEMENT_CONTROL_STRATEGY  defined */
