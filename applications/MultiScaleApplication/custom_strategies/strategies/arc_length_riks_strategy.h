//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $

#if !defined(KRATOS_ARCLENGTH_RIKS_STRATEGY )
#define  KRATOS_ARCLENGTH_RIKS_STRATEGY

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

#define ARCLEN_RIKS_INIT_VEC(X, newsize) \
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
class ArcLengthRiksStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(ArcLengthRiksStrategy);

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

    ArcLengthRiksStrategy(ModelPart& model_part,
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
		
		mFactor = 1.0;
		mCumulatedFactor = 1.0;
		mTotalFactor = 1.0;
		mDLambda_converged = 1.0;
    }

    ArcLengthRiksStrategy(ModelPart& model_part,
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
		
		mFactor = 1.0;
		mCumulatedFactor = 1.0;
		mTotalFactor = 1.0;
		mDLambda_converged = 1.0;
    }

    virtual ~ArcLengthRiksStrategy()
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
		//std::cout << "ALEN - prepare first solve - BEGIN\n";

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
			std::cout << "ArcLength Riks Strategy:\nError: Zero external load vector to be piloted by the arc length\n";
			exit(-1);
		}
		else
		{
			//std::cout << "ALEN - prepare first solve - calculating fhat\n";
			// init params
			mLambda = 1.0;
			mDLambda_converged = 1.0;
			mCycle = 0;
	
			// initialize arc length vectors
			ARCLEN_RIKS_INIT_VEC(mpFhat, sys_size);
			ARCLEN_RIKS_INIT_VEC(mpDu_converged, sys_size);
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
				std::cout << "ArcLength Strategy:\nError: Zero external load vector to be piloted by the arc length\n";
				exit(-1);
			}
		}
		
		// do not finalize this step. this is a pre-calculation to obtain Fhat!
		//mSolutionStepIsInitialized = false;
		this->FinalizeSolutionStep(*mpA, *mpDx, *mpb);

		//std::cout << "ALEN - prepare first solve - END\n";
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
		TSystemVectorType& mFhat         = *mpFhat; // reference external load to be piloted by the arc-length
		TSystemVectorType& mDu_converged = *mpDu_converged;
		
		unsigned int sys_size = TSparseSpace::Size(mDx);
		
		TSystemVectorType Du   (sys_size, 0.0);
		double current_lambda = mLambda;
		double Dlambda1 = 0.0;
		double Dlambda = 0.0;
		double current_factor = 1.0;
		//std::cout << "ALEN - solve - BEGIN\n";

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
			TSystemVectorType Fext (sys_size, 0.0);
			TSystemVectorType Du1  (sys_size, 0.0);
			TSystemVectorType ddu1 (sys_size, 0.0);
			TSystemVectorType ddu2 (sys_size, 0.0);
			TSystemVectorType ddu  (sys_size, 0.0);
			
			current_lambda = mLambda;
			KRATOS_WATCH(current_lambda);
			KRATOS_WATCH(mFactor);
			KRATOS_WATCH(mCumulatedFactor);
			KRATOS_WATCH(mFactor*mCumulatedFactor);
			
			/* 
			===============================================================================
			Begin the N-R equilibrium iteration loop
			===============================================================================
			*/

			unsigned int iteration_number = 0;
			model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

			if(mCycle == 1)
			{
				//std::cout << "ALEN - solve - First Cycle\n";
				//K,fint = assembleTangentStiffness( props, globdat )      
				//Da1    = globdat.dofs.solve( K , globdat.lam*fhat )
				//Dlam1  = globdat.lam
				TSparseSpace::SetToZero(Du1);
				TSparseSpace::SetToZero(mA);
				//mpBuilderAndSolver->BuildLHS(mpScheme, model, mA); // <<<<<<<<<<<<<<<<<<<<<< INT ERROR
				TSparseSpace::SetToZero(mb);
				mpBuilderAndSolver->Build(mpScheme, model, mA, mb); // <<<<<<<<<<<<<<<<<<<<<<<<<<<< INT ERROR
				noalias(Fext) = current_lambda*mFhat;
				mpBuilderAndSolver->SystemSolve(mA, Du1, Fext);
				Dlambda1 = current_lambda;
				//std::cout << "ALEN - solve - First Cycle End\n";
			}
			else
			{
				// Da1    = self.factor * self.Daprev
				// Dlam1  = self.factor * self.Dlamprev
				// globdat.lam += Dlam1
				//std::cout << "ALEN - solve - next Cycle\n";
				/*double aux = mTotalFactor*norm_2(mDu_converged);
				if(aux > 1.0) {
					mTotalFactor = 1.0;
				}*/
				current_factor = mFactor * mCumulatedFactor;
				noalias(Du1) = current_factor * mDu_converged;
				Dlambda1 = current_factor * mDLambda_converged;
				current_lambda += Dlambda1;
				//std::cout << "ALEN - solve - next Cycle end\n";
			}
			
			//KRATOS_WATCH(norm_2(Du1));

			// a [:] += Da1[:]
			// Da[:] =  Da1[:]
			// Dlam = Dlam1
			//std::cout << "ALEN - solve - prediction update\n";
			noalias(Du) = Du1;
			Dlambda = Dlambda1;
			rDofSet = mpBuilderAndSolver->GetDofSet();
			mpScheme->Update(model, rDofSet, mA, Du1, Fext);
			if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
			//std::cout << "ALEN - solve - prediction update - END\n";
			
			//std::cout << "ALEN - solve - prediction build\n";
			// K,fint = assembleTangentStiffness( props, globdat )  
			// res = globdat.lam*fhat-fint
			TSparseSpace::SetToZero(mb);
			TSparseSpace::SetToZero(mA);
			mpBuilderAndSolver->Build(mpScheme, model, mA, mb); // <<<<<<<<<<<<<<<<<<<<<<<<<<<< INT ERROR
			noalias(Fext) = current_lambda*mFhat;
			noalias(mb) += Fext;
			//std::cout << "ALEN - solve - prediction build - END\n";
			
			while (true)
			{
				/* 
				===============================================================================
				Initialize the current nonlinear iteration
				===============================================================================
				*/

				iteration_number++;
				//std::cout << "ALEN - solve - iteratiion: " << iteration_number << "\n";
				model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

				mpScheme->InitializeNonLinIteration(model, mA, mDx, mb);
				mIsConverged = mpConvergenceCriteria->PreCriteria(model, rDofSet, mA, mDx, mb); // occhio : usare residual norm criteria!

				//std::cout << "ALEN - solve - begin double solve\n";
				// d1 = globdat.dofs.solve( K , fhat )
				// d2 = globdat.dofs.solve( K , res )
				TSparseSpace::SetToZero(ddu1);
				TSparseSpace::SetToZero(ddu2);
				mpBuilderAndSolver->SystemSolve(mA, ddu1, mFhat);
				mpBuilderAndSolver->SystemSolve(mA, ddu2, mb);
				//std::cout << "ALEN - solve - end double solve\n";
				
				// ddlam = -dot(Da1,d2)/dot(Da1,d1)
				// dda   = ddlam*d1 + d2
				double ddlam = -inner_prod(Du1,ddu2)/inner_prod(Du1,ddu1);
				noalias(ddu) = ddlam*ddu1 + ddu2;
	
				// Dlam        += ddlam
				// globdat.lam += ddlam
				Dlambda        += ddlam;
				current_lambda += ddlam;
				
				// Da[:] += dda[:]
				// a [:] += dda[:]
				//std::cout << "ALEN - solve - update with ddu\n";
				noalias(Du) += ddu;
				rDofSet = mpBuilderAndSolver->GetDofSet();
				mpScheme->Update(model, rDofSet, mA, ddu, mb);
				if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
				mpScheme->FinalizeNonLinIteration(model, mA, ddu, mb);
				//std::cout << "ALEN - solve - update with ddu end\n";
				
				/* 
				===============================================================================
				Build system matrices and vectors, abort if errors are detected
				===============================================================================
				*/
				
				// K,fint = assembleTangentStiffness( props, globdat )
				// res = globdat.lam*fhat-fint
				// error  = globdat.dofs.norm( res ) / globdat.dofs.norm( globdat.lam*fhat )
				//std::cout << "ALEN - build new\n";
				TSparseSpace::SetToZero(mb);
				TSparseSpace::SetToZero(mA);
				mpBuilderAndSolver->Build(mpScheme, model, mA, mb);
				if(this->CheckIntegrationErrors() != 0.0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				noalias(Fext) = current_lambda*mFhat;
				noalias(mb) += Fext;
				//std::cout << "ALEN - build new END\n";
				
				/* 
				===============================================================================
				Check convergence criteria.
				Also check for divergence and if so, abort solution step.
				===============================================================================
				*/

				model.GetProcessInfo()[ITERATION_CONVERGENCE_FLAG] = 0;

				mIsConverged = mpConvergenceCriteria->PostCriteria(model, rDofSet, mA, ddu, mb);

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
		
		//std::cout << "ALEN - solve converged\n";

		if(mIsConverged)
		{
			mLambda = current_lambda;
			mDLambda_converged = Dlambda;
			noalias(mDu_converged) = Du;
			mCumulatedFactor = 1.0;
		}
		else
		{
			mCumulatedFactor = current_factor;
		}
		
		model.GetProcessInfo()[LAMBDA_OUTPUT] = mLambda;

		//std::cout << "ALEN - solve END\n";
		
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

    ArcLengthRiksStrategy(const ArcLengthRiksStrategy& Other)
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
	TSystemVectorPointerType mpDu_converged;
	double mLambda;
	double mDLambda_converged;
	double mFactor;
	double mCumulatedFactor;
	double mTotalFactor;
	unsigned int mCycle;

    bool            mReformDofSetAtEachStep;
	bool            mKeepSystemConstantDuringIterations;
	unsigned int    mMaxIterationNumber;
    bool            mCalculateReactionsFlag;
    bool            mInitializeWasPerformed;
	bool            mSolutionStepIsInitialized;
	bool            mIsConverged;
	
	bool            mAllowNonConvergence;

}; /* Class ArcLengthRiksStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_ARCLENGTH_RIKS_STRATEGY  defined */
