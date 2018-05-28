//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $

#if !defined(KRATOS_ARCLENGTH_DISS_NRG_STRATEGY )
#define  KRATOS_ARCLENGTH_DISS_NRG_STRATEGY

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

#define ARCLEN_DISS_NRG_INIT_VEC(X, newsize) \
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
class ArcLengthDissNrgStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(ArcLengthDissNrgStrategy);

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

    ArcLengthDissNrgStrategy(ModelPart& model_part,
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
		
		// initialize arc length parameters
		mLambda = 0.0;
		mDLambda = 1.0;
		mdTau = 0.0;
		mFactor = 1.0;
		mCumulatedFactor = 1.0;
		mTotalFactor = 1.0;
		mIsNrgControlled = false; // start with force-control and then switch to dissipated energy control
		mDissipatedEnergy = 0.0;
    }

    ArcLengthDissNrgStrategy(ModelPart& model_part,
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
		
		// initialize arc length parameters
		mLambda = 0.0;
		mDLambda = 1.0;
		mdTau = 0.0;
		mFactor = 1.0;
		mCumulatedFactor = 1.0;
		mTotalFactor = 1.0;
		mIsNrgControlled = false; // start with force-control and then switch to dissipated energy control
		mDissipatedEnergy = 0.0;
    }

    virtual ~ArcLengthDissNrgStrategy()
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

		ModelPart& model_part = BaseType::GetModelPart();

		if (!mInitializeWasPerformed) this->Initialize();

		bool verbose_allowed = model_part.GetCommunicator().MyPID() == 0;
		bool verbose         = this->GetEchoLevel() > 0 && verbose_allowed;

		model_part.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 0;

		if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true)
		{
			mpBuilderAndSolver->SetUpDofSet(mpScheme, model_part);
			mpBuilderAndSolver->SetUpSystem(model_part);
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
	
			// initialize arc length vectors
			ARCLEN_DISS_NRG_INIT_VEC(mpFhat, sys_size);
			// determine Fhat
			// assuming that this is the first stage, or that the previous stage has converged
			// the RHS vector will contain only the external forces and 0 internal forces (balanced without external forces)
			// Here we further assume that LAMBDA is set to 1
			TSystemVectorType& mFhat = *mpFhat;
			mpBuilderAndSolver->BuildRHS(mpScheme, model_part, mFhat);
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
		
		/* 
		===============================================================================
		Solution step initialization
		===============================================================================
		*/
		
		mIsConverged = false;

		ModelPart& model_part = BaseType::GetModelPart();

		if (!mInitializeWasPerformed) this->Initialize();

		bool verbose_allowed = model_part.GetCommunicator().MyPID() == 0;
		bool verbose         = this->GetEchoLevel() > 0 && verbose_allowed;

		model_part.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 0;

		if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true)
		{
			mpBuilderAndSolver->SetUpDofSet(mpScheme, model_part);
			mpBuilderAndSolver->SetUpSystem(model_part);
		}

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();

		if (verbose) PromptCurrentTime();

		this->Predict(); 

		if (!mSolutionStepIsInitialized) this->InitializeSolutionStep();

		TSystemMatrixType& mA	         = *mpA;
		TSystemVectorType& mDx	         = *mpDx;
		TSystemVectorType& mb	         = *mpb;
		TSystemVectorType& mFhat         = *mpFhat; // reference external load to be piloted by the arc-length
		
		unsigned int sys_size = TSparseSpace::Size(mDx);
		
		// arc length variables
		double lambda  = mLambda;
		double lambda0 = mLambda;
		double Dlambda = mDLambda;
		double dTau    = mdTau;
		double current_factor = 1.0;

		if(!mIsNrgControlled)
		{
			Dlambda *= mTotalFactor;
			lambda += Dlambda;
		}
		else
		{
			current_factor = mFactor * mCumulatedFactor;
			dTau *= current_factor;
		}
		//dTau = mdTau*mTotalFactor;

		TSystemVectorType resi (sys_size, 0.0); // residual vector 
		
		if(sys_size == 0)
		{
			/* 
			===============================================================================
			If all DOFs are set, just update the RHS and move on
			===============================================================================
			*/

			mpBuilderAndSolver->BuildRHS(mpScheme, model_part, mb);
			if(this->CheckIntegrationErrors() != 0.0)
			{
				if(verbose) this->PromptIntegrationErrors();
				mIsConverged = false;
				model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = mMaxIterationNumber;
			}
			else
			{
				if (verbose) std::cout << "No Free DOFs - Finalizing solution step" << std::endl;
				mIsConverged = true;
			}
		}
		else
		{
			// allocate extra vectors
			TSystemVectorType Fext (sys_size, 0.0);
			TSystemVectorType u    (sys_size, 0.0);
			TSystemVectorType Du   (sys_size, 0.0);
			TSystemVectorType ddu1 (sys_size, 0.0);
			TSystemVectorType ddu2 (sys_size, 0.0);
			TSystemVectorType h    (sys_size, 0.0);
			TSystemVectorType u_Du (sys_size, 0.0);
			TSystemVectorType aux  (sys_size, 0.0);
			
			// extract the current displacement vector
			for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
			{
				if(i_dof->IsFree())
				{
					//i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
					u[i_dof->EquationId()] = i_dof->GetSolutionStepValue();
				}
			}

			/* 
			===============================================================================
			Begin the N-R equilibrium iteration loop
			===============================================================================
			*/

			unsigned int iteration_number = 0;
			model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

			// K,fint = assembleTangentStiffness( props, globdat )
			TSparseSpace::SetToZero(mDx);
			TSparseSpace::SetToZero(mA);
			TSparseSpace::SetToZero(mb);
			mpBuilderAndSolver->Build(mpScheme, model_part, mA, mb); // <<<<<<<<<<<<<<<<<<<<<<<<<<<< INT ERROR
			
			while (true)
			{
				/* 
				===============================================================================
				Initialize the current nonlinear iteration
				===============================================================================
				*/

				iteration_number++;
				model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

				mpScheme->InitializeNonLinIteration(model_part, mA, mDx, mb);
				mIsConverged = mpConvergenceCriteria->PreCriteria(model_part, rDofSet, mA, mDx, mb);

				if(!mIsNrgControlled)
				{
					// FORCE-CONTROL
					// da = globdat.dofs.solve( K, globdat.lam*fhat-fint )
					KRATOS_WATCH(lambda);
					KRATOS_WATCH(Dlambda);
					KRATOS_WATCH(mFactor);
					KRATOS_WATCH(mTotalFactor);
					noalias(resi) = lambda*mFhat + mb;
					mpBuilderAndSolver->SystemSolve(mA, mDx, resi);
				}
				else
				{
					// DISSIPATED ENERGY CONTROL
					/*
					h  =  0.5 * lam0 * fhat
					w  = -0.5 * dot ( (a-Da) , fhat )
					g  =  0.5 * dot ( ( lam0 * Da - self.Dlam * ( a[:] - Da[:] ) ) , fhat ) - globdat.dtau
			  
					d1 = globdat.dofs.solve( K , globdat.lam*fhat - fint )
					d2 = globdat.dofs.solve( K , -1.0*fhat )

					denom  = dot ( h , d2 ) - w

					da     = d1 - ( d2 * ( dot( h , d1 ) + g ) ) / denom
					dlam   = -g - ( dot( -1.0*h , d1 ) - g * ( 1.0 + denom ) ) / denom;
			 
					self.Dlam   += dlam
					globdat.lam += dlam
					*/

					KRATOS_WATCH(lambda);
					KRATOS_WATCH(Dlambda);
					KRATOS_WATCH(dTau);
					KRATOS_WATCH(mFactor);
					KRATOS_WATCH(mTotalFactor);

					noalias(u_Du) = u - Du;
					noalias(aux)  = lambda0 * Du;
					noalias(aux) -= Dlambda * u_Du;
					
					noalias(h) =  0.5 * lambda0 * mFhat;
					double w   = -0.5 * inner_prod( u_Du, mFhat );
					double g   =  0.5 * inner_prod(  aux, mFhat ) - dTau;

					noalias(resi) = lambda*mFhat + mb;
					noalias(aux)  = -mFhat;
					TSparseSpace::SetToZero(ddu1);
					TSparseSpace::SetToZero(ddu2);
					mpBuilderAndSolver->SystemSolve(mA, ddu1, resi);
					mpBuilderAndSolver->SystemSolve(mA, ddu2, aux);
					
					double denom  = inner_prod( h , ddu2 ) - w;
					
					noalias(mDx) = ddu1 - ( ddu2 * ( inner_prod(h,ddu1) + g ) ) / denom;
					double dlam  = -g - ( inner_prod( -h, ddu1 ) - g * ( 1.0 + denom ) ) / denom;
					
					Dlambda += dlam;
					lambda  += dlam;
				}
				
				// update displacements
				// Da[:] += da[:]
				// a [:] += da[:]
				noalias(u)  += mDx;
				noalias(Du) += mDx;
				rDofSet = mpBuilderAndSolver->GetDofSet();
				mpScheme->Update(model_part, rDofSet, mA, mDx, resi);
				if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
				mpScheme->FinalizeNonLinIteration(model_part, mA, mDx, mb);
				
				
				// K,fint = assembleTangentStiffness( props, globdat )
				// error  = globdat.dofs.norm( globdat.lam*fhat-fint ) / globdat.dofs.norm( globdat.lam*fhat )
				TSparseSpace::SetToZero(mb);
				TSparseSpace::SetToZero(mA);
				mpBuilderAndSolver->Build(mpScheme, model_part, mA, mb);
				if(this->CheckIntegrationErrors() != 0.0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				noalias(resi) = lambda*mFhat + mb;
				
				
				/* 
				===============================================================================
				Check convergence criteria.
				Also check for divergence and if so, abort solution step.
				===============================================================================
				*/

				model_part.GetProcessInfo()[ITERATION_CONVERGENCE_FLAG] = 0;

				mIsConverged = mpConvergenceCriteria->PostCriteria(model_part, rDofSet, mA, mDx, resi);
				if(iteration_number == 1)
					mIsConverged = false;

				if(model_part.GetProcessInfo()[ITERATION_CONVERGENCE_FLAG] != 0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}

				if(mIsConverged || (iteration_number > mMaxIterationNumber) ) break;
			}

			// refesh RHS (we DON'T do this in BuilderAndSolver to check for integration errors!!!!)
			TSparseSpace::SetToZero(mb);
			mpBuilderAndSolver->BuildRHS(mpScheme, model_part, mb);
			if(this->CheckIntegrationErrors() != 0.0)
			{
				if(verbose) this->PromptIntegrationErrors();
				mIsConverged = false;
				iteration_number = mMaxIterationNumber;
				model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
			}
			noalias(resi) = lambda*mFhat + mb;
			
			// If converged, calculate the amount of energy that has been dissipated in the previous step.
			// dissnrg = 0.5 * dot( ( lam0 * Da - (globdat.lam-lam0) * ( a - Da ) ),fhat )
			if(mIsConverged)
			{
				noalias(u_Du) = u - Du;
				noalias(aux)  = lambda0 * Du;
				noalias(aux) -= (lambda-lambda0) * u_Du;
				double dissnrg = 0.5 * inner_prod( aux, mFhat );
				std::cout << "DISSIPATED ENERGY = " << dissnrg << std::endl;
				mDissipatedEnergy = dissnrg;

				/*
				if self.method == 'force-controlled':
				  if dissnrg > self.switchEnergy:
					print '   Switch to nrg diss. arc-length'
					self.method       = 'nrg-controlled'
					globdat.dtau = 0.25*self.switchEnergy
				  else:
					globdat.lam += self.Dlam
				else:
				  self.Dlam = 0.
				  globdat.dtau *= pow(0.5,0.25*(globdat.iiter-self.optiter))
				  if globdat.dtau > self.maxdTau:
					globdat.dtau = self.maxdTau
				*/
			}

			// check whether we have to switch to energy control arc length
			if(!mIsNrgControlled && mIsConverged)
			{
				if(mDissipatedEnergy > 1000.0)
				{
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					std::cout << "   SWITCH TO DISSIPATED ENERGY ARC LENGTH\n";
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					std::cout << "****************************************************************\n";
					mIsNrgControlled = true;
					dTau = 0.25*1000.0;
					mdTau = dTau;
				}
			}

			// reset dlambda in energy control
			if(mIsNrgControlled)
			{
				mDLambda = 0.0;
			}
		}

		/* 
		===============================================================================
		Perform the finalization of this solution step
		based on the STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL
		===============================================================================
		*/

		model_part.GetProcessInfo()[STRATEGY_SOLUTION_STEP_SOLVED] = 1;

		if(mIsConverged)
		{
			if(verbose) PromptSeparator();
			this->FinalizeSolutionStep(mA, mDx, resi);
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
			mLambda  = lambda;
			mDLambda = Dlambda;
			mdTau    = dTau;
			mCumulatedFactor = 1.0;
		}
		else
		{
			mCumulatedFactor = current_factor;
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
		ModelPart & model_part = BaseType::GetModelPart();

        if(!mpScheme->SchemeIsInitialized())
			mpScheme->Initialize(model_part);

        if(!mpScheme->ElementsAreInitialized())
			mpScheme->InitializeElements(model_part);

        if(!mpConvergenceCriteria->mConvergenceCriteriaIsInitialized)
			mpConvergenceCriteria->Initialize(model_part);

		BaseType::mStiffnessMatrixIsBuilt = false;

        mInitializeWasPerformed = true;

		mIsConverged = false;
    }

    void InitializeSolutionStep()
    {
		ModelPart & model_part = BaseType::GetModelPart();

        mpBuilderAndSolver->ResizeAndInitializeVectors(mpScheme, 
					mpA, mpDx, mpb, 
					model_part);

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        mpBuilderAndSolver->InitializeSolutionStep(model_part, mA, mDx, mb);
        mpScheme->InitializeSolutionStep(model_part, mA, mDx, mb);

		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
		mpConvergenceCriteria->InitializeSolutionStep(model_part, rDofSet, mA, mDx, mb);
		
        mSolutionStepIsInitialized = true;

		mIsConverged = false;
    }

	void FinalizeSolutionStep(TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
	{
		ModelPart & model_part = BaseType::GetModelPart();

		if (mCalculateReactionsFlag) {
            mpBuilderAndSolver->CalculateReactions(mpScheme, model_part, mA, mDx, mb);
		}

        mpScheme->FinalizeSolutionStep(model_part, mA, mDx, mb);
		mpScheme->Clean();
		
        mpBuilderAndSolver->FinalizeSolutionStep(model_part, mA, mDx, mb);

		if(mReformDofSetAtEachStep) 
			this->Clear();

		mSolutionStepIsInitialized = false;
	}
	
	void AbortSolutionStep()
	{
		ModelPart & model_part = BaseType::GetModelPart();

		if(BaseType::MoveMeshFlag())
		{
			for (ModelPart::NodeIterator i = model_part.NodesBegin(); i != model_part.NodesEnd(); ++i)
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
		ModelPart & model_part = BaseType::GetModelPart();
		std::vector<double> error_codes;
		for(ModelPart::ElementIterator elem_iter = model_part.ElementsBegin(); elem_iter != model_part.ElementsEnd(); elem_iter++)
		{
			Element& elem = *elem_iter;
			elem.GetValueOnIntegrationPoints(CONSTITUTIVE_INTEGRATION_ERROR_CODE, error_codes, model_part.GetProcessInfo());
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

		ModelPart & model_part = BaseType::GetModelPart();

        mpBuilderAndSolver->Check(model_part);
        mpScheme->Check(model_part);
        mpConvergenceCriteria->Check(model_part);

        KRATOS_CATCH("")
		
		return 0;
    }

protected:

    ArcLengthDissNrgStrategy(const ArcLengthDissNrgStrategy& Other)
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
	double mDLambda;
	double mdTau;
	double mFactor;
	double mCumulatedFactor;
	double mTotalFactor;
	bool mIsNrgControlled;
	double mDissipatedEnergy;

    bool            mReformDofSetAtEachStep;
	bool            mKeepSystemConstantDuringIterations;
	unsigned int    mMaxIterationNumber;
    bool            mCalculateReactionsFlag;
    bool            mInitializeWasPerformed;
	bool            mSolutionStepIsInitialized;
	bool            mIsConverged;
	
	bool            mAllowNonConvergence;

}; /* Class ArcLengthDissNrgStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_ARCLENGTH_DISS_NRG_STRATEGY  defined */
