//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $

#if !defined(KRATOS_STATIC_GENERAL_STRATEGY_KRYLOV_NEWTON )
#define  KRATOS_STATIC_GENERAL_STRATEGY_KRYLOV_NEWTON

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

// accelerator utility
#include "custom_utilities/krylov_subspace_accelerator_lapack_utilities.h"

namespace Kratos
{


template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class LineSearchOperator
{
public:

    typedef typename TSparseSpace::DataType           TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;
    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> TBuilderAndSolverType;
    typedef typename TBuilderAndSolverType::Pointer   TBuilderAndSolverPointerType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Scheme<TSparseSpace, TDenseSpace> TSchemeType;
    typedef typename TSchemeType::Pointer  TSchemePointerType;

    typedef ModelPart::DofType TDofType;
    typedef ModelPart::DofsArrayType DofsArrayType;

    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    LineSearchOperator(){};

    ~LineSearchOperator (){};

   double ExecuteLineSearch( TBuilderAndSolverPointerType pBuilderAndSolver,
			     TSchemePointerType pScheme,
			     ModelPart& r_model_part, 
			     TSystemMatrixType& A,
			     TSystemVectorType& Dx,
			     TSystemVectorType& b,
			     double rCurrentAlpha,
			     double rPreviousAlpha )
    {

      KRATOS_TRY

   	//bool linesearch_success=false;


   	//Save Current Displacement, Velocity, Acceleration) as history
   	DofsArrayType HistoricDofSet = pBuilderAndSolver->GetDofSet();

   	DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

   	TSystemVectorType  ReferenceDx = Dx;

   	//Timer::Start("LineSearch");
	     
   	//s0  (alpha=0)
   	double R0= inner_prod(Dx,b);

   	//s1  (alpha=1)
   	Dx = ReferenceDx *1.0;
	
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);

   	     //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	     pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);

   	double R1= inner_prod(ReferenceDx,b);
   	// ** Restore Current Displacement, Velocity, Acceleration
   	Dx *= (-1.0);
   	pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	b.clear();

	//std::cout<<" Initial Slope "<<R0<<" FinalSlope "<<R1<<std::endl;

   	if(R0*R1<0){

   	  //std::cout<<" Enters to the Linesearch iteration "<<R0*R1<<" < 0 "<<std::endl;

	  double R2 = R1;
	  if(fabs(R1)<fabs(R0))
	    R2=R0;
	  double R0start = R0;

	  double alpha = 0.0;
	  double nabla = 0.0;
	  double delta = 1.0;

  	  double CurrentAlpha  = 1.0; 
   	  //double PreviousAlpha = 0.0;

	  int iterations=0;
	  int max_iterations = 10;
	  //std::cout<<" [R1: "<<R1<<", R0: "<<R0<<"]"<<std::endl;

	  while(fabs(R2/R0start)>0.3 && iterations<max_iterations && (R1*R0)<0.0 && fabs(R1)>1.0e-7 && fabs(R0)>1.0e-7) {
	  

	    alpha = 0.5*(nabla+delta);


	    //std::cout<<" R2 = "<<fabs(R2)<<" ?> "<<0.5*R0start<<" =  0.5 * R0_start; R1 "<<R1<<" / ; R0 "<<R0<<std::endl;
	    //std::cout<<" [ CurrentAlpha: "<<CurrentAlpha<<" PreviousAlpha: "<<PreviousAlpha<<" ] --> Computed Alpha: "<<alpha<<std::endl;


	    //alpha_k-1 = alpha_k
	    //PreviousAlpha = CurrentAlpha;

	    //alpha_k  = alpha_k+1
	    CurrentAlpha  = alpha;

  
	    //compute s(alpha_k+1)
	    Dx = ReferenceDx * CurrentAlpha;
	
	    pScheme->Update(r_model_part,rDofSet,A,Dx,b);

	    //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	    pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);
	     
	    //slope_k = slope_k+1
	    R2 = inner_prod(ReferenceDx,b);
	    // ** Restore Current Displacement, Velocity, Acceleration
	    Dx *= (-1.0);
	    pScheme->Update(r_model_part,rDofSet,A,Dx,b);
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
	  //std::cout<<" CurrentSlope = "<<R2<<" ?> "<<0.8*fabs(R0start)<<"=  0.8*InitialSlope;  PreviousSlope "<<R1<<std::endl;
	   
	   
	  rPreviousAlpha = rCurrentAlpha;
	  rCurrentAlpha  = CurrentAlpha;
    
   	}
   	// else{
   	//   Step.CurrentAlpha  = Step.PreviousAlpha;
   	// }
  

   	if(rCurrentAlpha>1 || rCurrentAlpha<=0)
   	  rCurrentAlpha=1;

   	//Restore Current Displacement, Velocity, Acceleration
   	Dx      = ReferenceDx;
   	//rDofSet = HistoricDofSet;
	

   	//Timer::Stop("LineSearch");

   	KRATOS_CATCH( "" )

	return rCurrentAlpha;
      
    }

};



/// Short class definition.
/**   Detail class definition.
 */
template<class TSparseSpace,
         class TDenseSpace, 
         class TLinearSolver
         >
class StaticGeneralStrategyKrylovNewton
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(StaticGeneralStrategyKrylovNewton);

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

    StaticGeneralStrategyKrylovNewton(ModelPart& model_part,
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

    StaticGeneralStrategyKrylovNewton(ModelPart& model_part,
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

    virtual ~StaticGeneralStrategyKrylovNewton()
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
		/* 
		===============================================================================
		Solution step initialization
		===============================================================================
		*/

		mIsConverged = false;

		ModelPart& model = BaseType::GetModelPart();

		if (!mInitializeWasPerformed) this->Initialize();

		int  finalize_level  = model.GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL];
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

			// create the accelerator
			KrylovSubspaceAcceleratorUtilties::KrylovAcceleratorLapack krylov_acc;
			krylov_acc.BeginSolutionStep(TSparseSpace::Size(mDx));
			
			unsigned int iteration_number = 0;
			model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

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

				/* 
				===============================================================================
				Build system matrices and vectors
				===============================================================================
				*/
				TSparseSpace::SetToZero(mDx);
				TSparseSpace::SetToZero(mb);
				if(iteration_number == 1)
				{
					// build the LHS only once...
					TSparseSpace::SetToZero(mA);
					mpBuilderAndSolver->Build(mpScheme, model, mA, mb);
				}
				else
				{
					// ... or if the krylov method requires it
					if(krylov_acc.BuildSystemMatrix())
						mpBuilderAndSolver->Build(mpScheme, model, mA, mb);
					else
						mpBuilderAndSolver->BuildRHS(mpScheme, model, mb);
				}
				
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
					mpBuilderAndSolver->SystemSolve(mA, mDx, mb);
					// check solution validity
					if(!this->CheckSolutionValidity(mDx))
					{
						mIsConverged = false;
						iteration_number = mMaxIterationNumber;
						model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
						break;
					}
				}
				
				/* 
				===============================================================================
				Solve the Least Square problem for the Krylov accelerator
				===============================================================================
				*/

				if(krylov_acc.LeastSquares(mDx)<0)
				{
					if(verbose) this->PromptIntegrationErrors();
					mIsConverged = false;
					iteration_number = mMaxIterationNumber;
					model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
					break;
				}
				
				/* 
				===============================================================================
				Finalize the current nonlinear iteration
				Note: update with the krylov accelerated vector, not with mDx!
				===============================================================================
				*/
				
				krylov_acc.AccelerateSolution(mDx);
				rDofSet = mpBuilderAndSolver->GetDofSet();
				mpScheme->Update(model, rDofSet, mA, mDx, mb);
				if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
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

    StaticGeneralStrategyKrylovNewton(const StaticGeneralStrategyKrylovNewton& Other)
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

}; /* Class StaticGeneralStrategyKrylovNewton */

} /* namespace Kratos.*/

#endif /* KRATOS_STATIC_GENERAL_STRATEGY_KRYLOV_NEWTON  defined */
