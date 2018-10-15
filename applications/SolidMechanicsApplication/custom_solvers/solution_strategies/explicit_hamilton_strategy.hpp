//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EXPLICIT_HAMILTON_STRATEGY)
#define  KRATOS_EXPLICIT_HAMILTON_STRATEGY

/* System includes */


/* External includes */
//#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"

//default builder and solver
#include "custom_solvers/solution_builders_and_solvers/explicit_hamilton_builder_and_solver.hpp"


namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ExplicitHamiltonStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    /** Counted pointer of ClassName */

    KRATOS_CLASS_POINTER_DEFINITION(ExplicitHamiltonStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::BuilderAndSolverType BuilderAndSolverType;

    typedef TLinearSolver LinearSolverType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::SchemeType SchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::SystemMatrixType SystemMatrixType;

    typedef typename BaseType::SystemVectorType SystemVectorType;

    typedef typename BaseType::SystemMatrixPointerType SystemMatrixPointerType;

    typedef typename BaseType::SystemVectorPointerType SystemVectorPointerType;


    /** Constructors.
     */
    ExplicitHamiltonStrategy(
        ModelPart& model_part,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
    }

    ExplicitHamiltonStrategy(
        ModelPart& model_part,
        typename SchemeType::Pointer pScheme,
        typename LinearSolverType::Pointer pNewLinearSolver,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

        //set flags to default values
        mCalculateReactionsFlag = CalculateReactions;

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

	//saving the scheme
        mpScheme = pScheme;

        //saving the linear solver
        mpLinearSolver = pNewLinearSolver; //Not used in explicit strategies

        //setting up the default builder and solver
        mpBuilderAndSolver = typename BuilderAndSolverType::Pointer
                             (
                                 new ExplicitHamiltonBuilderAndSolver <TSparseSpace, TDenseSpace, TLinearSolver > (mpLinearSolver)
                             );

        //set flags to start correcty the calculations
        mSolutionStepIsInitialized  = false;

        mInitializeWasPerformed     = false;

        //set EchoLevel to the deffault value (only time is displayed)
        SetEchoLevel(1);

	ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

	rCurrentProcessInfo[ALPHA_TRAPEZOIDAL_RULE] = 0.5; //Alpha [0,1]  trapezoidal rule
	rCurrentProcessInfo[POSITION_UPDATE_LABEL] = false;
	rCurrentProcessInfo[ROTATION_UPDATE_LABEL] = false;
	rCurrentProcessInfo[MOMENTUM_UPDATE_LABEL] = false;

        BaseType::SetRebuildLevel(0);   //# if (BaseType::mRebuildLevel > 0) the mass matrix is rebuild

        KRATOS_CATCH( "" )
    }

    /** Destructor.
     */
    virtual ~ExplicitHamiltonStrategy()
    {
    }

    /** Destructor.
     */

    //Set and Get Scheme ... containing Builder, Update and other

    void SetScheme(typename SchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    typename SchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

     //Set and Get the BuilderAndSolver

    void SetBuilderAndSolver(typename BuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    typename BuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    //Ser and Get Flags

    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
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

    }

    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    //level of echo for the solving strategy
    // 0 -> mute... no echo at all
    // 1 -> printing time and basic informations
    // 2 -> printing linear solver data
    // 3 -> Print of debug informations:
    //    Echo of stiffness matrix, Dx, b...

    void SetEchoLevel(int Level)
    {
        BaseType::mEchoLevel = Level;
	mpBuilderAndSolver->SetEchoLevel(Level);
    }

    //**********************************************************************
    //**********************************************************************

    void Initialize()
    {
        KRATOS_TRY

        //pointers needed in the solution
        typename SchemeType::Pointer pScheme                     = GetScheme();
        typename BuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        ModelPart& r_model_part                                   = BaseType::GetModelPart();

        SystemMatrixType mA  = SystemMatrixType();   //dummy initialization. Not used in builder and solver.
	SystemVectorType mDx = SystemVectorType();
	SystemVectorType mb  = SystemVectorType();

        //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
        if (pScheme->SchemeIsInitialized() == false)
            pScheme->Initialize(BaseType::GetModelPart());

        //Initialize The Elements - OPERATIONS TO BE DONE ONCE
        if (pScheme->ElementsAreInitialized() == false)
            pScheme->InitializeElements(BaseType::GetModelPart());

        //Initialize The Conditions- OPERATIONS TO BE DONE ONCE
        if (pScheme->ConditionsAreInitialized() == false)
            pScheme->InitializeConditions(BaseType::GetModelPart());

        pBuilderAndSolver->BuildLHS(pScheme, r_model_part, mA); //calculate and store nodal variables: NODAL_MASS and INERTIA_DYADIC

        mInitializeWasPerformed = true;

        KRATOS_CATCH( "" )
    }

    //**********************************************************************
    //**********************************************************************

    void InitializeSolutionStep()
    {
        KRATOS_TRY

        typename SchemeType::Pointer pScheme                     = GetScheme();
        typename BuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        ModelPart& r_model_part                                   = BaseType::GetModelPart();

        SystemMatrixType mA = SystemMatrixType();   //dummy initialization. Not used in builder and solver.
        SystemVectorType mDx = SystemVectorType();
        SystemVectorType mb = SystemVectorType();

        //initial operations ... things that are constant over the Solution Step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        if(BaseType::mRebuildLevel > 0)
        {
	  pBuilderAndSolver->BuildLHS(pScheme, r_model_part, mA); //calculate and store nodal variables: NODAL_MASS and INERTIA_DYADIC
	}

	// if(BaseType::GetModelPart().GetProcessInfo()[TIME]<=BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]){
	//   std::cout<<" INITIAL CONDITIONS:: initial exernal forces "<<std::endl;
	//   pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb); //fills FORCE_RESIDUAL and MOMENT_RESIDUAL nodal variables
	//   //Position_Momentum and Rotation_Momentum Update (assign nodal contributions):
	//   DofsArrayType rDofSet; //dummy initialization. Not used in builder and solver
	//   BaseType::GetModelPart().GetProcessInfo()[POSITION_UPDATE_LABEL] = false;
	//   BaseType::GetModelPart().GetProcessInfo()[ROTATION_UPDATE_LABEL] = false;
	//   BaseType::GetModelPart().GetProcessInfo()[MOMENTUM_UPDATE_LABEL] = true;
	//   pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb); //SECOND CALL: updates Nodal  Momentum Variables ->POSITION_MOMENTUM,ROTATION_MOMENTUM :: (FORCE_RESIDUAL,MOMENT_RESIDUAL) needed -> (Scheme::mUpdateMomentumFlag = true)
	// }

        mSolutionStepIsInitialized = true;

        KRATOS_CATCH( "" )
    }


    //**********************************************************************
    //**********************************************************************
    /*
                      SOLUTION OF THE PROBLEM OF INTEREST
     */
    //**********************************************************************


    double Solve()
    {
        KRATOS_TRY

        DofsArrayType rDofSet; //dummy initialization. Not used in builder and solver
        SystemMatrixType mA  = SystemMatrixType();
        SystemVectorType mDx = SystemVectorType();
        SystemVectorType mb  = SystemVectorType();

        //pointers needed in the solution
        typename SchemeType::Pointer pScheme = GetScheme();
        typename BuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if(mInitializeWasPerformed == false)
            Initialize();

        //prints informations about the current time
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << " " << std::endl;
            std::cout << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        //initialize solution step
        if(mSolutionStepIsInitialized == false)
	  InitializeSolutionStep();

	//1) Position explicit strategy:

	//1.1) Update nodal displacements:
	bool update_at_start = false;
	if( update_at_start ){
	  BaseType::GetModelPart().GetProcessInfo()[POSITION_UPDATE_LABEL] = true;
	  BaseType::GetModelPart().GetProcessInfo()[ROTATION_UPDATE_LABEL] = false;
	  BaseType::GetModelPart().GetProcessInfo()[MOMENTUM_UPDATE_LABEL] = false;
	  pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb); //FIST CALL: updates Nodal Movement Variables ->DISPLACEMENT,VELOCITY,ACCELERATION (FORCE_RESIDUAL,MOMENT_RESIDUAL,POSITION_MOMENTUM) needed -> (Scheme::mUpdatePositionFlag = true)
	}
	//1.2) Update nodal positions:	(at end)
	//if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

	//2) Rotation explicit strategy (nodally implicit):

	//2.1) Build RHS and Solve nodal implicit scheme:
        pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb); // Implicitly integrates the equation of rotation  ->fills STEP_ROTATION (MOMENT_RESIDUAL,ROTATION_MOMENTUM) needed

	//2.2) Update nodal rotations:
	BaseType::GetModelPart().GetProcessInfo()[POSITION_UPDATE_LABEL] = false;
	BaseType::GetModelPart().GetProcessInfo()[ROTATION_UPDATE_LABEL] = true;
	BaseType::GetModelPart().GetProcessInfo()[MOMENTUM_UPDATE_LABEL] = false;
        pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb); //FIST CALL: updates Nodal Movement Variables ->DISPLACEMENT,VELOCITY,ACCELERATION,ROTATION,ANGULAR_VELOCITY,ANGULAR_ACCELERATION  (FORCE_RESIDUAL,MOMENT_RESIDUAL,POSITION_MOMENTUM) needed -> (Scheme::mUpdateRotationFlag = true)

	if( !update_at_start ){
	  BaseType::GetModelPart().GetProcessInfo()[POSITION_UPDATE_LABEL] = true;
	  BaseType::GetModelPart().GetProcessInfo()[ROTATION_UPDATE_LABEL] = false;
	  BaseType::GetModelPart().GetProcessInfo()[MOMENTUM_UPDATE_LABEL] = false;
	  pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb); //FIST CALL: updates Nodal Movement Variables ->DISPLACEMENT,VELOCITY,ACCELERATION (FORCE_RESIDUAL,MOMENT_RESIDUAL,POSITION_MOMENTUM) needed -> (Scheme::mUpdatePositionFlag = true)
	}

	//calculate reactions if required
	//if (mCalculateReactionsFlag == true)
	//{
	//   pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
	//}

	//ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

	//3) Update momentum equations:

	//3.1) Build RHS (redual forces and moments update):
	pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb); //fills FORCE_RESIDUAL and MOMENT_RESIDUAL nodal variables

	//3.2) Position_Momentum and Rotation_Momentum Update (assign nodal contributions):
	BaseType::GetModelPart().GetProcessInfo()[POSITION_UPDATE_LABEL] = false;
	BaseType::GetModelPart().GetProcessInfo()[ROTATION_UPDATE_LABEL] = false;
	BaseType::GetModelPart().GetProcessInfo()[MOMENTUM_UPDATE_LABEL] = true;
	pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb); //SECOND CALL: updates Nodal  Momentum Variables ->POSITION_MOMENTUM,ROTATION_MOMENTUM :: (FORCE_RESIDUAL,MOMENT_RESIDUAL) needed -> (Scheme::mUpdateMomentumFlag = true)

	if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

 	//Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation
        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb); //FinalizeSolution step of elements and conditions


        //Cleaning memory after the solution
        pScheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        return 0.00;

        KRATOS_CATCH( "" )

    }

    //**********************************************************************
    //**********************************************************************

    void Clear()
    {
        KRATOS_TRY
        std::cout << "Explicit strategy Clear function used" << std::endl;

       //setting to zero the internal flag to ensure that the dof sets are recalculated
        //GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        //GetBuilderAndSolver()->Clear();

        GetScheme()->Clear();

        KRATOS_CATCH( "" )
    }




    /*@} */
    /**@name Operators
     */
    /*@{ */

    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */

    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */



    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

protected:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    typename SchemeType::Pointer mpScheme;

    typename LinearSolverType::Pointer mpLinearSolver;

    typename BuilderAndSolverType::Pointer mpBuilderAndSolver;

    SystemVectorPointerType mpDx;
    SystemVectorPointerType mpb;
    SystemMatrixPointerType mpA;


    /**
    Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    - true  => reforme at each time step
    - false => form just one (more efficient)

    Default = false
     */
    bool mReformDofSetAtEachStep;

    /**
    Flag telling if it is needed or not to compute the reactions

    default = true
     */
    bool mCalculateReactionsFlag;

    bool mSolutionStepIsInitialized;

    bool mInitializeWasPerformed;

    bool mComputeTime;


    /*@} */
    /**@name Private Operators*/
    /*@{ */

      //**********************************************************************
    //**********************************************************************


    void CalculateReactions()
    {

    }

    //**********************************************************************
    //**********************************************************************

    /**
     * function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */

    int Check()
    {
        KRATOS_TRY

        BaseType::Check();

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH( "" )
    }


//***************************************************************************
//***************************************************************************


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /** Copy constructor.
     */
    ExplicitHamiltonStrategy(const ExplicitHamiltonStrategy& Other)
    {
    };


    /*@} */

}; /* Class ExplicitHamiltonStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_HAMILTON_BEAM_STRATEGY  defined */

