//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_RESIDUALBASED_LINEAR_STRATEGY )
#define  KRATOS_RESIDUALBASED_LINEAR_STRATEGY

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Enum's
///@{

///@}
///@name Functions
///@{

///@}
///@name Kratos Classes 
///@{

/**   
 * @class ResidualBasedLinearStrategy
 * @ingroup KratosCore
 * @brief This is a very simple strategy to solve linearly the problem
 * @details As a linear strategy the check on the convergence is not done and just one non linear iteration will be performed
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedLinearStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions */
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedLinearStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;
    
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default constructor 
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param CalculateReactionFlag The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param CalculateNormDxFlag The flag sets if the norm of Dx is computed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    ResidualBasedLinearStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool CalculateReactionFlag = false,
        bool ReformDofSetAtEachStep = false,
        bool CalculateNormDxFlag = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY

        mCalculateReactionsFlag = CalculateReactionFlag;
        mReformDofSetAtEachStep = ReformDofSetAtEachStep;
        mCalculateNormDxFlag = CalculateNormDxFlag;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                 new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (mpLinearSolver)
                             );

        // Set flag to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        this->SetEchoLevel(1);

        // By default the matrices are rebuilt at each solution step
        BaseType::SetRebuildLevel(1);

        KRATOS_CATCH("")
    }

    /**
     * Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param CalculateReactionFlag The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param CalculateNormDxFlag The flag sets if the norm of Dx is computed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    ResidualBasedLinearStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        bool CalculateReactionFlag = false,
        bool ReformDofSetAtEachStep = false,
        bool CalculateNormDxFlag = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY

        mCalculateReactionsFlag = CalculateReactionFlag;
        mReformDofSetAtEachStep = ReformDofSetAtEachStep;
        mCalculateNormDxFlag = CalculateNormDxFlag;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the  builder and solver
        mpBuilderAndSolver = pNewBuilderAndSolver;

        // Set flag to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        //set EchoLevel to the default value (only time is displayed)
        this->SetEchoLevel(1);

        // By default the matrices are rebuilt at each solution step
        BaseType::SetRebuildLevel(1);

        KRATOS_CATCH("")
    }

    /** 
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~ResidualBasedLinearStrategy() override
    {
      this->Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    /**
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);
    }

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }
    
    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY
        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        //if(mInitializeWasPerformed == false)
        //{
        //	Initialize();
        //	mInitializeWasPerformed = true;
        //}

        ////initialize solution step
        //if (mSolutionStepIsInitialized == false)
        //	InitializeSolutionStep();


        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        this->GetScheme()->Predict(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        KRATOS_CATCH("")
    }

    /**
     * @brief The problem of interest is solved
     * @details a double containing norm(Dx) is returned if CalculateNormDxFlag == true, else 0 is returned
     * @return norm(Dx)
     */

    double Solve() override
    {
        KRATOS_TRY

        //pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        ProcessInfo& pCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
        {
            Initialize();
            mInitializeWasPerformed = true;
        }

        //prints informations about the current time
        if (BaseType::GetEchoLevel() != 0 && rank == 0)
        {
            KRATOS_INFO("CurrentTime") << "\nCurrentTime = " << pCurrentProcessInfo[TIME] << std::endl;
        }

        //initialize solution step
        if (mSolutionStepIsInitialized == false)
        {
            InitializeSolutionStep();
            mSolutionStepIsInitialized = true;
        }

        //updates the database with a prediction of the solution
        Predict();

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

	pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);


        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);
            // passing smart pointers instead of references here
            // to prevent dangling pointer to system matrix when
            // reusing ml preconditioners in the trilinos tpl
            pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
            BaseType::mStiffnessMatrixIsBuilt = true;
        }
        else
        {
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);
            pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }

        if (BaseType::GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << mA << std::endl;
            KRATOS_INFO("Dx") << "Solution obtained = " << mDx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << mb << std::endl;
        }
        if (this->GetEchoLevel() == 4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] <<  ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), mA, false);
            
            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), mb);
        }

        //update results
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();
        pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        //move the mesh if needed
        if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

        //calculate if needed the norm of Dx
        double normDx = 0.00;
        if (mCalculateNormDxFlag == true)
        {
            normDx = TSparseSpace::TwoNorm(mDx);
        }

        // Calculate reactions if required
        if (mCalculateReactionsFlag == true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        
        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation
        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        
        //deallocate the systemvectors if needed
        if (mReformDofSetAtEachStep == true)
        {
            if (rank == 0 && BaseType::GetEchoLevel() > 0) std::cout << "Clearing System" << std::endl;
            this->Clear();
        }

        //Cleaning memory after the solution
        pScheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        return normDx;

        KRATOS_CATCH("")

    }

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
    }
    
    /**
     * @brief This method directly sets the input as the LHS
     * @param A The LHS matrix
     */
    void GetDirectSystemMatrix(TSystemMatrixType& A)
    {
        A = *mpA;
    }

    /**
     * @brief This method returns the residual norm
     * @return The residual norm
     */
    double GetResidualNorm() override
    {
        if (TSparseSpace::Size(*mpb) != 0)
            return TSparseSpace::TwoNorm(*mpb);
        else
            return 0.0;
    }

    /**
     * @brief This operations should be called before printing the results when non trivial results (e.g. stresses)
    need to be calculated given the solution of the step
     *@details This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
    void CalculateOutputData() override
    {
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
        GetScheme()->CalculateOutputData(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    typename TLinearSolver::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employed

    TSystemVectorPointerType mpDx; /// The incremement in the solution
    TSystemVectorPointerType mpb; /// The RHS vector of the system of equations
    TSystemMatrixPointerType mpA; /// The LHS matrix of the system of equations

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    bool mCalculateNormDxFlag; /// Calculates if required the norm of the correction term Dx

    /**
     * @brief Flag telling if it is needed or not to compute the reactions
     * @details default = true
     */
    bool mCalculateReactionsFlag;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    ///@}
    ///@name Private Operators*/
    ///@{

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        
        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            KRATOS_INFO("Entering Initialize") << "Entering in the Initialize of the ResidualBasedLinearStrategy" << std::endl;

        //pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();

        //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
        if (pScheme->SchemeIsInitialized() == false)
            pScheme->Initialize(BaseType::GetModelPart());

        //Initialize The Elements - OPERATIONS TO BE DONE ONCE
        if (pScheme->ElementsAreInitialized() == false)
            pScheme->InitializeElements(BaseType::GetModelPart());

        //Initialize The Conditions - OPERATIONS TO BE DONE ONCE                                                          
        if (pScheme->ConditionsAreInitialized() == false)
            pScheme->InitializeConditions(BaseType::GetModelPart());

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            KRATOS_INFO("Exiting Initialize") << "Exiting the  Initialize of the ResidualBasedLinearStrategy" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     * @todo Boost dependencies should be replaced by std equivalent
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = GetScheme();

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            KRATOS_INFO("Entering InitializeSolutionStep") << "Entering in the InitializeSolutionStep of the ResidualBasedLinearStrategy" << std::endl;


        // Loop to reform the dofset
        boost::timer system_construction_time;
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
        {
            boost::timer setup_dofs_time;
            //setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());
            if (BaseType::GetEchoLevel() > 0 && rank == 0)
                KRATOS_INFO("setup_dofs_time") << "setup_dofs_time : " << setup_dofs_time.elapsed() << std::endl;

            //shaping correctly the system
            boost::timer setup_system_time;
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
            if (BaseType::GetEchoLevel() > 0 && rank == 0)
                KRATOS_INFO("setup_system_time") << "setup_system_time : " << setup_system_time.elapsed() << std::endl;

            //setting up the Vectors involved to the correct size
            boost::timer system_matrix_resize_time;
            pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA, mpDx, mpb, BaseType::GetModelPart().Elements(), BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
            if (BaseType::GetEchoLevel() > 0 && rank == 0)
                KRATOS_INFO("system_matrix_resize_time") << "system_matrix_resize_time : " << system_matrix_resize_time.elapsed() << std::endl;
        }
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
            KRATOS_INFO("Construction time") << "System Construction Time : " << system_construction_time.elapsed() << std::endl;


        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        //initial operations ... things that are constant over the Solution Step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        //initial operations ... things that are constant over the Solution Step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            KRATOS_INFO("Exiting InitializeSolutionStep") << "Exiting the InitializeSolutionStep of the ResidualBasedLinearStrategy" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method prints information after reach the max number of interations
     */

    void MaxIterationsExceeded()
    {
        std::cout << "***************************************************" << std::endl;
        std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
        std::cout << "***************************************************" << std::endl;
    }

    /**
     * @brief Clears the internal storage
     * @note NULL could be changed to nullptr in the future (c++11)
     */
    void Clear() override
    {
        KRATOS_TRY;
        // If the preconditioner is saved between solves, it
        // should be cleared here.
        GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpA != NULL)
            SparseSpaceType::Clear(mpA);
        if (mpDx != NULL)
            SparseSpaceType::Clear(mpDx);
        if (mpb != NULL)
            SparseSpaceType::Clear(mpb);

        // Setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        KRATOS_CATCH("");
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Private Operations*/
    ///@{


    ///@}
    ///@name Private  Access */
    ///@{


    ///@}
    ///@name Private Inquiry */
    ///@{


    ///@}
    ///@name Un accessible methods */
    ///@{

    /** Copy constructor.
     */
    ResidualBasedLinearStrategy(const ResidualBasedLinearStrategy& Other);

    ///@}

}; /* Class ResidualBasedLinearStrategy */

///@}

///@name Type Definitions */
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_LINEAR_STRATEGY  defined */

