// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main author:    Manuel Messmer
//

#if !defined(KRATOS_PREBUCKLING_STRATEGY)
#define KRATOS_PREBUCKLING_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"

// Application includes
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class PrebucklingStrategy
 * @ingroup StructuralMechanicsApplication
 * @brief Strategy for linearized prebuckling analysis
 * @author Manuel Messmer
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver>
class PrebucklingStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PrebucklingStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType::Pointer SchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> ConvergenceCriteriaType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pEigenSolver The generalized eigenvalue problem solver employed
     * @param pBuilderAndSolver The builder and solver employed
     * @param pConvergenceCriteria The convergence criteria employed
     * @param MaxIteration The maximum number of non-linear iterations
     * @param InitialLoadIncrement Load increment of the first load step
     * @param SmallLoadIncrement Load increment of the small load step
     * @param PathFollowingStep Load increment of the big load step
     * @param ConvergenceRatio Convergence ratio for the computed eigenvalues
     */
    PrebucklingStrategy(
        ModelPart &rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pEigenSolver,
        BuilderAndSolverPointerType pBuilderAndSolver,
        typename ConvergenceCriteriaType::Pointer pConvergenceCriteria,
        int MaxIteration,
        double InitialLoadIncrement,
        double SmallLoadIncrement,
        double PathFollowingStep,
        double ConvergenceRatio )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpEigenSolver = pEigenSolver;

        mpBuilderAndSolver = pBuilderAndSolver;

        mpConvergenceCriteria = pConvergenceCriteria;

        mMaxIteration = MaxIteration;

        mInitialLoadIncrement = InitialLoadIncrement;

        mSmallLoadIncrement = SmallLoadIncrement;

        mPathFollowingStep = PathFollowingStep;

        mConvergenceRatio = ConvergenceRatio;

        // Set Eigensolver flags
        mpEigenSolver->SetDofSetIsInitializedFlag(false);
        mpEigenSolver->SetReshapeMatrixFlag(false);
        mpEigenSolver->SetCalculateReactionsFlag(false);
        // Set Static Builder and Solver flags
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        mpBuilderAndSolver->SetCalculateReactionsFlag(false);
        mpBuilderAndSolver->SetReshapeMatrixFlag(false);

        // Set EchoLevel to the default value (only time is displayed)
        this->SetEchoLevel(1);

        // Default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);

        // Set Matrices and Vectors to empty pointers
        mpStiffnessMatrix = TSparseSpace::CreateEmptyMatrixPointer();
        mpStiffnessMatrixPrevious = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpRHS = TSparseSpace::CreateEmptyVectorPointer();

        rModelPart.GetProcessInfo()[TIME] = 1.0;

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    PrebucklingStrategy(const PrebucklingStrategy &Other) = delete;

    /**
     * @brief Destructor.
     */
    ~PrebucklingStrategy() override
    {
        this->Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Get method for the scheme
     * @return mpScheme: The pointer to the scheme considered
     */
    SchemePointerType &pGetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Get method for the generalized eigenvalue problme solver
     * @return mpEigenSolver: The pointer to the eigen solver considered
     */
    BuilderAndSolverPointerType &pGetEigenSolver()
    {
        return mpEigenSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    BuilderAndSolverPointerType &pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    /**
     * @brief Get method for the convergence criteria
     * @return mpConvergenceCriteria: The pointer to the convergence criteria considered
     */
    ConvergenceCriteriaType &GetConvergenceCriteria()
    {
        return mpConvergenceCriteria;
    }

    /**
     * @brief Get method for solution found flag
     * @return mSolutionFound: Flag that indicates if solution is found
     */
    bool GetSolutionFoundFlag()
    {
        return mSolutionFound;
    }

    /**
     * @brief This sets the level of echo for the solving strategy
     * @param Level of echo for the solving strategy
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic informations
     * 2 -> Printing linear solver data
     * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
     * }
     */
    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        this->pGetEigenSolver()->SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /**
     * @brief Initialization of member variables and prior operation
     */
    void Initialize() override
    {
        KRATOS_TRY

        ModelPart &rModelPart = BaseType::GetModelPart();

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
            << "Entering Initialize" << std::endl;

        if (!mInitializeWasPerformed)
        {
            SchemePointerType &pScheme = this->pGetScheme();

            if ( !pScheme->SchemeIsInitialized() )
                pScheme->Initialize(rModelPart);

            if ( !pScheme->ElementsAreInitialized() )
                pScheme->InitializeElements(rModelPart);

            if ( !pScheme->ConditionsAreInitialized() )
                pScheme->InitializeConditions(rModelPart);
        }
        // Initialization of the convergence criteria
        mpConvergenceCriteria->Initialize(BaseType::GetModelPart());

        mInitializeWasPerformed = true;

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
            << "Exiting Initialize" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY

        BuilderAndSolverPointerType &pBuilderAndSolver = this->pGetBuilderAndSolver();
        pBuilderAndSolver->GetLinearSystemSolver()->Clear();
        BuilderAndSolverPointerType &pEigenSolver = this->pGetEigenSolver();
        pEigenSolver->GetLinearSystemSolver()->Clear();

        if (mpStiffnessMatrix != nullptr)
            mpStiffnessMatrix = nullptr;

        if (mpStiffnessMatrixPrevious != nullptr)
            mpStiffnessMatrixPrevious = nullptr;

        if (mpRHS != nullptr)
            mpRHS = nullptr;

        if (mpDx != nullptr)
            mpDx = nullptr;

        // Re-setting internal flag to ensure that the dof sets are recalculated
        pBuilderAndSolver->SetDofSetIsInitializedFlag(false);
        pEigenSolver->SetDofSetIsInitializedFlag(false);

        pBuilderAndSolver->Clear();
        pEigenSolver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;
        mSolutionStepIsInitialized = false;


        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY
        if (!mSolutionStepIsInitialized){
            ModelPart &rModelPart = BaseType::GetModelPart();

            KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
                << "Entering InitializeSolutionStep" << std::endl;

            // Solver
            BuilderAndSolverPointerType &pBuilderAndSolver = this->pGetBuilderAndSolver();
            // Scheme
            SchemePointerType &pScheme = this->pGetScheme();
            // Convergence Criteria
            typename ConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
            // Member Matrices and Vectors
            SparseMatrixType& rStiffnessMatrix  = *mpStiffnessMatrix;
            SparseVectorType& rRHS  = *mpRHS;
            SparseVectorType& rDx  = *mpDx;

            // Initialize dummy vectors
            SparseVectorPointerType _pDx = SparseSpaceType::CreateEmptyVectorPointer();
            SparseVectorPointerType _pb = SparseSpaceType::CreateEmptyVectorPointer();

            // Reset solution dofs
            BuiltinTimer system_construction_time;
            if ( !pBuilderAndSolver->GetDofSetIsInitializedFlag() ||
                pBuilderAndSolver->GetReshapeMatrixFlag() )
            {
                // Set up list of dofs
                BuiltinTimer setup_dofs_time;

                pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);

                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 )
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                // Set global equation ids
                BuiltinTimer setup_system_time;

                pBuilderAndSolver->SetUpSystem(rModelPart);

                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 )
                    << setup_system_time.ElapsedSeconds() << std::endl;

                // Resize and initialize system matrices
                BuiltinTimer system_matrix_resize_time;

                // Elastic Stiffness matrix; Solution Vector; RHS Vector
                pBuilderAndSolver->ResizeAndInitializeVectors(
                    pScheme, mpStiffnessMatrix, mpDx, mpRHS, rModelPart);
                // Previous Stiffness Matrix
                pBuilderAndSolver->ResizeAndInitializeVectors(
                    pScheme, mpStiffnessMatrixPrevious, _pDx, _pb, rModelPart);

                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 )
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;

            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 )
                << system_construction_time.ElapsedSeconds() << std::endl;

            // Initial operations ... things that are constant over the solution
            // step
            pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),
                                                rStiffnessMatrix, rDx, rRHS);

            // Initial operations ... things that are constant over the solution
            // step
            pScheme->InitializeSolutionStep(BaseType::GetModelPart(), rStiffnessMatrix, rDx, rRHS);

            pConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

            mSolutionStepIsInitialized = true;

            KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
                << "Exiting InitializeSolutionStep" << std::endl;
        }
        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     * @param mLambda Smallest eigenvalue of current solution step
     * @param mLambdaPrev Smallest eigenvalue of previous solution step
     * @param delta_load_multiplier Load increment between two load steps
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& rModelPart = BaseType::GetModelPart();
        SchemePointerType& pScheme = this->pGetScheme();
        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        SparseMatrixType& rStiffnessMatrix  = *mpStiffnessMatrix;
        SparseMatrixType& rStiffnessMatrixPrevious = *mpStiffnessMatrixPrevious;
        SparseVectorType& rRHS  = *mpRHS;
        SparseVectorType& rDx  = *mpDx;

        typename ConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

        // Initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;

        // Update loadfactor increment
        double delta_load_multiplier = 0.0;
        if( mLoadStepIteration == 1) // inital load increment
        {
            delta_load_multiplier = mInitialLoadIncrement*(mLambda + mLambdaPrev);
        }
        else if( mLoadStepIteration % 2 == 1 ) //small load increment
        {
            delta_load_multiplier = mSmallLoadIncrement*(mLambdaPrev );
        }

        BuiltinTimer system_solve_time;
        // Initialize nonlinear iteration
        this->pGetScheme()->InitializeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
        pConvergenceCriteria->InitializeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
        is_converged = mpConvergenceCriteria->PreCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

        TSparseSpace::SetToZero(rStiffnessMatrix);
        TSparseSpace::SetToZero(rRHS);
        TSparseSpace::SetToZero(rDx);
        // Build system and solve system
        pBuilderAndSolver->BuildAndSolve(pScheme, rModelPart, rStiffnessMatrix, rDx, rRHS);
        // Update internal variables
        pScheme->Update(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS );
        BaseType::MoveMesh();
        // Finalize nonlinear solution step
        pScheme->FinalizeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
        pConvergenceCriteria->FinalizeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

        if (is_converged){
            is_converged = mpConvergenceCriteria->PostCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
        }

        // Start iteration cycle
        while ( !is_converged &&
            iteration_number++ < mMaxIteration)
        {
            // Setting the number of iteration
            rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            // Initialize nonlinear iteration step
            pScheme->InitializeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
            pConvergenceCriteria->InitializeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
            is_converged = mpConvergenceCriteria->PreCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

            TSparseSpace::SetToZero(rStiffnessMatrix);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rRHS);
            // Build and solve system
            pBuilderAndSolver->BuildAndSolve(pScheme, rModelPart, rStiffnessMatrix,rDx, rRHS);

            // Update internal variables
            this->pGetScheme()->Update(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
            BaseType::MoveMesh();

            // Finalize nonlinear iteration step
            this->pGetScheme()->FinalizeNonLinIteration( rModelPart,rStiffnessMatrix, rDx, rRHS );
            pConvergenceCriteria->FinalizeNonLinearIteration(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);

            if (is_converged){
                is_converged = mpConvergenceCriteria->PostCriteria(rModelPart, pBuilderAndSolver->GetDofSet(), rStiffnessMatrix, rDx, rRHS);
            }
        }
        KRATOS_INFO_IF("Nonlinear Loadstep Time: ", BaseType::GetEchoLevel() > 0)
                << system_solve_time.ElapsedSeconds() << std::endl;

        if ( !is_converged ) {
            KRATOS_INFO_IF("Nonlinear Loadstep: ", this->GetEchoLevel() > 0)
                << "Convergence not achieved after ( " << mMaxIteration
                << " ) Iterations !" << std::endl;
        } else {
            KRATOS_INFO_IF("Nonlinear Loadstep: ", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << mMaxIteration << " iterations" << std::endl;
        }

        // Vector and matrix are initialized by eigensolver
        DenseVectorType Eigenvalues;
        DenseMatrixType Eigenvectors;

        if( mLoadStepIteration % 2 == 0 ){
            //Copy matrices after path following step and initial step
            rStiffnessMatrixPrevious = rStiffnessMatrix;

            if( mLoadStepIteration > 0){
                // Update mLambdaPrev
                mLambdaPrev = mLambdaPrev + mPathFollowingStep*mLambda;
            }
        }
        else if( mLoadStepIteration % 2 == 1 ){
            // Evaluate eigenvalue problem after small load increment
            //#####################################################################################
            //Method 1
            //Poy, Denise Lori-eng; On the Buckling Finite Element Analysis of Beam Structures
            //#####################################################################################
            // Solve Eigenvalue Problem
            // BuiltinTimer eigen_solve_time;
            // this->pGetEigenSolver()->GetLinearSystemSolver()->Solve(
            //     rStiffnessMatrix,
            //     rStiffnessMatrixPrevious,
            //     Eigenvalues,
            //     Eigenvectors);

            // mLambda = 1.0 / ( 1.0 - Eigenvalues(0) ) * delta_load_multiplier;
            // for(int i = 0; i < Eigenvalues.size(); i++ )
            // {
            //     Eigenvalues[i] = mLambdaPrev + 1.0 / ( 1.0 - Eigenvalues[i] )*delta_load_multiplier;
            // }
            // this->AssignVariables(Eigenvalues, Eigenvectors);
            //End Method 1#########################################################################

            //#####################################################################################
            //Method 2:
            // Jia, X.Mang, H. A.; Assessment of solutions from the consistently linearized eigenproblem
            // by means of finite difference approximations
            //#####################################################################################
            // Solve Eigenvalue Problem

            // This implementation of the eigenvalue problem relies on the difference of the stiffness matrix
            // between the current and previous step
            rStiffnessMatrix = rStiffnessMatrixPrevious - rStiffnessMatrix;

            this->pGetEigenSolver()->GetLinearSystemSolver()->Solve(
                rStiffnessMatrixPrevious,
                rStiffnessMatrix,
                Eigenvalues,
                Eigenvectors);

            // Update eigenvalues to loadfactors (Instead of dividing matrix by delta_load_multiplier, here eigenvalues are multiplied)
            mLambda = Eigenvalues(0)*delta_load_multiplier;
            for(unsigned int i = 0; i < Eigenvalues.size(); i++ )
            {
                Eigenvalues[i] = mLambdaPrev + Eigenvalues[i]*delta_load_multiplier;
            }

            this->AssignVariables(Eigenvalues, Eigenvectors);

            // Reset elastic stiffness matrix for next loadstep
            mpStiffnessMatrix = nullptr;
            pBuilderAndSolver->ResizeAndInitializeVectors(
                pScheme, mpStiffnessMatrix, mpDx, mpRHS, rModelPart);
            // End Method 2#########################################################################

            // Convergence criteria for buckling analysis
            if( std::abs(mLambda/mLambdaPrev) < mConvergenceRatio ){
                mSolutionFound = true;
                KRATOS_INFO_IF("Prebuckling Analysis: ", BaseType::GetEchoLevel() > 0)
                << "Convergence achieved in " << mLoadStepIteration + 1 << " Load Iterations!" << std::endl;
            }
        }

        mLoadStepIteration++;
        // Prepare load conditions for next solution step
        UpdateLoadConditions();

        return true;
        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
            << "Entering FinalizeSolutionStep" << std::endl;

        typename ConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
        // Member Matrices and Vectors
        SparseMatrixType& rStiffnessMatrix  = *mpStiffnessMatrix;
        SparseVectorType& rRHS  = *mpRHS;
        SparseVectorType& rDx  = *mpDx;

        pGetBuilderAndSolver()->FinalizeSolutionStep(
            BaseType::GetModelPart(), rStiffnessMatrix, rDx, rRHS);

        pGetScheme()->FinalizeSolutionStep(BaseType::GetModelPart(),
                                           rStiffnessMatrix, rDx, rRHS);

        pConvergenceCriteria->FinalizeSolutionStep(BaseType::GetModelPart(), pGetBuilderAndSolver()->GetDofSet(), rStiffnessMatrix, rDx, rRHS );

        // Cleaning memory after the solution
        pGetScheme()->Clean();
        // Reset flags for next step
        mSolutionStepIsInitialized = false;

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
            << "Exiting FinalizeSolutionStep" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        ModelPart &rModelPart = BaseType::GetModelPart();

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
            << "Entering Check" << std::endl;

        // Check the model part
        BaseType::Check();

        // Check the scheme
        this->pGetScheme()->Check(rModelPart);

        // Check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        KRATOS_INFO_IF("PrebucklingStrategy", BaseType::GetEchoLevel() > 2 )
            << "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

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

    SchemePointerType mpScheme;

    BuilderAndSolverPointerType mpEigenSolver;

    BuilderAndSolverPointerType mpBuilderAndSolver;

    // SparseMatrixPointerType mpMassMatrix;

    SparseMatrixPointerType mpStiffnessMatrix;

    SparseMatrixPointerType mpStiffnessMatrixPrevious;

    SparseVectorPointerType mpRHS;

    SparseVectorPointerType mpDx;

    typename ConvergenceCriteriaType::Pointer mpConvergenceCriteria;
    /// The pointer to the convergence criteria employed

    bool mInitializeWasPerformed = false;

    bool mSolutionStepIsInitialized = false;

    bool mSolutionFound = false;

    unsigned int mLoadStepIteration = 0;

    unsigned int mMaxIteration;

    std::vector<array_1d<double,3>> mpInitialLoads;

    double mInitialLoadIncrement;
    double mSmallLoadIncrement;
    double mPathFollowingStep;
    double mConvergenceRatio;

    double mLambda = 0.0;
    double mLambdaPrev = 1.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Updates load conditions
     * @param TIME is used as load multiplier of external loads
     */
    void UpdateLoadConditions()
    {
        ModelPart &rModelPart = BaseType::GetModelPart();

        if( mLoadStepIteration == 1){
            // Initial load increment
            rModelPart.GetProcessInfo()[TIME] = ( 1.0 + mInitialLoadIncrement );
        }
        else if( mLoadStepIteration % 2 == 0){
            // Do path following step
            rModelPart.GetProcessInfo()[TIME] = ( mLambdaPrev + mPathFollowingStep * mLambda );
        }
        else{
            // Do small load increment
            rModelPart.GetProcessInfo()[TIME] = (1 + mSmallLoadIncrement) * mLambdaPrev;
        }
    }

    /**
     * @brief Assign eigenvalues and eigenvectors to kratos variables
     * @detail Copied from eigensolver_strategy
     */
    void AssignVariables(DenseVectorType &rEigenvalues, DenseMatrixType &rEigenvectors )
    {
        ModelPart &rModelPart = BaseType::GetModelPart();
        const std::size_t NumEigenvalues = rEigenvalues.size();

        // Store eigenvalues in process info
        rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR] = rEigenvalues;

        for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
        {
            ModelPart::NodeType::DofsContainerType &NodeDofs = itNode->GetDofs();
            const std::size_t NumNodeDofs = NodeDofs.size();

            Matrix &rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
            if (rNodeEigenvectors.size1() != NumEigenvalues || rNodeEigenvectors.size2() != NumNodeDofs){
                rNodeEigenvectors.resize(NumEigenvalues, NumNodeDofs, false);
            }

            // TO BE VERIFIED!! In the current implmentation of Dofs there are nor reordered and only pushec back.
            // // the jth column index of EIGENVECTOR_MATRIX corresponds to the jth nodal dof. therefore,
            // // the dof ordering must not change.
            // if (NodeDofs.IsSorted() == false)
            // {
            //     NodeDofs.Sort();
            // }

            // Fill the EIGENVECTOR_MATRIX
            for (std::size_t i = 0; i < NumEigenvalues; i++){
                for (std::size_t j = 0; j < NumNodeDofs; j++){
                    auto itDof = std::begin(NodeDofs) + j;
                    if( !(*itDof)->IsFixed() ){
                        rNodeEigenvectors(i, j) = rEigenvectors(i, (*itDof)->EquationId());
                    }
                    else{
                        rNodeEigenvectors(i, j) = 0.0;
                    }
                }
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; /* Class PrebucklingStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos */

#endif /* KRATOS_EIGENSOLVER_STRATEGY  defined */
