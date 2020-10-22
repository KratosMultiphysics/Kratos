//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(MOR_OFFLINE_SECOND_ORDER_STRATEGY)
#define MOR_OFFLINE_SECOND_ORDER_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_complex_interface.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"

#include "mor_application_variables.h"
#include "custom_utilities/dirichlet_utility.hpp"

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
 * @class MorOfflineSecondOrderStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
          class TReducedSparseSpace,
          class TReducedDenseSpace
          >
class MorOfflineSecondOrderStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace>>
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorOfflineSecondOrderStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace>> BaseType;

    // typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;
    typedef BuilderAndSolver< TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace> > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef TLinearSolver TLinearSolverType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    // typedef typename TReducedSparseSpace::MatrixType TReducedSparseMatrixType;

    typedef typename TReducedDenseSpace::MatrixType TReducedDenseMatrixType;

    typedef typename TReducedDenseSpace::VectorType TReducedDenseVectorType;

    typedef typename TReducedSparseSpace::MatrixPointerType TReducedSparseMatrixPointerType;

    typedef typename TReducedSparseSpace::VectorPointerType TReducedSparseVectorPointerType;

    typedef typename TReducedDenseSpace::MatrixPointerType TReducedDenseMatrixPointerType;

    typedef typename TReducedDenseSpace::VectorPointerType TReducedDenseVectorPointerType;



    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor for symmetric problems
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     * @param UseDefinedOutputFlag If true, the solution is only computed at specified dofs
     */
    MorOfflineSecondOrderStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        bool UseDampingFlag = false)
        : SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace>>(rModelPart, false),
            mpScheme(pScheme),
            mpBuilderAndSolver(pBuilderAndSolver),
            mpLinearSolver(pNewLinearSolver),
            mUseDamping(UseDampingFlag),
            mSystemIsSymmetric(true)
    {
        KRATOS_TRY;

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;
        mImportedSystem = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        // GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(0);

        // mpK = TSparseSpace::CreateEmptyMatrixPointer();
        // mpKi = TSparseSpace::CreateEmptyMatrixPointer();
        // mpC = TSparseSpace::CreateEmptyMatrixPointer();
        // mpM = TSparseSpace::CreateEmptyMatrixPointer();
        // mpMi = TSparseSpace::CreateEmptyMatrixPointer();
        // mpRHS = TSparseSpace::CreateEmptyVectorPointer();
        // mpOV = TSparseSpace::CreateEmptyVectorPointer();

        // mpKr = TReducedDenseSpace::CreateEmptyMatrixPointer();
        // mpCr = TReducedDenseSpace::CreateEmptyMatrixPointer();
        // mpMr = TReducedDenseSpace::CreateEmptyMatrixPointer();
        // mpRHSr = TReducedDenseSpace::CreateEmptyVectorPointer();
        // mpOVr = TReducedDenseSpace::CreateEmptyVectorPointer();
        // mpBasis = TReducedDenseSpace::CreateEmptyMatrixPointer();
        // mpBasisLeft = TReducedDenseSpace::CreateEmptyMatrixPointer();

        KRATOS_CATCH("");
    }

    /**
     * Constructor for unsymmetric problems
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     * @param UseDefinedOutputFlag If true, the solution is only computed at specified dofs
     */
    MorOfflineSecondOrderStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        typename TLinearSolverType::Pointer pNewAdjointLinearSolver,
        bool UseDampingFlag = false)
        : SolvingStrategy<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace,TDenseSpace>>(rModelPart, false),
            mpScheme(pScheme),
            mpBuilderAndSolver(pBuilderAndSolver),
            mpLinearSolver(pNewLinearSolver),
            mpAdjointLinearSolver(pNewAdjointLinearSolver),
            mUseDamping(UseDampingFlag),
            mSystemIsSymmetric(false)
    {
        KRATOS_TRY;

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;
        mImportedSystem = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        // GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(0);

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorOfflineSecondOrderStrategy() override
    {
        Clear();
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
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetLinearSolver(typename TLinearSolverType::Pointer pNewLinearSolver)
    {
        mpLinearSolver = pNewLinearSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TLinearSolverType::Pointer GetLinearSolver()
    {
        return mpLinearSolver;
    };

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
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
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Initialization of member variables and prior operations
     */
    virtual void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            std::cout << "Hello I am initializing the offline strategy" << std::endl;

            mpK = TSparseSpace::CreateEmptyMatrixPointer();
            mpKi = TSparseSpace::CreateEmptyMatrixPointer();
            mpC = TSparseSpace::CreateEmptyMatrixPointer();
            mpM = TSparseSpace::CreateEmptyMatrixPointer();
            mpMi = TSparseSpace::CreateEmptyMatrixPointer();
            mpRHS = TSparseSpace::CreateEmptyVectorPointer();
            mpOV = TSparseSpace::CreateEmptyVectorPointer();

            mpKr = TReducedDenseSpace::CreateEmptyMatrixPointer();
            mpMr = TReducedDenseSpace::CreateEmptyMatrixPointer();
            mpRHSr = TReducedDenseSpace::CreateEmptyVectorPointer();
            mpOVr = TReducedDenseSpace::CreateEmptyVectorPointer();
            mpBasis = TReducedDenseSpace::CreateEmptyMatrixPointer();
            mpBasisLeft = TReducedDenseSpace::CreateEmptyMatrixPointer();
            mpCr = TReducedDenseSpace::CreateEmptyMatrixPointer();

            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    virtual void Clear() override
    {
        KRATOS_TRY;

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpK != nullptr)
            SparseSpaceType::Clear(mpK);
        if (mpKi != nullptr)
            SparseSpaceType::Clear(mpKi);
        if (mpM != nullptr)
            SparseSpaceType::Clear(mpM);
        if (mpMi != nullptr)
            SparseSpaceType::Clear(mpMi);
        if (mpC != nullptr)
            SparseSpaceType::Clear(mpC);
        if (mpRHS != nullptr)
            SparseSpaceType::Clear(mpRHS);
        if (mpOV != nullptr)
            SparseSpaceType::Clear(mpOV);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep() override
    {
        KRATOS_TRY;
        // std::cout << "offline initialize base class call\n";

        if (mSolutionStepIsInitialized == false)
        {
            KRATOS_WATCH(mUseDamping)
            // std::cout << "offline initialize base class perform\n";
            ModelPart& r_model_part = BaseType::GetModelPart();

            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

            const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            TSystemMatrixType& r_K  = *mpK;
            TSystemMatrixType& r_Ki  = *mpKi;
            TSystemMatrixType& r_M = *mpM;
            TSystemMatrixType& r_Mi = *mpMi;
            TSystemMatrixType& r_C = *mpC;
            TSystemVectorType& r_RHS  = *mpRHS;
            TSystemVectorType& r_OV  = *mpOV;

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            if ((p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true) &&
                mImportedSystem == false)
            {
                //setting up the list of the DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                //shaping correctly the system
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_system_time.ElapsedSeconds() << std::endl;

                //setting up the Vectors involved to the correct size
                BuiltinTimer system_matrix_resize_time;
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpK, mpOV, mpRHS,
                                                                 BaseType::GetModelPart());
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpM, mpOV, mpRHS,
                                                                 BaseType::GetModelPart());
                if( mUseDamping )
                {
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpC, mpOV, mpRHS,
                                                                    BaseType::GetModelPart());
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpKi, mpOV, mpRHS,
                                                                    BaseType::GetModelPart());
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpMi, mpOV, mpRHS,
                                                                    BaseType::GetModelPart());
                }

                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;

                //set up system matrices
                std::vector<unsigned int> fixed_dofs;
                DirichletUtility::GetDirichletConstraints<typename TBuilderAndSolverType::DofsArrayType>(this->GetBuilderAndSolver()->GetDofSet(), fixed_dofs);

                TSystemVectorType tmp(r_K.size1(), 0.0);

                //set up the stiffness matrix and rhs
                r_model_part.GetProcessInfo()[BUILD_LEVEL] = 1;
                TSparseSpace::SetToZero(r_RHS);
                p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_K, r_RHS);
                DirichletUtility::ApplyDirichletConditions<TSparseSpace>(r_K, r_RHS, fixed_dofs, 1.0);

                //set up the mass matrix
                r_model_part.GetProcessInfo()[BUILD_LEVEL] = 201;
                p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_M, tmp);
                DirichletUtility::ApplyDirichletConditions<TSparseSpace>(r_M, tmp, fixed_dofs, 0.0);

                //set up the damping matrix and imaginary components of stiffness and mass
                if( mUseDamping )
                {
                    r_model_part.GetProcessInfo()[BUILD_LEVEL] = 101;
                    p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_C, tmp);
                    DirichletUtility::ApplyDirichletConditions<TSparseSpace>(r_C, tmp, fixed_dofs, 0.0);

                    r_model_part.GetProcessInfo()[BUILD_LEVEL] = 111;
                    p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_Ki, tmp);
                    DirichletUtility::ApplyDirichletConditions<TSparseSpace>(r_Ki, tmp, fixed_dofs, 0.0);

                    r_model_part.GetProcessInfo()[BUILD_LEVEL] = 121;
                    p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_Mi, tmp);
                    DirichletUtility::ApplyDirichletConditions<TSparseSpace>(r_Mi, tmp, fixed_dofs, 0.0);
                }

                //create output vector
                r_model_part.GetProcessInfo()[BUILD_LEVEL] = 301;
                TSparseSpace::SetToZero(r_OV);
                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, r_OV);
                r_OV -= r_RHS;
            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;

            //initial operations ... things that are constant over the Solution Step
            // p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rA, rOV, rRHS);
            // p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rM, rOV, rRHS);

            // //initial operations ... things that are constant over the Solution Step
            // p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rA, rOV, rRHS);
            // p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rM, rOV, rRHS);


            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     * @todo Add support for unsymmetric problems
     */
    virtual void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        // project the system matrices onto the reduced subspace
        TSystemMatrixType& r_K  = *mpK;
        TSystemMatrixType& r_M = *mpM;
        TSystemVectorType& r_RHS  = *mpRHS;
        TSystemMatrixType& r_D  = *mpC;

        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_stiffness_matrix_reduced = this->GetKr();
        auto& r_mass_matrix_reduced = this->GetMr();
        auto& r_damping_matrix_reduced = this->GetDr();
        auto& r_output_vector = this->GetOutputVector();
        auto& r_output_vector_r = this->GetOVr();
        auto& r_basis = this->GetBasis();

        const size_t reduced_system_size = r_basis.size2();

        BuiltinTimer system_projection_time;

        r_force_vector_reduced.resize( reduced_system_size, false);
        noalias(r_force_vector_reduced) = prod( r_RHS, r_basis );

        r_output_vector_r.resize( reduced_system_size, false);
        noalias(r_output_vector_r) = prod( r_output_vector, r_basis );

        ProjectMatrix<TSystemMatrixType>(r_K, r_basis, r_stiffness_matrix_reduced);
        ProjectMatrix<TSystemMatrixType>(r_M, r_basis, r_mass_matrix_reduced);
        if (mUseDamping)
            ProjectMatrix<TSystemMatrixType>(r_D, r_basis, r_damping_matrix_reduced);

        KRATOS_INFO_IF("System Projection Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << system_projection_time.ElapsedSeconds() << std::endl;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), r_K, r_RHS, r_RHS);
        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), r_M, r_RHS, r_RHS);
        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), r_D, r_RHS, r_RHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), r_K, r_RHS, r_RHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), r_M, r_RHS, r_RHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), r_D, r_RHS, r_RHS);

        //Cleaning memory after the solution
        p_scheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            SparseSpaceType::Clear(mpK);
            SparseSpaceType::Clear(mpRHS);
            SparseSpaceType::Clear(mpM);
            SparseSpaceType::Clear(mpC);
            SparseSpaceType::Clear(mpOV);

            this->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_ERROR << "ERROR: Calling MOR offline strategy base class" << std::endl;
		return false;
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    virtual int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")

    }

    void ImportSystem(TSystemMatrixType& rK,
        TSystemMatrixType& rC,
        TSystemMatrixType& rM,
        TSystemVectorType& rRHS,
        TSystemVectorType& rOV)
    {
        Initialize();

        TSystemMatrixType& K = *mpK;
        K = TSystemMatrixType(rK);
        mpKi->resize(K.size1(), K.size2(), false);

        TSystemMatrixType& C = *mpC;
        C = TSystemMatrixType(rC);

        TSystemMatrixType& M = *mpM;
        M = TSystemMatrixType(rM);
        mpMi->resize(M.size1(), M.size2(), false);

        TSystemVectorType& RHS = *mpRHS;
        RHS = TSystemVectorType(rRHS);

        TSystemVectorType& OV = *mpOV;
        OV = TSystemVectorType(rOV);

        mImportedSystem = true;
    }

    void ImportSystem(TSystemMatrixType& rK,
        TSystemMatrixType& rKi,
        TSystemMatrixType& rC,
        TSystemMatrixType& rM,
        TSystemVectorType& rRHS,
        TSystemVectorType& rOV)
    {
        Initialize();
        std::cout << "hallo?\n";
        TSystemMatrixType& K = *mpK;
        K = TSystemMatrixType(rK);

        TSystemMatrixType& Ki = *mpKi;
        Ki = TSystemMatrixType(rKi);

        TSystemMatrixType& C = *mpC;
        C = TSystemMatrixType(rC);

        TSystemMatrixType& M = *mpM;
        M = TSystemMatrixType(rM);

        TSystemMatrixType& Mi = *mpMi;
        Mi.resize(M.size1(), M.size2(), false);

        TSystemVectorType& RHS = *mpRHS;
        RHS = TSystemVectorType(rRHS);

        TSystemVectorType& OV = *mpOV;
        OV = TSystemVectorType(rOV);

        mImportedSystem = true;
    }

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const std::string prefix)
    {
        TSystemMatrixType& rA = *mpK;
        TSystemMatrixType& rM = *mpM;
        TSystemMatrixType& rD = *mpC;
        TSystemVectorType& rRHS = *mpRHS;
        TSystemVectorType& rOV = *mpOV;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("Dx")  << "Mass Matrix = " << mpM << std::endl;
            KRATOS_INFO("Sx")  << "Damping Matrix = " << mpC << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        std::stringstream matrix_market_name;
        matrix_market_name << prefix << "_K.mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

        std::stringstream matrix_market_mass_name;
        matrix_market_mass_name << prefix << "_M.mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false);

        std::stringstream matrix_market_damping_name;
        matrix_market_damping_name << prefix << "_D.mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rD, false);

        std::stringstream matrix_market_vectname;
        matrix_market_vectname << prefix << "_RHS.mm";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rRHS);

        std::stringstream matrix_market_outvectname;
        matrix_market_outvectname << prefix << "_OV.mm";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_outvectname.str()).c_str(), rOV);
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

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    virtual TSystemMatrixType& GetK()
    {
        TSystemMatrixType& mK = *mpK;

        return mK;
    }

    virtual TSystemMatrixPointerType& pGetK()
    {
        return mpK;
    }

    virtual TSystemMatrixType& GetKi()
    {
        TSystemMatrixType& mKi = *mpKi;

        return mKi;
    }

    virtual TSystemMatrixPointerType& pGetKi()
    {
        return mpKi;
    }

    virtual TSystemMatrixType& GetM()
    {
        TSystemMatrixType& mM = *mpM;

        return mM;
    }

    virtual TSystemMatrixPointerType& pGetM()
    {
        return mpM;
    }

    virtual TSystemMatrixType& GetMi()
    {
        TSystemMatrixType& mMi = *mpMi;

        return mMi;
    }

    virtual TSystemMatrixPointerType& pGetMi()
    {
        return mpMi;
    }

    virtual TSystemMatrixType& GetC()
    {
        TSystemMatrixType& mS = *mpC;

        return mS;
    }

    virtual TSystemMatrixPointerType& pGetC()
    {
        return mpC;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    virtual TSystemVectorType& GetSystemVector() override
    {
        TSystemVectorType& mb = *mpRHS;

        return mb;
    }

    virtual TReducedDenseMatrixType& GetMr()
    {
        TReducedDenseMatrixType &mMr = *mpMr;

        return mMr;
    }

    virtual TReducedDenseMatrixPointerType& pGetMr()
    {
        return mpMr;
    }

    virtual TReducedDenseMatrixType& GetKr()
    {
        TReducedDenseMatrixType &mAr = *mpKr;

        return mAr;
    }

    virtual TReducedDenseMatrixPointerType& pGetKr()
    {
        return mpKr;
    }

    virtual TReducedDenseMatrixType& GetDr()
    {
        TReducedDenseMatrixType &mSr = *mpCr;

        return mSr;
    };

    virtual TReducedDenseMatrixPointerType& pGetDr()
    {
        return mpCr;
    }

    virtual TReducedDenseVectorType& GetRHSr()
    {

        TReducedDenseVectorType& mb = *mpRHSr;

        return mb;
    };

    virtual TReducedDenseVectorPointerType& pGetRHSr()
    {
        return mpRHSr;
    }

    // virtual TReducedSparseMatrixType& GetBasis()
    virtual TReducedDenseMatrixType& GetBasis()
    {
        TReducedDenseMatrixType &mBasis = *mpBasis;

        return mBasis;
    }

    virtual TReducedDenseMatrixPointerType& pGetBasis()
    {
        return mpBasis;
    }

    virtual TReducedDenseMatrixType& GetBasisLeft()
    {
        if( mSystemIsSymmetric ) {
            return *mpBasis;
        } else {
            return *mpBasisLeft;
        }
    }

    virtual TReducedDenseMatrixPointerType& pGetBasisLeft()
    {
        if( mSystemIsSymmetric ) {
            return mpBasis;
        } else {
            return mpBasisLeft;
        }
    }

    virtual TSystemVectorType& GetOutputVector()
    {
        TSystemVectorType &mOutputVector = *mpOV;

        return mOutputVector;
    }

    virtual TReducedDenseVectorType& GetOVr()
    {
        TReducedDenseVectorType &mOVr = *mpOVr;

        return mOVr;
    };

    std::size_t GetSystemSize()
    {
        if( mImportedSystem )
        {
            return mpK->size1();
        }
        else
        {
            return this->GetBuilderAndSolver()->GetEquationSystemSize();
        }
    }

    bool SystemIsSymmetric()
    {
        return mSystemIsSymmetric;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}

  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employe
    typename TLinearSolverType::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename TLinearSolverType::Pointer mpAdjointLinearSolver; /// The pointer to the linear solver for the adjoint problem considered

    TSystemVectorPointerType mpRHS; /// The RHS vector of the system of equations
    TSystemMatrixPointerType mpK; /// The Stiffness matrix of the system of equations
    TSystemMatrixPointerType mpKi; /// The imaginary part of the stiffness matrix of the system of equations
    TSystemMatrixPointerType mpM; /// The Mass matrix of the system of equations
    TSystemMatrixPointerType mpMi; /// The imaginary part of the Mass matrix of the system of equations
    TSystemMatrixPointerType mpC; /// The Damping matrix of the system of equations
    TSystemVectorPointerType mpOV; /// The output vector for single output

    // reduced matrices
    TReducedDenseVectorPointerType mpRHSr; //reduced RHS
    TReducedDenseMatrixPointerType mpKr;
    TReducedDenseMatrixPointerType mpMr;
    TReducedDenseMatrixPointerType mpBasis;
    TReducedDenseMatrixPointerType mpBasisLeft;
    TReducedDenseMatrixPointerType mpCr;
    TReducedDenseVectorPointerType mpOVr;

    // the reduced system size
    std::size_t mSystemSizeR;

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mUseDamping; /// Flag to set if damping should be considered

    bool mSystemIsSymmetric; /// Flag to set if the system is symmetric

    bool mImportedSystem; // Flag to set if an imported system is considered

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method prints information after reach the max number of iterations
     * @todo Replace by logger
     */

    virtual void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief computes X=W^H * A * V
     */
    template <typename MatrixType>
    void ProjectMatrix(MatrixType& rA, ComplexMatrix& rV, ComplexMatrix& rW, ComplexMatrix& rX)
    {
        const size_t system_size = rV.size1();
        const size_t reduced_system_size = rV.size2();
        ComplexVector v_col = ComplexZeroVector(system_size);
        ComplexVector tmp_1 = ComplexZeroVector(system_size);
        ComplexVector tmp_2 = ComplexZeroVector(reduced_system_size);
        ComplexMatrix rWH = ComplexZeroMatrix(reduced_system_size, system_size);
        noalias(rWH) = herm(rW);

        #pragma omp parallel for firstprivate(v_col, tmp_1, tmp_2) schedule(static)
        for( int i=0; i<static_cast<int>(reduced_system_size); ++i )
        {
            noalias(v_col) = column(rV,i);

            axpy_prod(rA, v_col, tmp_1, true);      // A*V  (=T)
            axpy_prod(rWH, tmp_1, tmp_2, true);     // W^H * T
            noalias(column(rX, i)) = tmp_2;         // X = W^H A V
        }
    }

    /**
     * @brief computes X=W^T * A * V
     */
    template <typename MatrixType>
    void ProjectMatrix(MatrixType& rA, Matrix& rV, Matrix& rW, Matrix& rX)
    {
        const size_t reduced_system_size = rV.size2();
        Vector v_col = ZeroVector(rV.size1());
        Vector tmp_1 = ZeroVector(rV.size1());
        Vector tmp_2 = ZeroVector(rV.size2());

        #pragma omp parallel for firstprivate(v_col, tmp_1, tmp_2) schedule(static)// nowait
        for( int i=0; i<static_cast<int>(reduced_system_size); ++i )
        {
            noalias(v_col) = column(rV,i);

            axpy_prod(rA, v_col, tmp_1, true);      // A*V  (=T)
            axpy_prod(tmp_1, rW, tmp_2, true);      // W' * T
            noalias(column(rX, i)) = tmp_2;         // X = W' A V
        }
    }

    /**
     * @brief computes X=V^H * A * V
     */
    template <typename MatrixType>
    void ProjectMatrix(MatrixType& rA, ComplexMatrix& rV, ComplexMatrix& rX)
    {
        const size_t system_size = rV.size1();
        const size_t reduced_system_size = rV.size2();
        ComplexVector v_col = ComplexZeroVector(system_size);
        ComplexVector tmp_1 = ComplexZeroVector(system_size);
        ComplexVector tmp_2 = ComplexZeroVector(reduced_system_size);
        ComplexMatrix rVH = ComplexZeroMatrix(reduced_system_size, system_size);
        noalias(rVH) = herm(rV);

        #pragma omp parallel for firstprivate(v_col, tmp_1, tmp_2) schedule(static)
        for( int i=0; i<static_cast<int>(reduced_system_size); ++i )
        {
            noalias(v_col) = column(rV,i);

            axpy_prod(rA, v_col, tmp_1, true);      // A*V  (=T)
            axpy_prod(rVH, tmp_1, tmp_2, true);     // V^H * T
            noalias(column(rX, i)) = tmp_2;         // X = V^H A V
        }
    }

    /**
     * @brief computes X=V^T * A * V
     */
    template <typename MatrixType>
    void ProjectMatrix(MatrixType& rA, Matrix& rV, Matrix& rX)
    {
        const size_t reduced_system_size = rV.size2();
        Vector v_col = ZeroVector(rV.size1());
        Vector tmp_1 = ZeroVector(rV.size1());
        Vector tmp_2 = ZeroVector(rV.size2());

        #pragma omp parallel for firstprivate(v_col, tmp_1, tmp_2) schedule(static)// nowait
        for( int i=0; i<static_cast<int>(reduced_system_size); ++i )
        {
            noalias(v_col) = column(rV,i);

            axpy_prod(rA, v_col, tmp_1, true);      // A*V  (=T)
            axpy_prod(tmp_1, rV, tmp_2, true);      // V' * T
            noalias(column(rX, i)) = tmp_2;         // X = V' A V
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    template <typename SparseSpaceType>
    void MatrixOutput(typename SparseSpaceType::MatrixType matrix, const std::string name)
    {
        std::stringstream matrix_market_name;
        matrix_market_name << name << ".mm";
        SparseSpaceType::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), matrix, false);
    }

    template <typename SparseSpaceType>
    void VectorOutput(typename SparseSpaceType::VectorType vector, const std::string name)
    {
        std::stringstream matrix_market_name;
        matrix_market_name << name << ".mm.rhs";
        SparseSpaceType::WriteMatrixMarketVector((char *)(matrix_market_name.str()).c_str(), vector);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * Copy constructor.
     */

    MorOfflineSecondOrderStrategy(const MorOfflineSecondOrderStrategy &Other){};

    ///@}

}; /* Class MorOfflineSecondOrderStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_OFFLINE_SECOND_ORDER_STRATEGY  defined */