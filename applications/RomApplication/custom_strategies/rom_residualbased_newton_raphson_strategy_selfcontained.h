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

#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "rom_application_variables.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"

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
 * @class RomResidualBasedNewtonRaphsonStrategy
 * @ingroup KratosCore
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 * @author Riccardo Rossi
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class RomResidualBasedNewtonRaphsonStrategy
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{

    // The size_t types
    typedef std::size_t SizeType;

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(RomResidualBasedNewtonRaphsonStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef RomResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

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

    typedef GlobalROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> GlobalROMBuilderAndSolverType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit RomResidualBasedNewtonRaphsonStrategy() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit RomResidualBasedNewtonRaphsonStrategy(ModelPart& rModelPart)
        : RomResidualBasedNewtonRaphsonStrategy(rModelPart, RomResidualBasedNewtonRaphsonStrategy::GetDefaultParameters())
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit RomResidualBasedNewtonRaphsonStrategy(ModelPart& rModelPart, Parameters ThisParameters)
        : BaseType(rModelPart),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            // Tells to the builder and solver if the reactions have to be Calculated or not
            p_builder_and_solver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

            // Tells to the Builder And Solver if the system matrix and vectors need to
            // be reshaped at each step or not
            p_builder_and_solver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
        } else {
            KRATOS_WARNING("RomResidualBasedNewtonRaphsonStrategy") << "BuilderAndSolver is not initialized. Please assign one before settings flags" << std::endl;
        }

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();
    }

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit RomResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, MoveMeshFlag),
          mpScheme(pScheme),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mReformDofSetAtEachStep(ReformDofSetAtEachStep),
          mCalculateReactionsFlag(CalculateReactions),
          mMaxIterationNumber(MaxIterations),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        KRATOS_TRY;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSolver));

        // Tells to the builder and solver if the reactions have to be Calculated or not
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit RomResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, MoveMeshFlag),
          mpScheme(pScheme),
          mpBuilderAndSolver(pNewBuilderAndSolver),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mReformDofSetAtEachStep(ReformDofSetAtEachStep),
          mCalculateReactionsFlag(CalculateReactions),
          mMaxIterationNumber(MaxIterations),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        KRATOS_TRY

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // Tells to the builder and solver if the reactions have to be Calculated or not
        p_builder_and_solver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        p_builder_and_solver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    KRATOS_DEPRECATED_MESSAGE("Constructor deprecated, please use the constructor without linear solver")
    explicit RomResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : RomResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY

        KRATOS_WARNING("RomResidualBasedNewtonRaphsonStrategy") << "This constructor is deprecated, please use the constructor without linear solver" << std::endl;

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // We check if the linear solver considered for the builder and solver is consistent
        auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
        KRATOS_ERROR_IF(p_linear_solver != pNewLinearSolver) << "Inconsistent linear solver in strategy and builder and solver. Considering the linear solver assigned to builder and solver :\n" << p_linear_solver->Info() << "\n instead of:\n" << pNewLinearSolver->Info() << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Constructor with Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param Settings Settings used in the strategy
     */
    RomResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters Settings)
        : BaseType(rModelPart),
          mpScheme(pScheme),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        KRATOS_TRY;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSolver));

        // Tells to the builder and solver if the reactions have to be Calculated or not
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Constructor specifying the builder and solver and using Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param Settings Settings used in the strategy
     */
    RomResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters Settings)
        : BaseType(rModelPart),
          mpScheme(pScheme),
          mpBuilderAndSolver(pNewBuilderAndSolver),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        KRATOS_TRY
        // Validate and assign defaults
        Settings = this->ValidateAndAssignParameters(Settings, this->GetDefaultParameters());
        this->AssignSettings(Settings);
        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // Tells to the builder and solver if the reactions have to be Calculated or not
        p_builder_and_solver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        p_builder_and_solver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    /**
     * @brief Constructor specifying the builder and solver and using Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param Parameters Settings used in the strategy
     */
    KRATOS_DEPRECATED_MESSAGE("Constructor deprecated, please use the constructor without linear solver")
    RomResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters Settings)
        : RomResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, Settings)
    {
        KRATOS_TRY

        KRATOS_WARNING("RomResidualBasedNewtonRaphsonStrategy") << "This constructor is deprecated, please use the constructor without linear solver" << std::endl;

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // We check if the linear solver considered for the builder and solver is consistent
        auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
        KRATOS_ERROR_IF(p_linear_solver != pNewLinearSolver) << "Inconsistent linear solver in strategy and builder and solver. Considering the linear solver assigned to builder and solver :\n" << p_linear_solver->Info() << "\n instead of:\n" << pNewLinearSolver->Info() << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~RomResidualBasedNewtonRaphsonStrategy() override
    {
        // If the linear solver has not been deallocated, clean it before
        // deallocating mpA. This prevents a memory error with the the ML
        // solver (which holds a reference to it).
        // NOTE: The linear solver is hold by the B&S
        auto p_builder_and_solver = this->GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->Clear();
        }

        // Deallocating system vectors to avoid errors in MPI. Clear calls
        // TrilinosSpace::Clear for the vectors, which preserves the Map of
        // current vectors, performing MPI calls in the process. Due to the
        // way Python garbage collection works, this may happen after
        // MPI_Finalize has already been called and is an error. Resetting
        // the pointers here prevents Clear from operating with the
        // (now deallocated) vectors.
        mpA.reset();
        mpDx.reset();
        mpb.reset();

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
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
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
     * @brief This method sets the flag mFullUpdateFlag
     * @param UseOldStiffnessInFirstIterationFlag The flag that tells if
     */
    void SetUseOldStiffnessInFirstIterationFlag(bool UseOldStiffnessInFirstIterationFlag)
    {
        mUseOldStiffnessInFirstIteration = UseOldStiffnessInFirstIterationFlag;
    }

    /**
     * @brief This method returns the flag mFullUpdateFlag
     * @return The flag that tells if
     */
    bool GetUseOldStiffnessInFirstIterationFlag()
    {
        return mUseOldStiffnessInFirstIteration;
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
     * @brief This method sets the flag mMaxIterationNumber
     * @param MaxIterationNumber This is the maximum number of on linear iterations
     */
    void SetMaxIterationNumber(unsigned int MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }

    /**
     * @brief This method gets the flag mMaxIterationNumber
     * @return mMaxIterationNumber: This is the maximum number of on linear iterations
     */
    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
        mpConvergenceCriteria->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY
        const DataCommunicator &r_comm = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator();
        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
            Initialize();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        DofsArrayType& r_dof_set = GetBuilderAndSolver()->GetDofSet();

        GetScheme()->Predict(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

        // Applying constraints if needed
        auto& r_constraints_array = BaseType::GetModelPart().MasterSlaveConstraints();
        const int local_number_of_constraints = r_constraints_array.size();
        const int global_number_of_constraints = r_comm.SumAll(local_number_of_constraints);
        if(global_number_of_constraints != 0) {
            const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();

            block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
                rConstraint.ResetSlaveDofs(r_process_info);
            });
            block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
                rConstraint.Apply(r_process_info);
            });

            // The following is needed since we need to eventually compute time derivatives after applying
            // Master slave relations
            TSparseSpace::SetToZero(rDx);
            this->GetScheme()->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        }

        // Move the mesh if needed
        if (this->MoveMeshFlag() == true)
            BaseType::MoveMesh();

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TConvergenceCriteriaType::Pointer p_convergence_criteria = mpConvergenceCriteria;

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            //initialisation of the convergence criteria
            if (p_convergence_criteria->IsInitialized() == false)
                p_convergence_criteria->Initialize(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // Clearing the system of equations
        if (mpA != nullptr)
            SparseSpaceType::Clear(mpA);
        if (mpDx != nullptr)
            SparseSpaceType::Clear(mpDx);
        if (mpb != nullptr)
            SparseSpaceType::Clear(mpb);

        // Clearing scheme
        auto p_scheme = GetScheme();
        if (p_scheme != nullptr) {
            GetScheme()->Clear();
        }

        mInitializeWasPerformed = false;
        
        // Clear variables that in the past pertained to the builder&solver
        mSlaveIds.clear();
        mMasterIds.clear();
        mInactiveSlaveDofs.clear();
        mT.resize(0,0,false);
        mConstantVector.resize(0,false);

        KRATOS_CATCH("");
    }

    /**
     * @brief This should be considered as a "post solution" convergence check which is useful for coupled analysis - the convergence criteria used is the one used inside the "solve" step
     */
    bool IsConverged() override
    {
        KRATOS_TRY;

        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb = *mpb;

        TSystemVectorType b_rom(GetNumberOfROMModes(BaseType::GetModelPart()));

        TSparseSpace::SetToZero(b_rom);

        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            TSparseSpace::SetToZero(rb);
            TSparseSpace::SetToZero(b_rom);

            dynamic_cast<GlobalROMBuilderAndSolverType*>(GetBuilderAndSolver().get())->BuildRomRHS(GetScheme(), BaseType::GetModelPart(), rb, b_rom);
        }

        TSystemVectorType rDq = BaseType::GetModelPart().GetRootModelPart().GetValue(ROM_SOLUTION_INCREMENT);

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), rA, rDq, b_rom);

        KRATOS_CATCH("");
    }

    /**
     * @brief This operations should be called before printing the results when non trivial results
     * (e.g. stresses)
     * Need to be calculated given the solution of the step
     * @details This operations should be called only when needed, before printing as it can involve a non
     * negligible cost
     */
    void CalculateOutputData() override
    {
        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        GetScheme()->CalculateOutputData(BaseType::GetModelPart(),
                                         GetBuilderAndSolver()->GetDofSet(),
                                         rA, rDx, rb);
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        KRATOS_WATCH('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        KRATOS_WATCH('RomStrategy')
        KRATOS_WATCH('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')

        // Pointers needed in the solution
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        ModelPart& r_model_part = BaseType::GetModelPart();

        // Set up the system, operation performed just once unless it is required
        // to reform the dof set at each iteration
        BuiltinTimer system_construction_time;
        if (!p_builder_and_solver->GetDofSetIsInitializedFlag() || mReformDofSetAtEachStep)
        {
            // Setting up the list of the DOFs to be solved
            BuiltinTimer setup_dofs_time;
            p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
            KRATOS_INFO_IF("RomResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup Dofs Time: "
                << setup_dofs_time << std::endl;

            // Shaping correctly the system
            BuiltinTimer setup_system_time;
            p_builder_and_solver->SetUpSystem(r_model_part);
            KRATOS_INFO_IF("RomResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup System Time: "
                << setup_system_time << std::endl;

            // Setting up the Vectors involved to the correct size
            BuiltinTimer system_matrix_resize_time;
            p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb,
                                                                r_model_part);
            KRATOS_INFO_IF("RomResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Matrix Resize Time: "
                << system_matrix_resize_time << std::endl;
        }

        KRATOS_INFO_IF("RomResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Construction Time: "
            << system_construction_time << std::endl;

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;
        TSystemVectorType b_rom(GetNumberOfROMModes(r_model_part));


        // Initial operations ... things that are constant over the Solution Step
        p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);

        // Initial operations ... things that are constant over the Solution Step
        p_scheme->InitializeSolutionStep(r_model_part, rA, rDx, rb);

        // Initialisation of the convergence criteria
        if (mpConvergenceCriteria->GetActualizeRHSflag())
        {
            TSparseSpace::SetToZero(rb);
            TSparseSpace::SetToZero(b_rom);
            dynamic_cast<GlobalROMBuilderAndSolverType*>(GetBuilderAndSolver().get())->BuildRomRHS(GetScheme(), BaseType::GetModelPart(), rb, b_rom);
        }

        mpConvergenceCriteria->InitializeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, b_rom);

        if (mpConvergenceCriteria->GetActualizeRHSflag()) {
            TSparseSpace::SetToZero(rb);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& r_model_part = BaseType::GetModelPart();

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        p_scheme->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        p_builder_and_solver->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        //Cleaning memory after the solution
        p_scheme->Clean();

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            this->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
    * @brief Gets the current solution vector
    * @param rDofSet The set of degrees of freedom (DOFs)
    * @param this_solution The vector that will be filled with the current solution values
    * @details This method retrieves the current solution values for the provided DOF set.
    * The provided solution vector will be resized to match the size of the DOF set if necessary,
    * and will be filled with the solution values corresponding to each DOF. Each value is accessed
    * using the equation ID associated with each DOF.
    */
    void GetCurrentSolution(DofsArrayType& rDofSet, Vector& this_solution) {
        this_solution.resize(rDofSet.size());
        block_for_each(rDofSet, [&](const auto& r_dof) {
            this_solution[r_dof.EquationId()] = r_dof.GetSolutionStepValue();
        });
    }

    /**
     * @brief Gets the matrix of non-converged solutions and the DOF set
     * @return A tuple containing the matrix of non-converged solutions and the DOF set
     * @details This method returns a tuple where the first element is a matrix of non-converged solutions and the second element is the corresponding DOF set.
     */
    std::tuple<Matrix, DofsArrayType> GetNonconvergedSolutions(){
        return std::make_tuple(mNonconvergedSolutionsMatrix, GetBuilderAndSolver()->GetDofSet());
    }

    /**
    * @brief Sets the state for storing non-converged solutions.
    * @param state The state to set for storing non-converged solutions (true to enable, false to disable).
    * @details This method enables or disables the storage of non-converged solutions at each iteration
    * by setting the appropriate flag. When the flag is set to true, non-converged solutions will be stored.
    */
    void SetUpNonconvergedSolutionsFlag(bool state) {
        mStoreNonconvergedSolutionsFlag = state;
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();
        std::vector<Vector> NonconvergedSolutions;

        if (mStoreNonconvergedSolutionsFlag) {
            Vector initial;
            GetCurrentSolution(r_dof_set,initial);
            NonconvergedSolutions.push_back(initial);
        }

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;
        TSystemVectorType b_rom(GetNumberOfROMModes(r_model_part));

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            if (mUseOldStiffnessInFirstIteration){
                p_builder_and_solver->BuildAndSolveLinearizedOnPreviousIteration(p_scheme, r_model_part, rA, rDx, rb,BaseType::MoveMeshFlag());
            } else {
                p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
            }
        } else {
            TSparseSpace::SetToZero(rDx);  // Dx = 0.00;
            TSparseSpace::SetToZero(rb);

            TSparseSpace::SetToZero(b_rom);
            dynamic_cast<GlobalROMBuilderAndSolverType*>(GetBuilderAndSolver().get())->BuildRomRHS(GetScheme(), BaseType::GetModelPart(), rb, b_rom);
        }

        // Debugging info
        EchoInfo(iteration_number);

        // Updating the results stored in the database
        UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (mStoreNonconvergedSolutionsFlag) {
            Vector first;
            GetCurrentSolution(r_dof_set,first);
            NonconvergedSolutions.push_back(first);
        }

        if (is_converged) {

            if (mpConvergenceCriteria->GetActualizeRHSflag()) {

                TSparseSpace::SetToZero(rb);


                TSparseSpace::SetToZero(b_rom);

                dynamic_cast<GlobalROMBuilderAndSolverType*>(GetBuilderAndSolver().get())->BuildRomRHS(GetScheme(), BaseType::GetModelPart(), rb, b_rom);
            }

            TSystemVectorType rDq = BaseType::GetModelPart().GetRootModelPart().GetValue(ROM_SOLUTION_INCREMENT);

            is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDq, b_rom);
        }

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false &&
               iteration_number++ < mMaxIterationNumber)
        {
            //setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (GetKeepSystemConstantDuringIterations() == false)
                    {
                        //A = 0.00;
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
            }
            else
            {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            EchoInfo(iteration_number);

            // Updating the results stored in the database
            UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            if (mStoreNonconvergedSolutionsFlag == true){
                Vector ith;
                GetCurrentSolution(r_dof_set,ith);
                NonconvergedSolutions.push_back(ith);
            }

            residual_is_updated = false;

            if (is_converged == true)
            {
                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(rb);


                    TSparseSpace::SetToZero(b_rom);

                    dynamic_cast<GlobalROMBuilderAndSolverType*>(GetBuilderAndSolver().get())->BuildRomRHS(GetScheme(), BaseType::GetModelPart(), rb, b_rom);
                    residual_is_updated = true;
                }

                TSystemVectorType rDq = BaseType::GetModelPart().GetRootModelPart().GetValue(ROM_SOLUTION_INCREMENT);

                is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDq, b_rom);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber) {
            MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("RomResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << mMaxIterationNumber << " iterations" << std::endl;
        }

        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (residual_is_updated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(mb);
            //    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, mb);
        }

        //calculate reactions if required
        if (mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        if (mStoreNonconvergedSolutionsFlag) {
            mNonconvergedSolutionsMatrix = Matrix( r_dof_set.size(), NonconvergedSolutions.size() );
            for (std::size_t i = 0; i < NonconvergedSolutions.size(); ++i) {
                block_for_each(r_dof_set, [&](const auto& r_dof) {
                    mNonconvergedSolutionsMatrix(r_dof.EquationId(), i) = NonconvergedSolutions[i](r_dof.EquationId());
                });
            }
        }

        return is_converged;
    }

    /**
     * @brief Function to get the number of modes being used in the ROM solver.
     * @details It returns the information from the NUMBER_OF_ROM_MODES global variable form the RomApplication, which is set in the Builder&Solver.
     */
    SizeType GetNumberOfROMModes(ModelPart &rModelPart) const noexcept
    {
        auto& r_root_mp = rModelPart.GetRootModelPart();
        return r_root_mp.GetValue(NUMBER_OF_ROM_MODES);
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetScheme()->Check(BaseType::GetModelPart());

        mpConvergenceCriteria->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                : "newton_raphson_strategy",
            "echo_level"                          : 1
            "use_old_stiffness_in_first_iteration": false,
            "max_iteration"                       : 10,
            "reform_dofs_at_each_step"            : false,
            "compute_reactions"                   : false,
            "builder_and_solver_settings"         : {},
            "convergence_criteria_settings"       : {},
            "linear_solver_settings"              : {},
            "scheme_settings"                     : {}
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "newton_raphson_strategy";
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
    TSystemMatrixType &GetSystemMatrix() override
    {
        TSystemMatrixType &mA = *mpA;

        return mA;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    TSystemVectorType& GetSystemVector() override
    {
        TSystemVectorType& mb = *mpb;

        return mb;
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    TSystemVectorType& GetSolutionVector() override
    {
        TSystemVectorType& mDx = *mpDx;

        return mDx;
    }

    /**
     * @brief Set method for the flag mKeepSystemConstantDuringIterations
     * @param Value If we consider constant the system of equations during the iterations
     */
    void SetKeepSystemConstantDuringIterations(bool Value)
    {
        mKeepSystemConstantDuringIterations = Value;
    }

    /**
     * @brief Get method for the flag mKeepSystemConstantDuringIterations
     * @return True if we consider constant the system of equations during the iterations, false otherwise
     */
    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "RomResidualBasedNewtonRaphsonStrategy";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}


    //#########################################
    //ResidualBased Builder and Solver
    //#########################################

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param b The RHS vector
     */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // Assemble all elements
        const auto timer = BuiltinTimer();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++) {
                auto it_elem = el_begin + k;

                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_elem, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }

            }

            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++) {
                auto it_cond = cond_begin + k;

                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_cond, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >= 1) << "Build time: " << timer << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 2) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation chosen the size of the matrix could
     * be equal to the total number of Dofs or to the number of unrestrained dofs
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA
        ) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const auto it_elem_begin = rModelPart.ElementsBegin();
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        // Contributions to the system
        LocalSystemMatrixType lhs_contribution(0, 0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_id;

        // Assemble all elements
        const auto timer = BuiltinTimer();

        #pragma omp parallel firstprivate(nelements, nconditions, lhs_contribution, equation_id )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; ++k) {
                auto it_elem = it_elem_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_elem, lhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHS(rA, lhs_contribution, equation_id);
                }
            }

            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; ++k) {
                auto it_cond = it_cond_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_cond, lhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHS(rA, lhs_contribution, equation_id);
                }
            }
        }

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >= 1) << "Build time LHS: " << timer << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 2) << "Finished parallel building LHS" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Build a rectangular matrix of size n*N where "n" is the number of unrestrained degrees of freedom
     * and "N" is the total number of degrees of freedom involved.
     * @details This matrix is obtained by building the total matrix without the lines corresponding to the fixed
     * degrees of freedom (but keeping the columns!!)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        TSystemVectorType tmp(A.size1(), 0.0);
        this->Build(pScheme, rModelPart, A, tmp);

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void SystemSolve(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY
        double norm_b;
        if (TSparseSpace::Size(rb) != 0) {
            norm_b = TSparseSpace::TwoNorm(rb);
        } else {
            norm_b = 0.0;
        }

        if (norm_b != 0.0) {
            // Do solve
            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        } else {
            TSparseSpace::SetToZero(rDx);
        }

        if(mT.size1() != 0) { // If there are master-slave constraints
            // Recover solution of the original problem
            TSystemVectorType Dxmodified = rDx;

            // Recover solution of the original problem
            TSparseSpace::Mult(mT, Dxmodified, rDx);
        }

        // Prints information about the current time
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    virtual void SystemSolveWithPhysics(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        ModelPart& rModelPart
        )
    {
        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            TSystemVectorType Dxmodified(rb.size());

            // Initialize the vector
            TSparseSpace::SetToZero(Dxmodified);

            InternalSystemSolveWithPhysics(rA, Dxmodified, rb, rModelPart);

            //recover solution of the original problem
            TSparseSpace::Mult(mT, Dxmodified, rDx);
        } else {
            InternalSystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        }
    }

    /**
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void InternalSystemSolveWithPhysics(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        ModelPart& rModelPart
    )
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(rb) != 0)
            norm_b = TSparseSpace::TwoNorm(rb);
        else
            norm_b = 0.00;

        if (norm_b != 0.00) {
            //provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(rA, rDx, rb, BaseType::mDofSet, rModelPart);

            //do solve
            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        } else {
            KRATOS_WARNING_IF("ResidualBasedBlockBuilderAndSolver", mOptions.IsNot(SILENT_WARNINGS)) << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // Prints information about the current time
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            const auto timer_constraints = BuiltinTimer();
            Timer::Start("ApplyConstraints");
            ApplyConstraints(pScheme, rModelPart, A, b);
            Timer::Stop("ApplyConstraints");
            KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >=1) << "Constraints build time: " << timer_constraints << std::endl;
        }

        ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >=1) << "System solve time: " << timer << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time Linearizing with the database at the old iteration
     * @details It is ideally the fastest and safer function to use when it is possible to solve just after building
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     * @param MoveMesh tells if the update of the scheme needs  to be performed when calling the Update of the scheme
     */
    void BuildAndSolveLinearizedOnPreviousIteration(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh
        ) override
    {
        Timer::Start("Linearizing on Old iteration");

        KRATOS_INFO_IF("BlockBuilderAndSolver", this->GetEchoLevel() > 0) << "Linearizing on Old iteration" << std::endl;

        KRATOS_ERROR_IF(rModelPart.GetBufferSize() == 1) << "BlockBuilderAndSolver: \n"
                << "The buffer size needs to be at least 2 in order to use \n"
                << "BuildAndSolveLinearizedOnPreviousIteration \n"
                << "current buffer size for modelpart: " << rModelPart.Name() << std::endl
                << "is :" << rModelPart.GetBufferSize()
                << " Please set IN THE STRATEGY SETTINGS "
                << " UseOldStiffnessInFirstIteration=false " << std::endl;

        DofsArrayType fixed_dofs;
        for(auto& r_dof : BaseType::mDofSet){
            if(r_dof.IsFixed()){
                fixed_dofs.push_back(&r_dof);
                r_dof.FreeDof();
            }
        }

        //TODO: Here we need to take the vector from other ones because
        // We cannot create a trilinos vector without a communicator. To be improved!
        TSystemVectorType dx_prediction(rDx);
        TSystemVectorType rhs_addition(rb); //we know it is zero here, so we do not need to set it

        // Here we bring back the database to before the prediction,
        // but we store the prediction increment in dx_prediction.
        // The goal is that the stiffness is computed with the
        // converged configuration at the end of the previous step.
        block_for_each(BaseType::mDofSet, [&](Dof<double>& rDof){
            // NOTE: this is initialized to - the value of dx prediction
            dx_prediction[rDof.EquationId()] = -(rDof.GetSolutionStepValue() - rDof.GetSolutionStepValue(1));

        });

        // Use UpdateDatabase to bring back the solution to how it was at the end of the previous step
        pScheme->Update(rModelPart, BaseType::mDofSet, rA, dx_prediction, rb);
        if (MoveMesh) {
            VariableUtils().UpdateCurrentPosition(rModelPart.Nodes(),DISPLACEMENT,0);
        }

        Timer::Stop("Linearizing on Old iteration");

        Timer::Start("Build");

        this->Build(pScheme, rModelPart, rA, rb);

        Timer::Stop("Build");

        // Put back the prediction into the database
        TSparseSpace::InplaceMult(dx_prediction, -1.0); //change sign to dx_prediction
        TSparseSpace::UnaliasedAdd(rDx, 1.0, dx_prediction);

        // Use UpdateDatabase to bring back the solution
        // to where it was taking into account BCs
        // it is done here so that constraints are correctly taken into account right after
        pScheme->Update(rModelPart, BaseType::mDofSet, rA, dx_prediction, rb);
        if (MoveMesh) {
            VariableUtils().UpdateCurrentPosition(rModelPart.Nodes(),DISPLACEMENT,0);
        }

        // Apply rb -= A*dx_prediction
        TSparseSpace::Mult(rA, dx_prediction, rhs_addition);
        TSparseSpace::UnaliasedAdd(rb, -1.0, rhs_addition);

        for(auto& dof : fixed_dofs)
            dof.FixDof();

        if (!rModelPart.MasterSlaveConstraints().empty()) {
            const auto timer_constraints = BuiltinTimer();
            Timer::Start("ApplyConstraints");
            this->ApplyConstraints(pScheme, rModelPart, rA, rb);
            Timer::Stop("ApplyConstraints");
            KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >=1) << "Constraints build time: " << timer_constraints << std::endl;
        }
        this->ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();

        Timer::Start("Solve");

        this->SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >=1) << "System solve time: " << timer << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, rb);

        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyRHSConstraints");
            ApplyRHSConstraints(pScheme, rModelPart, rb);
            Timer::Stop("ApplyRHSConstraints");
        }

        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >=1) << "System solve time: " << timer << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        Timer::Start("BuildRHS");

        BuildRHSNoDirichlet(pScheme,rModelPart,b);

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& rDof){
            const std::size_t i = rDof.EquationId();

            if (rDof.IsFixed())
                b[i] = 0.0;
        });

        Timer::Stop("BuildRHS");

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
     * and condition its Dofs.
     * @details The list of dofs is stores inside the BuilderAndSolver as it is closely connected to the
     * way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
    ) override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the dofs" << std::endl;

        //Gets the array of elements from the modeler
        ElementsArrayType& r_elements_array = rModelPart.Elements();
        const int number_of_elements = static_cast<int>(r_elements_array.size());

        DofsVectorType dof_list, second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        unsigned int nthreads = ParallelUtilities::GetNumThreads();

        typedef std::unordered_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Number of threads" << nthreads << "\n" << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        /**
         * Here we declare three sets.
         * - The global set: Contains all the DoF of the system
         * - The slave set: The DoF that are not going to be solved, due to MPC formulation
         */
        set_type dof_global_set;
        dof_global_set.reserve(number_of_elements*20);

        #pragma omp parallel firstprivate(dof_list, second_dof_list)
        {
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // We cleate the temporal set and we reserve some space on them
            set_type dofs_tmp_set;
            dofs_tmp_set.reserve(20000);

            // Gets the array of elements from the modeler
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i) {
                auto it_elem = r_elements_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of conditions from the modeler
            ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
            const int number_of_conditions = static_cast<int>(r_conditions_array.size());
            #pragma omp for  schedule(guided, 512) nowait
            for (int i = 0; i < number_of_conditions; ++i) {
                auto it_cond = r_conditions_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of constraints from the modeler
            auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
            const int number_of_constraints = static_cast<int>(r_constraints_array.size());
            #pragma omp for  schedule(guided, 512) nowait
            for (int i = 0; i < number_of_constraints; ++i) {
                auto it_const = r_constraints_array.begin() + i;

                // Gets list of Dof involved on every element
                it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
            }
        }

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dof_global_set.size());
        for (auto it= dof_global_set.begin(); it!= dof_global_set.end(); it++)
        {
            Doftemp.push_back( *it );
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
        if (BaseType::GetCalculateReactionsFlag()) {
            for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                    KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " <<std::endl
                        << "Node : "<<dof_iterator->Id()<< std::endl
                        << "Dof : "<<(*dof_iterator)<<std::endl<<"Not possible to calculate reactions."<<std::endl;
            }
        }
#endif

        KRATOS_CATCH("");
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(
        ModelPart& rModelPart
    ) override
    {
        //int free_id = 0;
        BaseType::mEquationSystemSize = BaseType::mDofSet.size();

        IndexPartition<std::size_t>(BaseType::mDofSet.size()).for_each([&, this](std::size_t Index){
            typename DofsArrayType::iterator dof_iterator = this->mDofSet.begin() + Index;
            dof_iterator->SetEquationId(Index);
        });
    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
    ) override
    {
        KRATOS_TRY
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }

        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(pScheme, A, rModelPart);
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permitted."<<std::endl;
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, A, rModelPart);
            }
        }
        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        TSparseSpace::SetToZero(Dx);
        if (b.size() != BaseType::mEquationSystemSize) {
            b.resize(BaseType::mEquationSystemSize, false);
        }
        TSparseSpace::SetToZero(b);

        ConstructMasterSlaveConstraintsStructure(rModelPart);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        TSparseSpace::SetToZero(b);

        //refresh RHS to have the correct reactions
        BuildRHSNoDirichlet(pScheme, rModelPart, b);

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& rDof){
            const std::size_t i = rDof.EquationId();

            rDof.GetSolutionStepReactionValue() = -b[i];
        });
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation chosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver chosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        const std::size_t system_size = rA.size1();
        Vector scaling_factors (system_size);

        const auto it_dof_iterator_begin = BaseType::mDofSet.begin();

        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        IndexPartition<std::size_t>(BaseType::mDofSet.size()).for_each([&](std::size_t Index){
            auto it_dof_iterator = it_dof_iterator_begin + Index;
            if (it_dof_iterator->IsFixed()) {
                scaling_factors[Index] = 0.0;
            } else {
                scaling_factors[Index] = 1.0;
            }
        });

        // Detect if there is a line of all zeros and set the diagonal to a certain number (1 if not scale, some norms values otherwise) if this happens
        mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

        double* Avalues = rA.value_data().begin();
        std::size_t* Arow_indices = rA.index1_data().begin();
        std::size_t* Acol_indices = rA.index2_data().begin();

        IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
            const std::size_t col_begin = Arow_indices[Index];
            const std::size_t col_end = Arow_indices[Index+1];
            const double k_factor = scaling_factors[Index];
            if (k_factor == 0.0) {
                // Zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (Acol_indices[j] != Index )
                        Avalues[j] = 0.0;

                // Zero out the RHS
                rb[Index] = 0.0;
            } else {
                // Zero out the column which is associated with the zero'ed row
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        });
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix (RHS only)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS vector
     */
    void ApplyRHSConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            BuildMasterSlaveConstraints(rModelPart);

            // We compute the transposed matrix of the global relation matrix
            TSystemMatrixType T_transpose_matrix(mT.size2(), mT.size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, mT, 1.0);

            TSystemVectorType b_modified(rb.size());
            TSparseSpace::Mult(T_transpose_matrix, rb, b_modified);
            TSparseSpace::Copy(b_modified, rb);

            // Apply diagonal values on slaves
            IndexPartition<std::size_t>(mSlaveIds.size()).for_each([&](std::size_t Index){
                const IndexType slave_equation_id = mSlaveIds[Index];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    rb[slave_equation_id] = 0.0;
                }
            });
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            BuildMasterSlaveConstraints(rModelPart);

            // We compute the transposed matrix of the global relation matrix
            TSystemMatrixType T_transpose_matrix(mT.size2(), mT.size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, mT, 1.0);

            TSystemVectorType b_modified(rb.size());
            TSparseSpace::Mult(T_transpose_matrix, rb, b_modified);
            TSparseSpace::Copy(b_modified, rb);

            TSystemMatrixType auxiliar_A_matrix(mT.size2(), rA.size2());
            SparseMatrixMultiplicationUtility::MatrixMultiplication(T_transpose_matrix, rA, auxiliar_A_matrix); //auxiliar = T_transpose * rA
            T_transpose_matrix.resize(0, 0, false);                                                             //free memory

            SparseMatrixMultiplicationUtility::MatrixMultiplication(auxiliar_A_matrix, mT, rA); //A = auxilar * T   NOTE: here we are overwriting the old A matrix!
            auxiliar_A_matrix.resize(0, 0, false);                                              //free memory

            // Compute the scale factor value
            mScaleFactor = TSparseSpace::GetScaleNorm(rModelPart.GetProcessInfo(), rA, mScalingDiagonal);

            // Apply diagonal values on slaves
            IndexPartition<std::size_t>(mSlaveIds.size()).for_each([&](std::size_t Index){
                const IndexType slave_equation_id = mSlaveIds[Index];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    rA(slave_equation_id, slave_equation_id) = mScaleFactor;
                    rb[slave_equation_id] = 0.0;
                }
            });
        }

        KRATOS_CATCH("")
    }

  private:
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
  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    typename TSchemeType::Pointer mpScheme = nullptr; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver = nullptr; /// The pointer to the builder and solver employed
    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria = nullptr; /// The pointer to the convergence criteria employed

    TSystemVectorPointerType mpDx; /// The increment in the solution
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

    /**
     * @brief Flag telling if it is needed or not to compute the reactions
     * @details default = true
     */
    bool mCalculateReactionsFlag;

    /**
     * @brief Flag telling if a full update of the database will be performed at the first iteration
     * @details default = false
     */
    bool mUseOldStiffnessInFirstIteration = false;

    unsigned int mMaxIterationNumber; /// The maximum number of iterations, 30 by default

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mKeepSystemConstantDuringIterations; // Flag to allow keeping system matrix constant during iterations

    /**
     * @brief This matrix stores the non-converged solutions
     * @details The matrix is structured such that each column represents the solution vector at a specific non-converged iteration.
     */
    Matrix mNonconvergedSolutionsMatrix;

    /**
     * @brief Flag indicating whether to store non-converged solutions
     * @details Only when set to true (by calling the SetUpNonconvergedSolutionsGathering method) will the non-converged solutions at each iteration be stored.
     */
    bool mStoreNonconvergedSolutionsFlag = false;



    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief Here the database is updated
     * @param A The LHS matrix of the system of equations
     * @param Dx The incremement in the solution
     * @param b The RHS vector of the system of equations
     * @param MoveMesh The flag that allows to move the mesh
     */
    virtual void UpdateDatabase(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh)
    {
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        p_scheme->Update(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        // Move the mesh if needed
        if (MoveMesh == true)
            BaseType::MoveMesh();
    }

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber)
    {
        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("Dx")  << "Solution obtained = " << rDx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rb << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("Dx")  << "Solution obtained = " << rDx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rb << std::endl;
        }
        else if (this->GetEchoLevel() == 4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rb);

            std::stringstream matrix_market_dxname;
            matrix_market_dxname << "dx_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_dxname.str()).c_str(), rDx);

            std::stringstream dof_data_name;
            unsigned int rank=BaseType::GetModelPart().GetCommunicator().MyPID();
            dof_data_name << "dofdata_" << BaseType::GetModelPart().GetProcessInfo()[TIME]
                << "_" << IterationNumber << "_rank_"<< rank << ".csv";
            WriteDofInfo(dof_data_name.str(), rDx);
        }
    }

    /**
     * @brief This method prints information after reach the max number of iterations
     */
    virtual void MaxIterationsExceeded()
    {
        KRATOS_INFO_IF("RomResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
            << "ATTENTION: max iterations ( " << mMaxIterationNumber
            << " ) exceeded!" << std::endl;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
        mMaxIterationNumber = ThisParameters["max_iteration"].GetInt();
        mReformDofSetAtEachStep = ThisParameters["reform_dofs_at_each_step"].GetBool();
        mCalculateReactionsFlag = ThisParameters["compute_reactions"].GetBool();
        mUseOldStiffnessInFirstIteration = ThisParameters["use_old_stiffness_in_first_iteration"].GetBool();

        // Saving the convergence criteria to be used
        if (ThisParameters["convergence_criteria_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }

        // Saving the scheme
        if (ThisParameters["scheme_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }

        // Setting up the default builder and solver
        if (ThisParameters["builder_and_solver_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }

        AssignBuilderAndSolverSettings(ThisParameters)
    }

    /**
     * @brief Validate and assign defaults for the parameters of what used to be the Builder&Solver
     */
    explicit AssignBuilderAndSolverSettings(
        Parameters ThisParameters
        )
    {
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    void WriteDofInfo(std::string FileName, const TSystemVectorType& rDX)
    {
        std::ofstream out(FileName);

        out.precision(15);
        out << "EquationId,NodeId,VariableName,IsFixed,Value,coordx,coordy,coordz" << std::endl;
        for(const auto& rdof : GetBuilderAndSolver()->GetDofSet()) {
            const auto& coords = BaseType::GetModelPart().Nodes()[rdof.Id()].Coordinates();
            out << rdof.EquationId() << "," << rdof.Id() << "," << rdof.GetVariable().Name() << "," << rdof.IsFixed() << ","
                        << rdof.GetSolutionStepValue() << "," <<  "," << coords[0]  << "," << coords[1]  << "," << coords[2]<< "\n";
        }
        out.close();
    }

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
    ///@name Un accessible methods
    ///@{

    /**
     * Copy constructor.
     */

    RomResidualBasedNewtonRaphsonStrategy(const RomResidualBasedNewtonRaphsonStrategy &Other){};

    ///@}

}; /* Class RomResidualBasedNewtonRaphsonStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY  defined */