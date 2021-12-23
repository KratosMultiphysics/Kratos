//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Ignasi Pouplana
//

#if !defined(KRATOS_ARC_LENGTH_STRATEGY)
#define KRATOS_ARC_LENGTH_STRATEGY

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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
 * @class ArcLengthStrategy
 * @ingroup KratosCore
 * @brief This is the base ArcLengthStrategy
 * @details The theoretical details can be found in "Geometrical interpretation of the ARC-LENGTH method", M. Fafard and B. Massicotte,
 * Computers and Structures Vol 46 pp 603-615 (1993). (Ramm arc-length method)
 * @author Alejandro Cornejo
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ArcLengthStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
    public:
    ///@name Type Definitions
    ///@{

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(ArcLengthStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor
    ArcLengthStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewConvergenceCriteria, pNewBuilderAndSolver, rParameters, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            mDesiredIterations = rParameters["desired_iterations"].GetInt();
            mMaxRadiusFactor   = rParameters["max_radius_factor"].GetDouble();
            mMinRadiusFactor   = rParameters["min_radius_factor"].GetDouble();
            mInitializeArcLengthWasPerformed = false;
        }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~ArcLengthStrategy() override
    {
        BaseType::~BaseType();
    }


    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (!BaseType::mInitializeWasPerformed)
        {
            BaseType::Initialize();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // Setting to zero the internal flag to ensure that the dof sets are recalculated. Also clear the linear solver stored in the B&S
        auto p_builder_and_solver = GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->SetDofSetIsInitializedFlag(false);
            p_builder_and_solver->Clear();
        }

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
        mSolutionStepIsInitialized = false;

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

        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            TSparseSpace::SetToZero(rb);
            GetBuilderAndSolver()->BuildRHS(GetScheme(), BaseType::GetModelPart(), rb);
        }

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), rA, rDx, rb);

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

        if (!mSolutionStepIsInitialized) {
            // Pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
            ModelPart& r_model_part = BaseType::GetModelPart();

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
            {
                //setting up the list of the DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
                KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup Dofs Time: "
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                //shaping correctly the system
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(r_model_part);
                KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup System Time: "
                    << setup_system_time.ElapsedSeconds() << std::endl;

                //setting up the Vectors involved to the correct size
                BuiltinTimer system_matrix_resize_time;
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb,
                                                                 r_model_part);
                KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Matrix Resize Time: "
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;
            }

            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Construction Time: "
                << system_construction_time.ElapsedSeconds() << std::endl;

            TSystemMatrixType& rA  = *mpA;
            TSystemVectorType& rDx = *mpDx;
            TSystemVectorType& rb  = *mpb;

            // Initial operations ... things that are constant over the Solution Step
            p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);

            // Initial operations ... things that are constant over the Solution Step
            p_scheme->InitializeSolutionStep(r_model_part, rA, rDx, rb);

            // Initialisation of the convergence criteria
            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(rb);
                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            mpConvergenceCriteria->InitializeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                TSparseSpace::SetToZero(rb);

            mSolutionStepIsInitialized = true;
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

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            this->Clear();
        }

        KRATOS_CATCH("");
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

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

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

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }

        // Debugging info
        EchoInfo(iteration_number);

        // Updating the results stored in the database
        UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (is_converged) {
            if (mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
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

            residual_is_updated = false;

            if (is_converged == true)
            {
                if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber) {
            MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
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

        return is_converged;
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
        return "ResidualBasedNewtonRaphsonStrategy";
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

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    unsigned int mMaxIterationNumber; /// The maximum number of iterations, 30 by default

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mKeepSystemConstantDuringIterations; // Flag to allow keeping system matrix constant during iterations

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
        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
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

    ResidualBasedNewtonRaphsonStrategy(const ResidualBasedNewtonRaphsonStrategy &Other){};

    ///@}

}; /* Class ResidualBasedNewtonRaphsonStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_ARC_LENGTH_STRATEGY  defined */