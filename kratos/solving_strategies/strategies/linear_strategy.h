//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "schemes/new_scheme.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "utilities/builtin_timer.h"

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
 * @class LinearStrategy
 * @ingroup KratosCore
 * @brief This is a very simple strategy to solve linearly the problem
 * @details As a linear strategy the check on the convergence is not done and just one non linear iteration will be performed
 * @author Ruben Zorrilla
 */
//TODO: Add a base class
template<class TMatrixType, class TVectorType, class TLinearSolver>
class LinearStrategy
{
public:
    ///@name Type Definitions */
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(LinearStrategy);

    /// The definition of the current class
    using SolvingStrategyType = LinearStrategy<TMatrixType, TVectorType, TLinearSolver>; //TODO: This will be the eventual base class

    /// The definition of the current class
    using ClassType = LinearStrategy<TMatrixType, TVectorType, TLinearSolver>;

    /// Data type definition
    using DataType = typename TMatrixType::DataType;

    /// Index type definition
    using IndexType = typename TMatrixType::IndexType;

    /// Size type definition
    using SizeType = std::size_t;

    // Scheme type definition
    using SchemeType = NewScheme<TMatrixType, TVectorType>;

    // Scheme pointer type definition
    using SchemePointerType = typename SchemeType::Pointer;

    /// DoF array type definition
    using DofsArrayType = typename SchemeType::DofsArrayType;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Matrix pointer type definition
    using SystemMatrixPointerType = typename SystemMatrixType::Pointer;

    /// Vector type definition
    using SystemVectorType = TVectorType;

    /// Vector type definition
    using SystemVectorPointerType = typename SystemVectorType::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit LinearStrategy() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit LinearStrategy(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mpModelPart(rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // Set flags to start the calculations correctly
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        // FIXME:
        // GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        // FIXME:
        // GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param CalculateReactionFlag The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param CalculateNormDxFlag The flag sets if the norm of Dx is computed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit LinearStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool CalculateReactionFlag = false,
        bool ReformDofSetAtEachStep = false,
        bool CalculateNormDxFlag = false,
        bool MoveMeshFlag = false)
        : mpScheme(pScheme)
        , mpModelPart(rModelpart)
        , mReformDofSetAtEachStep(ReformDofSetAtEachStep)
        , mCalculateNormDxFlag(CalculateNormDxFlag)
        , mCalculateReactionsFlag(CalculateReactionFlag)
    {
        KRATOS_TRY

        // Set flag to start the calculations correctly
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        //FIXME:
        // Tells to the builder and solver if the reactions have to be Calculated or not
        // GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        //FIXME:
        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        // GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        this->SetEchoLevel(1);

        // By default the matrices are rebuilt at each solution step
        BaseType::SetRebuildLevel(1);

        KRATOS_CATCH("")
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~LinearStrategy()
    {
        //  If the linear solver has not been deallocated, clean it before
        //  deallocating mpA. This prevents a memory error with the the ML
        //  solver (which holds a reference to it).
        mpLinearSolver->Clear();

        //TODO: Update comment to current implementation
        // Deallocating system vectors to avoid errors in MPI. Clear calls
        // TrilinosSpace::Clear for the vectors, which preserves the Map of
        // current vectors, performing MPI calls in the process. Due to the
        // way Python garbage collection works, this may happen after
        // MPI_Finalize has already been called and is an error. Resetting
        // the pointers here prevents Clear from operating with the
        // (now deallocated) vectors.
        mpA.clear();
        mpdx.clear();
        mpb.clear();

        // Finally clear current class
        this->Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(SchemePointerType pScheme)
    {
        mpScheme = pScheme;
    }

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    SchemePointerType pGetScheme()
    {
        return mpScheme;
    }

    /**
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(const bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
        pGetScheme()->SetCalculateReactionsFlag(mCalculateReactionsFlag) //TODO: Most probably we can avoid saving the reactions flag in the scheme

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag() const
    {
        return mCalculateReactionsFlag;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param ReformDofSetAtEachStepFlag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(const bool ReformDofSetAtEachStepFlag)
    {
        mReformDofSetAtEachStep = ReformDofSetAtEachStepFlag;
        // FIXME: I don't think this is required to be set at the scheme
        // GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag() const
    {
        return mReformDofSetAtEachStep;
    }

    /**
     * @brief This method sets the flag mMoveMeshFlag
     * @param MoveMeshFlag The flag that tells if the mesh is to be updated
     */
    void SetMoveMeshFlag(const bool MoveMeshFlag)
    {
        mMoveMeshFlag = MoveMeshFlag;
    }

    /**
     * @brief This method returns the flag mMoveMeshFlag
     * @return The flag that tells if the mesh is to be updated
     */
    bool GetMoveMeshFlag() const
    {
        return mMoveMeshFlag;
    }

    //TODO: Move to the future base class
    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    };

    //TODO: Move to the future base class
    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    };

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int EchoLevel)
    {
        mEchoLevel = EchoLevel;
        pGetScheme()->SetEchoLevel(EchoLevel);
    }

    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict()
    {
        KRATOS_TRY

        // Internal solution loop check to avoid repetitions
        KRATOS_ERROR_IF(!mInitializeWasPerformed) << "Initialize needs to be performed. Call Initialize() once before the solution loop." << std::endl;
        KRATOS_ERROR_IF(!mSolutionStepIsInitialized) << "InitializeSolutionStep needs to be performed. Call InitializeSolutionStep() before Predict()." << std::endl;

        // Call the time scheme predict
        pGetScheme()->Predict(GetModelPart(), *mpA, *mpdx, *mpb);

        //TODO: Constraints implementation
        // // Applying constraints if needed
        // const auto &r_comm = GetModelPart().GetCommunicator().GetDataCommunicator();
        // auto& r_constraints_array = BaseType::GetModelPart().MasterSlaveConstraints();
        // const int local_number_of_constraints = r_constraints_array.size();
        // const int global_number_of_constraints = r_comm.SumAll(local_number_of_constraints);
        // if(global_number_of_constraints != 0) {
        //     const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();

        //     block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
        //         rConstraint.ResetSlaveDofs(r_process_info);
        //     });
        //     block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
        //         rConstraint.Apply(r_process_info);
        //     });

        //     //the following is needed since we need to eventually compute time derivatives after applying
        //     //Master slave relations
        //     TSparseSpace::SetToZero(rDx);
        //     this->pGetScheme()->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        // }

        // Mesh motion update
        if (GetMoveMeshFlag()) {
            MoveMesh();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize()
    {
        KRATOS_TRY

        if (!mInitializeWasPerformed) {
            // Get scheme
            auto p_scheme = pGetScheme();

            // Initialize scheme (this has to be done once)
            if (!p_scheme->SchemeIsInitialized()) {
                p_scheme->Initialize(GetModelPart());
            }

            // Initialize elements (this has to be done once)
            if (!p_scheme->ElementsAreInitialized()) {
                p_scheme->InitializeElements(GetModelPart()); //TODO: Think if this belongs to the scheme
            }

            // Initialize conditions (this has to be done once)
            if (!p_scheme->ConditionsAreInitialized()) {
                p_scheme->InitializeConditions(GetModelPart()); //TODO: Think if this belongs to the scheme
            }

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief The problem of interest is solved
     * @details a double containing norm(Dx) is returned if CalculateNormDxFlag == true, else 0 is returned
     * @return norm(Dx)
     */
    double Solve()
    {
        KRATOS_ERROR << "Solve method is no longer available. Call the solving steps one by one." << std::endl;

        return 0.0;
    }

    /**
     * @brief Clears the internal storage
     * @note NULL could be changed to nullptr in the future (c++11)
     */
    void Clear()
    {
        KRATOS_TRY;

        // Clear the linear solver
        mpLinearSolver->Clear();

        // Clearing the system of equations
        if (mpA != nullptr) {
            mpA->clear();
        }
        if (mpdx != nullptr) {
            mpdx->clear();
        }
        if (mpb != nullptr) {
            mpb->clear();
        }

        // Clearing scheme
        // Note that this resets the DOF set
        if (pGetScheme() != nullptr) {
            pGetScheme()->Clear();
        }

        // Reset initialization flags
        mInitializeWasPerformed = false;
        mSolutionStepIsInitialized = false;

        KRATOS_CATCH("");
    }

    /**
     * @brief This operations should be called before printing the results when non trivial results (e.g. stresses)
    need to be calculated given the solution of the step
     *@details This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
    void CalculateOutputData()
    {
        pGetScheme()->CalculateOutputData(GetModelPart(), *mpA, *mpdx, *mpb);
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     * @todo Boost dependencies should be replaced by std equivalent
     */
    void InitializeSolutionStep()
    {
        KRATOS_TRY

        if (!mSolutionStepIsInitialized) {
            // Get scheme pointer
            auto p_scheme = pGetScheme();

            // Set up the system
            // This operation is performed just once unless it is required to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            if (!mDofSetIsInitialized || mReformDofSetAtEachStep) {
                // Setting up the DOFs list to be solved
                BuiltinTimer setup_dofs_time;
                p_scheme->SetUpDofArray(mDofSet);
                mDofSetIsInitialized = true;
                KRATOS_INFO_IF("LinearStrategy", GetEchoLevel() > 0) << "Setup DOFs Time: " << setup_dofs_time << std::endl;

                // Shaping the system
                BuiltinTimer setup_system_time;
                p_scheme->SetUpSystem(GetModelPart());
                KRATOS_INFO_IF("LinearStrategy", GetEchoLevel() > 0) << "Set up system time: " << setup_system_time << std::endl;

                // Allocating the system vectors to their correct sizes
                BuiltinTimer system_matrix_resize_time;
                p_scheme->ResizeAndInitializeVectors(mpA, mpdx, mpb);
                KRATOS_INFO_IF("LinearStrategy", GetEchoLevel() > 0) << "System matrix resize time: " << system_matrix_resize_time << std::endl;
            }
            KRATOS_INFO_IF("LinearStrategy", GetEchoLevel() > 0) << "System construction time: " << system_construction_time << std::endl;

            // Call the scheme InitializeSolutionStep
            p_scheme->InitializeSolutionStep(GetModelPart(), *mpA, *mpDx, *mpb);

            // Set the flag to avoid calling this twice
            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep()
    {
        KRATOS_TRY;

        // Get scheme pointer
        auto p_scheme = pGetScheme();

        // Finalisation of the solution step (operations to be done after achieving convergence)
        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), *mpA, *mpdx, *mpb);

        // Cleaning memory after the solution
        p_scheme->Clean();

        // Reset flag for next step
        mSolutionStepIsInitialized = false;

        // Clear if needed (note that this deallocates the system arrays)
        if (mReformDofSetAtEachStep == true) {
            Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep()
    {
        // Get scheme pointer
        auto p_scheme = pGetScheme();

        // Initialize non-linear iteration (once as this is a linear strategy)
	    p_scheme->InitializeNonLinIteration(GetModelPart(), *mpA, *mpdx, *mpb);

        if (mRebuildLevel > 0 || mStiffnessMatrixIsBuilt == false) {
            // Initialize values
            mpA->SetValue(0.0);
            mpdx->SetValue(0.0);
            mpb->SetValue(0.0);

            p_scheme->Build(GetModelPart(), *mpA, *mpb);
            p_scheme->ApplyDirichletConditions(GetModelPart(), *mpA, *mpb);
            //TODO: Call the linear solver to solve in here!!!!


            p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
            mStiffnessMatrixIsBuilt = true;
        } else {
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);
            p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }

        // Debugging info
        EchoInfo();

        //update results
        p_scheme->Update(BaseType::GetModelPart(), GetDofSet(), rA, rDx, rb);

        //move the mesh if needed
        if (GetMoveMeshFlag()) {
            MoveMesh();
        }

        p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // Calculate reactions if required
        if (mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme,
                                                     BaseType::GetModelPart(),
                                                     rA, rDx, rb);

        return true;
    }

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    SystemMatrixType& GetSystemMatrix()
    {
        SystemMatrixType& mA = *mpA;

        return mA;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    SystemVectorType& GetSystemVector()
    {
        SystemVectorType& mb = *mpb;

        return mb;
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    SystemVectorType& GetSolutionVector()
    {
        SystemVectorType& mDx = *mpdx;

        return mDx;
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    DofsArrayType& GetDofSet()
    {
        return mDofSet;
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    const DofsArrayType& GetDofSet() const
    {
        return mDofSet;
    }

    /**
     * @brief This method returns the residual norm
     * @return The residual norm
     */
    double GetResidualNorm()
    {
        if (TSparseSpace::Size(*mpb) != 0)
            return TSparseSpace::TwoNorm(*mpb);
        else
            return 0.0;
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check()
    {
        KRATOS_TRY

        BaseType::Check();

        // FIXME:
        // GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        pGetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                         : "linear_strategy",
            "compute_norm_dx"              : false,
            "reform_dofs_at_each_step"     : false,
            "compute_reactions"            : false,
            "builder_and_solver_settings"  : {},
            "linear_solver_settings"       : {},
            "scheme_settings"              : {}
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
        return "linear_strategy";
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This returns the level of echo for the solving strategy
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic information
     * 2 -> Printing linear solver data
     * 3 -> Print of debug information: Echo of stiffness matrix, Dx, b...
     * }
     * @return Level of echo for the solving strategy
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "LinearStrategy";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

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

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters)
    {
        BaseType::AssignSettings(ThisParameters);
        mCalculateNormDxFlag = ThisParameters["compute_norm_dx"].GetBool();
        mReformDofSetAtEachStep = ThisParameters["reform_dofs_at_each_step"].GetBool();
        mCalculateReactionsFlag = ThisParameters["compute_reactions"].GetBool();

        // Saving the scheme
        if (ThisParameters["scheme_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }

    }

    //TODO: Move this to the future base class
    /**
     * @brief This function is designed to move the mesh
     * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
     */
    virtual void MoveMesh()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(GetModelPart().HasNodalSolutionStepVariable(DISPLACEMENT_X)) << "It is impossible to move the mesh since the DISPLACEMENT variable is not in the ModelPart. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;

        block_for_each(GetModelPart().Nodes(), [](Node& rNode){
            noalias(rNode.Coordinates()) = rNode.GetInitialPosition().Coordinates();
            noalias(rNode.Coordinates()) += rNode.FastGetSolutionStepValue(DISPLACEMENT);
        });

        KRATOS_INFO_IF("SolvingStrategy", this->GetEchoLevel() != 0) << "Mesh moved." << std::endl;

        KRATOS_CATCH("")
    }

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

    typename TSchemeType::Pointer mpScheme = nullptr; /// The pointer to the linear solver considered

    DofsArrayType mDofSet; /// The set containing the DoF of the system

    SystemVectorPointerType mpdx; /// The incremement in the solution

    SystemVectorPointerType mpb; /// The RHS vector of the system of equations

    SystemMatrixPointerType mpA; /// The LHS matrix of the system of equations

    int mEchoLevel;

    ModelPart* mpModelPart = nullptr;

    bool mMoveMeshFlag;

    bool mDofSetIsInitialized = false;

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
     * @brief This method returns the components of the system of equations depending of the echo level
     */
    virtual void EchoInfo()
    {
        SystemMatrixType& rA  = *mpA;
        SystemVectorType& rDx = *mpdx;
        SystemVectorType& rb  = *mpb;

        if (GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("Dx")  << "Solution obtained = " << rDx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rb << std::endl;
        }
        if (this->GetEchoLevel() == 4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] <<  ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rb);
        }
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
    LinearStrategy(const LinearStrategy& Other);

    ///@}

}; /* Class LinearStrategy */

///@}

///@name Type Definitions */
///@{


///@}

} /* namespace Kratos.*/

