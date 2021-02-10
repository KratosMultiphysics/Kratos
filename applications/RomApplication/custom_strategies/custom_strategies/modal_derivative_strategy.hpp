//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Altug Emiroglu, https://github.com/emiroglu
//
//

#if !defined(KRATOS_ROM_MODAL_DERIVATIVE_STRATEGY )
#define  KRATOS_ROM_MODAL_DERIVATIVE_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/timer.h"
#include "utilities/builtin_timer.h"
#include "utilities/atomic_utilities.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "rom_application_variables.h"

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

/** @brief Solving strategy base class
 * @details This is the base class from which we will derive all the strategies (line-search, NR, etc...)
 */

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ModalDerivativeStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

//     typedef std::set<Dof::Pointer,ComparePDof>                                    DofSetType;

    // Base type definition
    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TSchemeType::Pointer TSchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer TBuilderAndSolverPointerType;

    typedef typename BaseType::TDofType TDofType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::DofIteratorType DofIteratorType;

    typedef typename BaseType::DofConstantIteratorType DofConstantIteratorType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ModalDerivativeStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor.
    ModalDerivativeStrategy(
        ModelPart& rModelPart,
        TSchemePointerType pScheme,
        TBuilderAndSolverPointerType pBuilderAndSolver,
        Parameters InputParameters
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart),
        mpScheme(pScheme),
        mpBuilderAndSolver(pBuilderAndSolver)
    {
        KRATOS_TRY

        // Set derivative type: static or dynamic
        const auto& derivative_type = InputParameters["derivative_type"].GetString();
        if (derivative_type == "static")
            mDerivativeTypeFlag = false;
        else if (derivative_type == "dynamic")
            mDerivativeTypeFlag = true;
        else
            KRATOS_ERROR << "\"derivative_type\" can only be \"static\" or \"dynamic\""  << std::endl;

        mMassOrthonormalizeFlag = InputParameters["mass_orthonormalize"].GetBool();
        
        mFixedDofIndex = 0;

        mNumberInitialBasis = rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR].size();
        rModelPart.GetProcessInfo()[DERIVATIVE_INDEX] = mNumberInitialBasis;

        if ( mDerivativeTypeFlag)   // Dynamic derivatives are unsymmetric
            rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR].resize(mNumberInitialBasis*( mNumberInitialBasis + 1 ), true);
        else // Static derivatives are symmetric
            rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR].resize(mNumberInitialBasis + mNumberInitialBasis * ( mNumberInitialBasis + 1 ) / 2, true);
        
        // Set the ordering of DOFs in the input
        const std::vector<std::string> nodal_variable_names = InputParameters["rom_settings"]["nodal_unknowns"].GetStringArray();
        const int number_of_nodal_dofs = nodal_variable_names.size();
        // Setting up mapping: VARIABLE_KEY --> CORRECT_ROW_IN_BASIS
        for(int k=0; k<number_of_nodal_dofs; k++){
            if(KratosComponents<Variable<double>>::Has(nodal_variable_names[k]))
            {
                const auto& var = KratosComponents<Variable<double>>::Get(nodal_variable_names[k]);
                rModelPart.GetProcessInfo()[MAP_PHI][var.Key()] = k;
            }
            else
                KRATOS_ERROR << "variable \""<< nodal_variable_names[k] << "\" not valid" << std::endl;
        }

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build at each solution step)
        this->SetRebuildLevel(1);
        
        mInitializeWasPerformed = false;
        mSolutionStepIsInitialized = false;

        KRATOS_CATCH("")
    }

    // /// Deleted copy constructor.
    // ModalDerivativeStrategy(const ModalDerivativeStrategy& Other) = delete;

    /// Destructor.
    ~ModalDerivativeStrategy() override
    {
        // Clear() controls order of deallocation to avoid invalid memory access
        // in some special cases.
        // warning: BaseType::GetModelPart() may be invalid here.
        this->Clear();
    }
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialization to be performed once before using the strategy.
     */
    void Initialize() override
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        const auto rank = r_model_part.GetCommunicator().MyPID();

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Initialize" << std::endl;

        if (mInitializeWasPerformed == false)
        {
            auto& p_scheme = this->pGetScheme();

            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(r_model_part);

            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(r_model_part);

            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(r_model_part);

            mInitializeWasPerformed = true;
        }

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Initialize" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        auto& p_builder_and_solver = this->pGetBuilderAndSolver();
        // if the preconditioner is saved between solves, it
        // should be cleared here.
        p_builder_and_solver->GetLinearSystemSolver()->Clear();
        if (mpA != nullptr){
            TSparseSpace::Clear(mpA);
            mpA = nullptr;
        }
        if (mpStiffnessMatrix != nullptr){
            TSparseSpace::Clear(mpStiffnessMatrix);
            mpStiffnessMatrix = nullptr;
        }
        if (mpMassMatrix != nullptr){
            TSparseSpace::Clear(mpMassMatrix);
            mpMassMatrix = nullptr;
        }
        if (mpDx != nullptr)
        {
            TSparseSpace::Clear(mpDx);
            mpDx = nullptr;
        }
        if (mpb != nullptr)
        {
            TSparseSpace::Clear(mpb);
            mpb = nullptr;
        }

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        p_builder_and_solver->SetDofSetIsInitializedFlag(false);
        p_builder_and_solver->Clear();
        this->pGetScheme()->Clear();

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mSolutionStepIsInitialized == false)
        {
            //pointers needed in the solution
            auto& p_scheme = this->pGetScheme();
            auto& p_builder_and_solver = this->pGetBuilderAndSolver();
            auto& r_model_part = BaseType::GetModelPart();
            
            const auto rank = r_model_part.GetCommunicator().MyPID();

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                p_builder_and_solver->GetReshapeMatrixFlag() == true)
            {
                // Setting up the list of DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                // Set global equation ids
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(r_model_part);
                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_system_time.ElapsedSeconds() << std::endl;

                // Setting up the Vectors involved to the correct size
                BuiltinTimer system_matrix_resize_time;
                // System matrix necessary in both static and dynamic cases
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb,
                                                                 r_model_part);
                
                if (mDerivativeTypeFlag)
                {   // Dynamic derivatives need both mass and stiffness to be built separately
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpMassMatrix, mpDx, mpb,
                                                                 r_model_part);
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpStiffnessMatrix, mpDx, mpb,
                                                                 r_model_part);
                } 
                else if (!mDerivativeTypeFlag && mMassOrthonormalizeFlag)
                {   // Static derivatives need mass matrix only if mass orthonormalization is necessary
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpMassMatrix, mpDx, mpb,
                                                                 r_model_part);
                }                
                
                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;
            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;

            auto& rA  = *mpA;
            auto& rDx = *mpDx;
            auto& rb  = *mpb;

            //initial operations ... things that are constant over the Solution Step
            p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);

            //initial operations ... things that are constant over the Solution Step
            p_scheme->InitializeSolutionStep(r_model_part, rA, rDx, rb);

            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        const auto rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering FinalizeSolutionStep" << std::endl;

        auto& p_scheme = this->pGetScheme();
        auto& p_builder_and_solver = this->pGetBuilderAndSolver();
        auto& r_model_part = BaseType::GetModelPart();

        auto& rA  = *mpA;
        auto& rDx = *mpDx;
        auto& rb  = *mpb;

        //Finalisation of the solution step
        p_builder_and_solver->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        p_scheme->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        
        //Cleaning memory after the solution
        p_scheme->Clean();

        this->Clear();

        mInitializeWasPerformed = false;
        mSolutionStepIsInitialized = false;

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting FinalizeSolutionStep" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step for static derivatives case
     * 
     */
    bool SolveSolutionStepStaticDerivatives()
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        auto& r_current_process_info = r_model_part.GetProcessInfo();
        auto& p_scheme = this->pGetScheme();
        auto& p_builder_and_solver = this->pGetBuilderAndSolver();
        auto& rA = *mpA;
        auto& rb = *mpb;
        auto& rDx = *mpDx;
        TSystemVectorType basis(p_builder_and_solver->GetDofSet().size());
        const bool master_slave_constraints_defined = r_model_part.MasterSlaveConstraints().size() != 0;

        /*
            BUILD_LEVEL = 1 : Scheme builds M
            BUILD_LEVEL = 2 : Scheme builds K
        */

        // Build system matrices
        // Build mass matrix
        if (mMassOrthonormalizeFlag)
        {
            BuiltinTimer time_mass_matrix;
            auto& rMassMatrix = *mpMassMatrix;        
            r_current_process_info[BUILD_LEVEL] = 1;
            p_builder_and_solver->BuildLHS(p_scheme,r_model_part,rMassMatrix);
            KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Build time MassMatrix: " << time_mass_matrix.ElapsedSeconds() << std::endl;
        }

        // Build stiffness matrix directly into the system matrix
        BuiltinTimer time_system_matrix;
        r_current_process_info[BUILD_LEVEL] = 2;
        p_builder_and_solver->BuildLHS(p_scheme,r_model_part,rA);
        KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Build time SystemMatrix: " << time_system_matrix.ElapsedSeconds() << std::endl;

        // Static derivative LHS is constant. Apply constraints and conditions
        // Master-Slave constraints
        BuiltinTimer time_constraints_and_conditions;
        if (master_slave_constraints_defined)
        {
            p_builder_and_solver->ApplyConstraints(p_scheme, r_model_part, rA, rb);
        }
        // Dirichlet conditions
        p_builder_and_solver->ApplyDirichletConditions_LHS(p_scheme, r_model_part, rA, rDx);
        KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Apply LHS constraints and conditions: " << time_system_matrix.ElapsedSeconds() << std::endl;

        // Get eigenvalues vector
        auto& r_eigenvalues = r_current_process_info[EIGENVALUE_VECTOR];

        for (std::size_t basis_i = 0; basis_i < mNumberInitialBasis; basis_i++)
        {
            // Set the BASIS_I counter for use in the scheme
            r_current_process_info[BASIS_I] = basis_i;

            // Get basis_i
            TSparseSpace::SetToZero(basis);
            this->GetBasis(basis_i, basis);
            // Derivative wrt basis_j (starts from basis_i due to symmetry of static derivatives)
            for (std::size_t basis_j = basis_i; basis_j < mNumberInitialBasis; basis_j++)
            {
                BuiltinTimer time_derivative;
                // Set the BASIS_J counter for using in the scheme
                r_current_process_info[BASIS_J] = basis_j;

                // Reset RHS and solution vector at each step
                TSparseSpace::SetToZero(rb);
                TSparseSpace::SetToZero(rDx);

                // Start building RHS
                BuiltinTimer time_build_rhs;
                // Build RHS: -dK/dp . basis
                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);

                // Compute the derivative of the eigenvalue : basis^T . dK/dp . basis
                const double deigenvalue_i_dbasis_j = -inner_prod(basis, rb);
                r_eigenvalues[r_current_process_info[DERIVATIVE_INDEX]] = deigenvalue_i_dbasis_j;

                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Build time RHS: " << time_build_rhs.ElapsedSeconds() << std::endl;

                // Apply the master-slave constraints to RHS only
                BuiltinTimer time_rhs_constraints_and_conditions;
                if (master_slave_constraints_defined)
                {
                    p_builder_and_solver->ApplyRHSConstraints(p_scheme, r_model_part, rb);
                }

                // Apply Dirichlet conditions
                p_builder_and_solver->ApplyDirichletConditions_RHS(p_scheme, r_model_part, rDx, rb);
                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Apply RHS constraints and conditions: " << time_rhs_constraints_and_conditions.ElapsedSeconds() << std::endl;
                
                // Start Solve
                BuiltinTimer solve;
                // Compute particular solution
                p_builder_and_solver->SystemSolve(rA, rDx, rb);

                // Reconstruct slave DOF solution
                if (master_slave_constraints_defined)
                {
                    p_builder_and_solver->ReconstructSlaveSolution(p_scheme, r_model_part, rA, rDx, rb);
                }
                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "System solve time: " << solve.ElapsedSeconds() << std::endl;

                // Mass orthonormalization
                if (mMassOrthonormalizeFlag)
                {
                    this->MassOrthonormalize(rDx);
                }
                
                // Assign solution to ROM_BASIS
                this->AssignVariables(rDx);

                // Update the derivative index
                r_current_process_info[DERIVATIVE_INDEX] += 1;

                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Compute dBasis_" << basis_i+1 << "_dalpha_" << basis_j+1 << " time: " << time_derivative.ElapsedSeconds() << std::endl;
            }
        }

        KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Eigenvalues and derivatives: " << r_current_process_info[EIGENVALUE_VECTOR] << std::endl;

        return true;

        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step for dynamic derivatives case
     * 
     */
    bool SolveSolutionStepDynamicDerivatives()
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        auto& r_current_process_info = r_model_part.GetProcessInfo();
        auto& p_scheme = this->pGetScheme();
        auto& p_builder_and_solver = this->pGetBuilderAndSolver();
        auto& rA = *mpA;
        auto& rStiffnessMatrix = *mpStiffnessMatrix;
        auto& rMassMatrix = *mpMassMatrix;
        auto& rb = *mpb;
        auto& rDx = *mpDx;
        TSystemVectorType basis(p_builder_and_solver->GetDofSet().size());
        const bool master_slave_constraints_defined = r_model_part.MasterSlaveConstraints().size() != 0;
        bool dynamic_constraint_added = false;

        /*
            BUILD_LEVEL = 1 : Scheme builds M
            BUILD_LEVEL = 2 : Scheme builds K
        */

        // Build system matrices
        // Build mass matrix
        BuiltinTimer time_mass_matrix;
        r_current_process_info[BUILD_LEVEL] = 1;
        p_builder_and_solver->BuildLHS(p_scheme,r_model_part,rMassMatrix);
        KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Build time MassMatrix: " << time_mass_matrix.ElapsedSeconds() << std::endl;

        // Build stiffness matrix separately
        BuiltinTimer time_stiffness_matrix;
        r_current_process_info[BUILD_LEVEL] = 2;
        p_builder_and_solver->BuildLHS(p_scheme,r_model_part,rStiffnessMatrix);
        KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Build time StiffnessMatrix: " << time_stiffness_matrix.ElapsedSeconds() << std::endl;
        
        // Get eigenvalues vector
        LocalSystemVectorType& r_eigenvalues = r_current_process_info[EIGENVALUE_VECTOR];

        // Derivative of basis_i
        for (std::size_t basis_i = 0; basis_i < mNumberInitialBasis; basis_i++)
        {
            // Set the BASIS_I counter for use in the scheme
            r_current_process_info[BASIS_I] = basis_i;

            // Get basis_i
            TSparseSpace::SetToZero(basis);
            this->GetBasis(basis_i, basis);

            // If dynamic derivatives then build system matrix for each eigenvalue
            BuiltinTimer time_system_matrix;
            const auto eigenvalue_i = r_eigenvalues[basis_i];            
            rA = rStiffnessMatrix - (eigenvalue_i * rMassMatrix);
            KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Build time SystemMatrix: " << time_system_matrix.ElapsedSeconds() << std::endl;

            // Dynamic derivative LHS varies for each eigenvalue. Apply constraints and conditions each time
            BuiltinTimer time_constraints_and_conditions;
            // Master-Slave constraints
            if (master_slave_constraints_defined)
            {
                p_builder_and_solver->ApplyConstraints(p_scheme, r_model_part, rA, rb);
            }
            // Add dynamic derivative dirichlet constraint
            this->AddDynamicDerivativeConstraint(basis);
            dynamic_constraint_added = true;
            // Dirichlet conditions
            p_builder_and_solver->ApplyDirichletConditions_LHS(p_scheme, r_model_part, rA, rDx);
            KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Apply LHS constraints and conditions: " << time_system_matrix.ElapsedSeconds() << std::endl;
            
            // Derivative wrt basis_j
            for (std::size_t basis_j = 0; basis_j < mNumberInitialBasis; basis_j++)
            {
                BuiltinTimer time_derivative;

                // Set the BASIS_J counter for using in the scheme
                r_current_process_info[BASIS_J] = basis_j;

                // Reset RHS and solution vector at each step
                TSparseSpace::SetToZero(rb);
                TSparseSpace::SetToZero(rDx);

                // Start building RHS
                BuiltinTimer time_build_rhs;

                // Build RHS partially : -dK/dp . basis
                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);

                // Compute the derivative of the eigenvalue : basis^T . dK/dp . basis
                const double deigenvalue_i_dbasis_j = -inner_prod(basis, rb);
                r_eigenvalues[r_current_process_info[DERIVATIVE_INDEX]] = deigenvalue_i_dbasis_j;

                // Add dynamic part of RHS
                rb += deigenvalue_i_dbasis_j * prod(rMassMatrix, basis);

                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Build time RHS: " << time_build_rhs.ElapsedSeconds() << std::endl;

                // Apply the master-slave constraints to RHS only
                BuiltinTimer time_rhs_constraints_and_conditions;
                if (master_slave_constraints_defined)
                {
                    p_builder_and_solver->ApplyRHSConstraints(p_scheme, r_model_part, rb);
                }

                // Apply Dirichlet conditions
                if (!dynamic_constraint_added)
                {
                    this->AddDynamicDerivativeConstraint(basis);
                    dynamic_constraint_added = true;
                }
                p_builder_and_solver->ApplyDirichletConditions_RHS(p_scheme, r_model_part, rDx, rb);
                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Apply RHS constraints and conditions: " << time_rhs_constraints_and_conditions.ElapsedSeconds() << std::endl;

                // Start Solve
                BuiltinTimer time_solve;
                // Compute particular solution
                p_builder_and_solver->SystemSolve(rA, rDx, rb);

                // Remove dynamic derivative constraint
                this->RemoveDynamicDerivativeConstraint();
                dynamic_constraint_added = false;

                // Reconstruct slave DOF solution
                if (master_slave_constraints_defined)
                {
                    p_builder_and_solver->ReconstructSlaveSolution(p_scheme, r_model_part, rA, rDx, rb);
                }

                // Compute and add null space solution for dynamic derivatives
                this->ComputeAndAddNullSpaceSolution(rDx, basis);

                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "System solve time: " << time_solve.ElapsedSeconds() << std::endl;

                // Mass orthonormalization
                if (mMassOrthonormalizeFlag)
                {
                    this->MassOrthonormalize(rDx);
                }

                // Assign solution to ROM_BASIS
                this->AssignVariables(rDx);

                // Update the derivative index
                r_current_process_info[DERIVATIVE_INDEX] += 1;

                KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) << "Compute dBasis_" << basis_i+1 << "_dalpha_" << basis_j+1 << " time: " << time_derivative.ElapsedSeconds() << std::endl;
            }
            
        }

        KRATOS_INFO_IF("ModalDerivativeStrategy", this->GetEchoLevel() >= 1) 
        << "Eigenvalues and derivatives: " << r_model_part.GetProcessInfo()[EIGENVALUE_VECTOR] << std::endl;

        return true;

        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        if (!mDerivativeTypeFlag)
        {   // Static derivatives
            return this->SolveSolutionStepStaticDerivatives();
        }
        else
        {   // Dynamic derivatives
            return this->SolveSolutionStepDynamicDerivatives();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function applies mass orthonormalization of the given vector
     * @details
     * { 
     * Given vector is mass-orthonormalized w.r.t. the existing and/or previously computed 
     * ROM basis vectors using the modified Gram-Schmidt method
     * } 
     */  
    void MassOrthonormalize(TSystemVectorType& rDx)
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        auto& rMassMatrix = *mpMassMatrix;
        TSystemVectorType basis(rMassMatrix.size1());
        const auto current_basis_index = r_model_part.GetProcessInfo()[DERIVATIVE_INDEX];

        // Mass-Orthogonalization using modified Gram-Schmidt method
        for (std::size_t basis_index = 0; basis_index < current_basis_index; basis_index++){

            // Get previous bases
            TSparseSpace::SetToZero(basis);
            this->GetBasis(basis_index, basis);

            // Apply Mass-Orthogonalization w.r.t. previous bases
            rDx -= (inner_prod(prod(rMassMatrix, basis), rDx) * basis);

        }

        // Mass-Normalization
        rDx /= std::sqrt(inner_prod(prod(rMassMatrix, rDx), rDx));

        KRATOS_CATCH("")        
    }

    /**
     * @brief This function retrieves the basis of given index
     * @details
     * { 
     * } 
     */  
    void GetBasis(std::size_t BasisIndex, TSystemVectorType& rBasis)
    {
        KRATOS_TRY

        const auto& r_model_part = BaseType::GetModelPart();
        const auto& r_current_process_info = r_model_part.GetProcessInfo();
        const auto& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();
        
        block_for_each(
            r_model_part.Nodes(),[&](ModelPart::NodeType& rNode)
            {
                auto& node_dofs = rNode.GetDofs();
                const auto& r_nodal_basis = rNode.GetValue(ROM_BASIS);

                for (const auto& rp_dof : node_dofs)
                {
                    bool is_active = !(r_dof_set.find(*rp_dof) == r_dof_set.end());
                    if (rp_dof->IsFree() && is_active)
                    {
                        const std::size_t dof_index = r_current_process_info[MAP_PHI].at(rp_dof->GetVariable().Key());
                        rBasis[rp_dof->EquationId()] = r_nodal_basis(dof_index,BasisIndex);
                    }
                }
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This function adds the dynamic derivative constraint
     * @details
     * { 
     * The dynamic derivative LHS is singular since A = K - lambda_i*M,
     * thus the dynamic derivative constraint has to be applied.
     * This function adds a Dirichlet constraint on the DOF with the maximum absolute value
     * } 
     */  
    void AddDynamicDerivativeConstraint(TSystemVectorType& rBasis)
    {
        KRATOS_TRY

        auto& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        auto it_max_dof = std::max_element(r_dof_set.begin(), r_dof_set.end(), 
        [&rBasis](const TDofType& rDof1, const TDofType& rDof2){
            return (std::abs(rBasis[rDof1.EquationId()]) < std::abs(rBasis[rDof2.EquationId()]));
        });

        // Fix the found DOF
        it_max_dof->FixDof();
        mFixedDofIndex = it_max_dof->EquationId();

        KRATOS_CATCH("")
    }

    /**
     * @brief This function frees the constrained DOF due to the dynamic derivative constraint
     * @details
     * { 
     * } 
     */  
    void RemoveDynamicDerivativeConstraint()
    {

        auto& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();

        // Free the DOF that is previously fixed
        auto p_dof = r_dof_set.begin()+mFixedDofIndex;
        p_dof->FreeDof();
    }

    /**
     * @brief This function computes and adds the null space solution
     * @details
     * { 
     * } 
     */ 
    void ComputeAndAddNullSpaceSolution(TSystemVectorType& rDx, TSystemVectorType& rBasis)
    {
        KRATOS_TRY

        auto& rMassMatrix = *mpMassMatrix;

        // Component c for the null space solution
        double c = -inner_prod(rDx, prod(rMassMatrix, rBasis));

        rDx += c*rBasis;

        KRATOS_CATCH("")
    }

    /**
     * @brief This function reduces the mass and stiffness matrices using the computed basis
     * @details
     * { 
     * This is implemented purely for debugging purposes
     * } 
     */  
    void ComputeReducedMatrices()
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        std::size_t num_basis = r_model_part.GetProcessInfo()[EIGENVALUE_VECTOR].size();

        Matrix Mred(num_basis,num_basis);
        Matrix Kred(num_basis,num_basis);
        TSystemVectorType basis_i(mpA->size1());
        TSystemVectorType basis_j(mpA->size1());
        TSystemMatrixType& rA = *mpA;
        TSystemMatrixType& rStiffnessMatrix = *mpStiffnessMatrix;
        TSystemMatrixType& rMassMatrix = *mpMassMatrix;

        for (std::size_t i = 0; i < num_basis; i++)
        {
            this->GetBasis(i, basis_i);
            for (std::size_t j = 0; j < num_basis; j++)
            {
                this->GetBasis(j, basis_j);
                Mred(i,j) = inner_prod(basis_i, prod(rMassMatrix, basis_j));
                if (mDerivativeTypeFlag)
                    Kred(i,j) = inner_prod(basis_i, prod(rStiffnessMatrix, basis_j));
                else
                    Kred(i,j) = inner_prod(basis_i, prod(rA, basis_j));
            }
        }

        KRATOS_WATCH(Mred)
        KRATOS_WATCH(Kred)

        KRATOS_CATCH("")
    }

    /**
     * @brief This function assigns the solution vector to ROM_BASIS nodal variable
     * @details
     * { 
     * } 
     */  
    void AssignVariables(TSystemVectorType& rDx)
    {
        KRATOS_TRY
        
        auto& r_model_part = BaseType::GetModelPart();
        
        block_for_each(r_model_part.Nodes(), [&](ModelPart::NodeType& r_node) {

            const auto derivative_index = r_model_part.GetProcessInfo()[DERIVATIVE_INDEX];
            const auto& r_dof_set = this->pGetBuilderAndSolver()->GetDofSet();
            const auto& node_dofs = r_node.GetDofs();

            // fill the ROM_BASIS
            std::size_t dof_ctr = 0;
            for (const auto& rp_node_dof : node_dofs)
            {
                bool is_active = !(r_dof_set.find(*rp_node_dof) == r_dof_set.end());
                if (rp_node_dof->IsFree() && is_active) {
                    r_node.GetValue(ROM_BASIS)(dof_ctr,derivative_index) = rDx(rp_node_dof->EquationId());
                }
                else {
                    r_node.GetValue(ROM_BASIS)(dof_ctr,derivative_index) = 0.0;
                }
                dof_ctr += 1;
            }
        });

        KRATOS_CATCH("")
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
    void SetEchoLevel(const int Level) override
    {
        BaseType::SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    TBuilderAndSolverPointerType& pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    TSchemePointerType& pGetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        const int rank = r_model_part.GetCommunicator().MyPID();

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Entering Check" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(r_model_part);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(r_model_part);

        KRATOS_INFO_IF("ModalDerivativeStrategy", BaseType::GetEchoLevel() > 2 && rank == 0)
            <<  "Exiting Check" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ModalDerivativeStrategy";
    }

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

private:

    ///@}
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    TSchemePointerType mpScheme;

    TBuilderAndSolverPointerType mpBuilderAndSolver;

    bool mDerivativeTypeFlag;

    bool mMassOrthonormalizeFlag;

    std::size_t mNumberInitialBasis;

    std::size_t mFixedDofIndex;

    TSystemMatrixPointerType mpMassMatrix;

    TSystemMatrixPointerType mpStiffnessMatrix;

    TSystemMatrixPointerType mpA;

    TSystemVectorPointerType mpDx;
    
    TSystemVectorPointerType mpb;

    bool mSolutionStepIsInitialized;

    bool mInitializeWasPerformed;

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
    ///@name Un accessible methods
    ///@{

    /** Copy constructor.
     */
    ModalDerivativeStrategy(const ModalDerivativeStrategy& Other);


    ///@}

}; /* Class ModalDerivativeStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_MODAL_DERIVATIVE_STRATEGY  defined */
