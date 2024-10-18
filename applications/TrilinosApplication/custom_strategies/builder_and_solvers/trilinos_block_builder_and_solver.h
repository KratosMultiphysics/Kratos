//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Vicente Mataix
//

#pragma once

// System includes
#include <unordered_set>

// External includes

/* Trilinos includes */
#include <Epetra_FECrsGraph.h>
#include <Epetra_IntVector.h>

// Project includes
#include "trilinos_space.h"
#include "custom_utilities/trilinos_assembling_utilities.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/timer.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"

#if !defined(START_TIMER)
#define START_TIMER(label, rank) \
    if (mrComm.MyPID() == rank)  \
        Timer::Start(label);
#endif
#if !defined(STOP_TIMER)
#define STOP_TIMER(label, rank) \
    if (mrComm.MyPID() == rank) \
        Timer::Stop(label);
#endif

namespace Kratos {

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
 * @class TrilinosBlockBuilderAndSolver
 * @ingroup TrilinosApplication
 * @brief Current class provides an implementation for trilinos builder and
 * solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the
 * residual already contains this information. Calculation of the reactions
 * involves a cost very similar to the calculation of the total residual
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz (MPC)
 * @note Should be TrilinosResidualBasedBlockBuilderAndSolver?
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosBlockBuilderAndSolver
    : public BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the flags
    KRATOS_DEFINE_LOCAL_FLAG( SILENT_WARNINGS );

    /// Definition of the pointer
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosBlockBuilderAndSolver);

    /// Definition of the base class
    using BaseType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// The size_t types
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    /// Definition of the classes from the base class
    using TSchemeType = typename BaseType::TSchemeType;
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// Epetra definitions
    using EpetraCommunicatorType = Epetra_MpiComm;

    /// DoF types definition
    using NodeType = Node;

    /// Defining the sparse matrices and vectors
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    /// Defining the local matrices and vectors
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    /// Definition of the pointer types
    using TSystemMatrixPointerType = typename BaseType::TSystemMatrixPointerType;
    using TSystemVectorPointerType = typename BaseType::TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor (empty)
     */
    explicit TrilinosBlockBuilderAndSolver() = default;

    /**
     * @brief Default constructor.
     */
    explicit TrilinosBlockBuilderAndSolver(EpetraCommunicatorType& rComm,
                                  int GuessRowSize,
                                  typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver),
          mrComm(rComm),
          mGuessRowSize(GuessRowSize)
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit TrilinosBlockBuilderAndSolver(
        EpetraCommunicatorType& rComm,
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver),
            mrComm(rComm)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * Copy constructor
     */
    TrilinosBlockBuilderAndSolver(const TrilinosBlockBuilderAndSolver& rOther) = delete;

    /**
     * Assignment operator
     */
    TrilinosBlockBuilderAndSolver& operator=(const TrilinosBlockBuilderAndSolver& rOther) = delete;

    // TODO: In order to create a Create method, the name of the DataCommunicator that is used needs to be passed in the settings (see DistributedImportModelPartUtility). Then an EpetraComm can be constructed from the MPI_Comm in the DataCommunicator

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to perform the build the system matrix and the residual
     * vector
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void Build(typename TSchemeType::Pointer pScheme,
               ModelPart& rModelPart,
               TSystemMatrixType& rA,
               TSystemVectorType& rb) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_ids_vector;
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        BuiltinTimer build_timer;

        // Assemble all elements
        for (auto it = rModelPart.Elements().ptr_begin(); it < rModelPart.Elements().ptr_end(); it++) {
            // Detect if the element is active or not. If the user did not make any choice the element is active by default
            if ((*it)->IsActive()) {
                // Calculate elemental contribution
                pScheme->CalculateSystemContributions(**it, LHS_Contribution, RHS_Contribution, equation_ids_vector, r_current_process_info);

                // Assemble the elemental contribution
                TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
                TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);
            }
        }

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // Assemble all conditions
        for (auto it = rModelPart.Conditions().ptr_begin(); it < rModelPart.Conditions().ptr_end(); it++) {
            // Detect if the condition is active or not. If the user did not make any choice the condition is active by default
            if ((*it)->IsActive()) {
                // Calculate condition contribution
                pScheme->CalculateSystemContributions(**it, LHS_Contribution, RHS_Contribution, equation_ids_vector, r_current_process_info);

                // Assemble the condition contribution
                TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
                TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);
            }
        }

        // Finalizing the assembly
        rA.GlobalAssemble();
        rb.GlobalAssemble();

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() >= 1) << "Build time: " << build_timer << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation choosen the size of the matrix
     * could be equal to the total number of Dofs or to the number of
     * unrestrained dofs
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     */
    void BuildLHS(typename TSchemeType::Pointer pScheme,
                  ModelPart& rModelPart,
                  TSystemMatrixType& rA) override
    {
        KRATOS_TRY
        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_ids_vector;
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        BuiltinTimer build_timer;

        // Assemble all elements
        for (auto it = rModelPart.Elements().ptr_begin(); it < rModelPart.Elements().ptr_end(); it++) {
            pScheme->CalculateLHSContribution(**it, LHS_Contribution, equation_ids_vector, r_current_process_info);

            // Assemble the elemental contribution
            TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
        }

        LHS_Contribution.resize(0, 0, false);

        // Assemble all conditions
        for (auto it = rModelPart.Conditions().ptr_begin(); it < rModelPart.Conditions().ptr_end(); it++) {
            // calculate elemental contribution
            pScheme->CalculateLHSContribution(**it, LHS_Contribution, equation_ids_vector, r_current_process_info);

            // Assemble the elemental contribution
            TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
        }

        // Finalizing the assembly
        rA.GlobalAssemble();

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() >= 1) << "Build time LHS: " << build_timer << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Build a rectangular matrix of size n*N where "n" is the number of
     * unrestrained degrees of freedom and "N" is the total number of degrees of
     * freedom involved.
     * @details This matrix is obtained by building the total matrix without the
     * lines corresponding to the fixed degrees of freedom (but keeping the
     * columns!!)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS_CompleteOnFreeRows(typename TSchemeType::Pointer pScheme,
                                     ModelPart& rModelPart,
                                     TSystemMatrixType& A) override
    {
        KRATOS_ERROR << "Method BuildLHS_CompleteOnFreeRows not implemented in "
                        "Trilinos Builder And Solver"
                     << std::endl;
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

        // Reference to T
        const TSystemMatrixType& r_T = GetConstraintRelationMatrix();

        // If there are master-slave constraints
        if(TSparseSpace::Size1(r_T) != 0) {
            // Recover solution of the original problem
            TSystemVectorType Dxmodified(rDx);

            // Recover solution of the original problem
            TSparseSpace::Mult(r_T, Dxmodified, rDx);
        }

        // Prints information about the current time
        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void SystemSolveWithPhysics(TSystemMatrixType& rA,
                                TSystemVectorType& rDx,
                                TSystemVectorType& rb,
                                ModelPart& rModelPart
                                )
    {
        KRATOS_TRY

        // If considering MPC
        if (rModelPart.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() > 0) {
            TSystemVectorType Dxmodified(rb);

            // Initialize the vector
            TSparseSpace::SetToZero(Dxmodified);

            InternalSystemSolveWithPhysics(rA, Dxmodified, rb, rModelPart);

            // Reference to T
            const TSystemMatrixType& r_T = GetConstraintRelationMatrix();

            // Recover solution of the original problem
            TSparseSpace::Mult(r_T, Dxmodified, rDx);
        } else {
            InternalSystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void InternalSystemSolveWithPhysics(TSystemMatrixType& rA,
                                TSystemVectorType& rDx,
                                TSystemVectorType& rb,
                                ModelPart& rModelPart)
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(rb) != 0)
            norm_b = TSparseSpace::TwoNorm(rb);
        else
            norm_b = 0.00;

        if (norm_b != 0.00) {
            if (BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded()) {
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(rA, rDx, rb, BaseType::mDofSet, rModelPart);
            }

            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        } else {
            TSparseSpace::SetToZero(rDx);
            KRATOS_WARNING("TrilinosResidualBasedBlockBuilderAndSolver") << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // Prints informations about the current time
        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", (BaseType::GetEchoLevel() > 1)) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is
     * possible to solve just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildAndSolve(typename TSchemeType::Pointer pScheme,
                       ModelPart& rModelPart,
                       TSystemMatrixType& rA,
                       TSystemVectorType& rDx,
                       TSystemVectorType& rb) override
    {
        KRATOS_TRY

        START_TIMER("Build", 0)

        Build(pScheme, rModelPart, rA, rb);

        STOP_TIMER("Build", 0)

        if(rModelPart.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() > 0) {
            const auto timer_constraints = BuiltinTimer();
            START_TIMER("ApplyConstraints", 0)
            ApplyConstraints(pScheme, rModelPart, rA, rb);
            STOP_TIMER("ApplyConstraints", 0)
            KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() >=1) << "Constraints build time: " << timer_constraints << std::endl;
        }

        // Apply dirichlet conditions
        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nBefore the solution of the system"
            << "\nSystem Matrix = " << rA << "\nunknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        START_TIMER("Solve", 0)

        const auto solve_timer = BuiltinTimer();

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() >=1) << "System solve time: " << solve_timer << std::endl;

        STOP_TIMER("Solve", 0)

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nAfter the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;
        KRATOS_CATCH("")
    }

    /**
     * @brief Corresponds to the previous, but the System's matrix is considered
     * already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildRHSAndSolve(typename TSchemeType::Pointer pScheme,
                          ModelPart& rModelPart,
                          TSystemMatrixType& rA,
                          TSystemVectorType& rDx,
                          TSystemVectorType& rb) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, rb);

        if(rModelPart.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() > 0) {
            START_TIMER("ApplyRHSConstraints", 0)
            ApplyRHSConstraints(pScheme, rModelPart, rb);
            STOP_TIMER("ApplyRHSConstraints", 0)
        }

        const auto solve_timer = BuiltinTimer();
        START_TIMER("Solve", 0)

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        STOP_TIMER("Solve", 0)

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() >=1) << "System solve time: " << solve_timer << std::endl;

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the
     * number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(typename TSchemeType::Pointer pScheme,
                  ModelPart& rModelPart,
                  TSystemVectorType& rb) override
    {
        KRATOS_TRY

        START_TIMER("BuildRHS ", 0)
        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_ids_vector;
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Assemble all elements
        for (auto it = rModelPart.Elements().ptr_begin(); it < rModelPart.Elements().ptr_end(); it++) {
            // Calculate elemental Right Hand Side Contribution
            pScheme->CalculateRHSContribution(**it, RHS_Contribution, equation_ids_vector, r_current_process_info);

            // Assemble the elemental contribution
            TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);
        }

        RHS_Contribution.resize(0, false);

        // Assemble all conditions
        for (auto it = rModelPart.Conditions().ptr_begin(); it < rModelPart.Conditions().ptr_end(); it++) {
            // calculate elemental contribution
            pScheme->CalculateRHSContribution(**it, RHS_Contribution, equation_ids_vector, r_current_process_info);

            // Assemble the elemental contribution
            TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);
        }

        // Finalizing the assembly
        rb.GlobalAssemble();

        STOP_TIMER("BuildRHS ", 0)

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking"
     * to each element and condition its Dofs.
     * @details The list of dofs is stores inside the BuilderAndSolver as it is
     * closely connected to the way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        ) override
    {
        KRATOS_TRY

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() > 2) << "Setting up the dofs" << std::endl;

        using DofsVectorType = Element::DofsVectorType;

        // Gets the array of elements from the modeler
        DofsVectorType dof_list, second_dof_list;
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        DofsArrayType temp_dofs_array;
        IndexType guess_num_dofs = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes() * 3;
        temp_dofs_array.reserve(guess_num_dofs);
        BaseType::mDofSet = DofsArrayType();

        // Taking dofs of elements
        auto& r_elements_array = rModelPart.GetCommunicator().LocalMesh().Elements();
        for (auto& r_elem : r_elements_array) {
            pScheme->GetDofList(r_elem, dof_list, r_current_process_info);
            for (auto i_dof = dof_list.begin(); i_dof != dof_list.end(); ++i_dof)
                temp_dofs_array.push_back(*i_dof);
        }

        // Taking dofs of conditions
        auto& r_conditions_array = rModelPart.GetCommunicator().LocalMesh().Conditions();
        for (auto& r_cond : r_conditions_array) {
            pScheme->GetDofList(r_cond, dof_list, r_current_process_info);
            for (auto i_dof = dof_list.begin(); i_dof != dof_list.end(); ++i_dof)
                temp_dofs_array.push_back(*i_dof);
        }

        // Taking dofs of constraints
        auto& r_constraints_array = rModelPart.GetCommunicator().LocalMesh().MasterSlaveConstraints();
        for (auto& r_const : r_constraints_array) {
            r_const.GetDofList(dof_list, second_dof_list, r_current_process_info);
            for (auto i_dof = dof_list.begin(); i_dof != dof_list.end(); ++i_dof)
                temp_dofs_array.push_back(*i_dof);
            for (auto i_dof = second_dof_list.begin(); i_dof != second_dof_list.end(); ++i_dof)
                temp_dofs_array.push_back(*i_dof);
        }

        temp_dofs_array.Unique();
        BaseType::mDofSet = temp_dofs_array;

        // Throws an exception if there are no Degrees of freedom involved in
        // the analysis
        const SizeType number_of_dofs = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(BaseType::mDofSet.size());
        KRATOS_ERROR_IF(number_of_dofs == 0) << "No degrees of freedom!" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have
        // reactions defined This is to be done only in debug mode
        if (BaseType::GetCalculateReactionsFlag()) {
            for (auto dof_iterator = BaseType::mDofSet.begin();
                 dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction())
                    << "Reaction variable not set for the following : " << std::endl
                    << "Node : " << dof_iterator->Id() << std::endl
                    << "Dof : " << (*dof_iterator) << std::endl
                    << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif
        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() > 2) << "Number of degrees of freedom:" << number_of_dofs << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() > 2) << "Finished setting up the dofs" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Organizes the DoF set in order to speed up the building phase
     *          Sets equation id for degrees of freedom
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(ModelPart& rModelPart) override
    {
        int free_size = 0;
        auto& r_comm = rModelPart.GetCommunicator();
        const auto& r_data_comm = r_comm.GetDataCommunicator();
        const int current_rank = r_data_comm.Rank();
        const int world_size = r_data_comm.Size();

        // Calculating number of fixed and free dofs
        for (const auto& r_dof : BaseType::mDofSet) {
            if (r_dof.GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                free_size++;
            }
        }

        // Calculating the total size and required offset
        int free_offset;
        int global_size;

        // The corresponding offset by the sum of the sizes in thread with
        // inferior current_rank
        free_offset = r_data_comm.ScanSum(free_size);

        // The total size by the sum of all size in all threads
        global_size = r_data_comm.SumAll(free_size);

        // Finding the offset for the beginning of the partition
        free_offset -= free_size;

        // Now setting the equation id with .
        for (auto& r_dof : BaseType::mDofSet) {
            if (r_dof.GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                r_dof.SetEquationId(free_offset++);
            }
        }

        BaseType::mEquationSystemSize = global_size;
        mLocalSystemSize = free_size;
        KRATOS_INFO_IF_ALL_RANKS("TrilinosBlockBuilderAndSolver", BaseType::GetEchoLevel() > 1)
            << "\n    BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize
            << "\n    mLocalSystemSize = " << mLocalSystemSize
            << "\n    free_offset = " << free_offset << std::endl;

        // by Riccardo ... it may be wrong!
        mFirstMyId = free_offset - mLocalSystemSize;

        // For MPC we fill mFirstMyIds
        if (r_comm.GlobalNumberOfMasterSlaveConstraints() > 0) {
            mFirstMyIds.resize(world_size);
            std::vector<int> send_first_ids(1, mFirstMyId);
            r_data_comm.AllGather(send_first_ids, mFirstMyIds);
        }

        // Synchronize DoFs
        r_comm.SynchronizeDofs();
    }

    /**
     * @brief Resizes the system matrix and the vector according to the number
     * of dos in the current rModelPart. This function also decides on the
     * sparsity pattern and the graph of the Trilinos csr matrix
     * @param pScheme The integration scheme considered
     * @param rpA The LHS matrix
     * @param rpDx The Unknowns vector
     * @param rpd The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& rpA,
        TSystemVectorPointerType& rpDx,
        TSystemVectorPointerType& rpb,
        ModelPart& rModelPart
        ) override
    {
        KRATOS_TRY

        // Resizing the system vectors and matrix
        if (rpA == nullptr || TSparseSpace::Size1(*rpA) == 0 || BaseType::GetReshapeMatrixFlag()) { // If the matrix is not initialized
            ConstructMatrixStructure(pScheme, rpA, rpDx, rpb, rModelPart);
        } else if (BaseType::mpReactionsVector == nullptr && this->mCalculateReactionsFlag) {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(rpDx->Map()));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        } else {
            if (TSparseSpace::Size1(*rpA) == 0 ||
                TSparseSpace::Size1(*rpA) != BaseType::mEquationSystemSize ||
                TSparseSpace::Size2(*rpA) != BaseType::mEquationSystemSize) {
                KRATOS_ERROR << "It should not come here resizing is not allowed this way!!!!!!!! ... " << std::endl;
            }
        }

        ConstructMasterSlaveConstraintsStructure(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief It computes the reactions of the system
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unknowns
     * @param rb The RHS vector of the system of equations
     */
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        TSparseSpace::SetToZero(rb);

        // Refresh RHS to have the correct reactions
        BuildRHS(pScheme, rModelPart, rb);

        // Initialize the Epetra importer
        // TODO: this part of the code has been pasted until a better solution
        // is found
        int system_size = TSparseSpace::Size(rb);
        int number_of_dofs = BaseType::mDofSet.size();
        std::vector<int> index_array(number_of_dofs);

        // Filling the array with the global ids
        int counter = 0;
        int id = 0;
        for (const auto& dof : BaseType::mDofSet) {
            id = dof.EquationId();
            if (id < system_size) {
                index_array[counter++] = id;
            }
        }

        std::sort(index_array.begin(), index_array.end());
        std::vector<int>::iterator NewEnd = std::unique(index_array.begin(), index_array.end());
        index_array.resize(NewEnd - index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        rb.Comm().SumAll(&tot_update_dofs, &check_size, 1);
        KRATOS_ERROR_IF(check_size < system_size)
            << "Dof count is not correct. There are less dofs than expected.\n"
            << "Expected number of active dofs = " << system_size
            << " dofs found = " << check_size << std::endl;

        // Defining a map as needed
        Epetra_Map dof_update_map(-1, index_array.size(),
                                  &(*(index_array.begin())), 0, rb.Comm());

        // Defining the importer class
        Kratos::shared_ptr<Epetra_Import> pDofImporter = Kratos::make_shared<Epetra_Import>(dof_update_map, rb.Map());

        // Defining a temporary vector to gather all of the values needed
        Epetra_Vector temp_RHS(pDofImporter->TargetMap());

        // Importing in the new temp_RHS vector the values
        int ierr = temp_RHS.Import(rb, *pDofImporter, Insert);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found - error code: " << ierr << std::endl;

        double* temp_RHS_values; // DO NOT make delete of this one!!
        temp_RHS.ExtractView(&temp_RHS_values);

        rb.Comm().Barrier();

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        // Store the RHS values in the reaction variable
        // NOTE: dofs are assumed to be numbered consecutively in the
        // BlockBuilderAndSolver
        for (int k = 0; k < ndofs; k++) {
            auto dof_iterator = BaseType::mDofSet.begin() + k;

            const int i = (dof_iterator)->EquationId();
            // (dof_iterator)->GetSolutionStepReactionValue() = -(*b[i]);
            const double react_val = temp_RHS[pDofImporter->TargetMap().LID(i)];
            (dof_iterator->GetSolutionStepReactionValue()) = -react_val;
        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy
     * or completely unexpensive depending on the implementation choosen and on
     * how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation
     * the user should refer to the particular Builder And Solver choosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // loop over all dofs to find the fixed ones
        std::vector<int> global_ids(BaseType::mDofSet.size());
        std::vector<int> is_dirichlet(BaseType::mDofSet.size());

        IndexType i = 0;
        for (const auto& dof : BaseType::mDofSet) {
            const int global_id = dof.EquationId();
            global_ids[i] = global_id;
            is_dirichlet[i] = dof.IsFixed();
            ++i;
        }

        // here we construct and fill a vector "fixed local" which cont
        Epetra_Map localmap(-1, global_ids.size(), global_ids.data(), 0, rA.Comm());
        Epetra_IntVector fixed_local(Copy, localmap, is_dirichlet.data());

        Epetra_Import dirichlet_importer(rA.ColMap(), fixed_local.Map());

        // defining a temporary vector to gather all of the values needed
        Epetra_IntVector fixed(rA.ColMap());

        // Detect if there is a line of all zeros and set the diagonal to a certain number (1 if not scale, some norms values otherwise) if this happens
        const auto& r_process_info = rModelPart.GetProcessInfo();
        mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(r_process_info, rA, rb, mScalingDiagonal);

        // Importing in the new temp vector the values
        int ierr = fixed.Import(fixed_local, dirichlet_importer, Insert);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;

        for (int i = 0; i < rA.NumMyRows(); i++) {
            int numEntries; // Number of non-zero entries
            double* vals;   // Row non-zero values
            int* cols;      // Column indices of row non-zero values
            rA.ExtractMyRowView(i, numEntries, vals, cols);

            const int row_gid = rA.RowMap().GID(i);
            const int row_lid = localmap.LID(row_gid);

            if (fixed_local[row_lid] == 0) { // Not a dirichlet row
                for (int j = 0; j < numEntries; j++) {
                    if (fixed[cols[j]] == true)
                        vals[j] = 0.0;
                }
            } else { // This IS a dirichlet row
                // Set to zero the rhs
                rb[0][i] = 0.0; // note that the index of i is expected to be
                                // coherent with the rows of A

                // Set to zero the whole row
                for (int j = 0; j < numEntries; j++) {
                    int col_gid = rA.ColMap().GID(cols[j]);
                    if (col_gid != row_gid)
                        vals[j] = 0.0;
                }
            }
        }

        KRATOS_CATCH("");
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

        if (rModelPart.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() > 0) {
            BuildMasterSlaveConstraints(rModelPart);

            // Reference to T
            const TSystemMatrixType& r_T = GetConstraintRelationMatrix();

            // Compute T' b
            const TSystemVectorType copy_b(rb);
            TSparseSpace::TransposeMult(r_T, copy_b, rb);

            // Apply diagonal values on slaves
            IndexPartition<std::size_t>(mSlaveIds.size()).for_each([&](std::size_t Index){
                const IndexType slave_equation_id = mSlaveIds[Index];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(rb, slave_equation_id, 0.0);
                }
            });

            // Global assembly
            rb.GlobalAssemble();
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

        if (rModelPart.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() > 0) {
            BuildMasterSlaveConstraints(rModelPart);

            // Reference to T
            const TSystemMatrixType& r_T = GetConstraintRelationMatrix();

            // Compute T' A T
            const TSystemMatrixType copy_A(rA);
            TSparseSpace::BtDBProductOperation(rA, copy_A, r_T, true, true);

            // Compute T' b
            const TSystemVectorType copy_b(rb);
            TSparseSpace::TransposeMult(r_T, copy_b, rb);

            // Compute the scale factor value
            const auto& r_process_info = rModelPart.GetProcessInfo();
            mScaleFactor = TSparseSpace::GetScaleNorm(r_process_info, rA, mScalingDiagonal);

            // Apply diagonal values on slaves
            IndexPartition<std::size_t>(mSlaveIds.size()).for_each([&](std::size_t Index){
                const IndexType slave_equation_id = mSlaveIds[Index];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(rA, slave_equation_id, slave_equation_id, mScaleFactor);
                    TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(rb, slave_equation_id, 0.0);
                }
            });

            // Global assembly
            rb.GlobalAssemble();
            rA.GlobalAssemble();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();

        mSlaveIds.clear();
        mMasterIds.clear();
        mInactiveSlaveDofs.clear();
        TSparseSpace::Clear(mpT);
        TSparseSpace::Clear(mpConstantVector);
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model part of the problem to solve
     * @return 0 all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "trilinos_block_builder_and_solver",
            "guess_row_size"                       : 45,
            "block_builder"                        : true,
            "diagonal_values_for_dirichlet_dofs"   : "use_max_diagonal",
            "silent_warnings"                      : false
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
        return "trilinos_block_builder_and_solver";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns constraint relation (T) matrix
     * @return The constraint relation (T) matrix
     */
    typename TSparseSpace::MatrixType& GetConstraintRelationMatrix() override
    {
        auto& r_T = *mpT;
        return r_T;
    }

    /**
     * @brief This method returns constraint constant vector
     * @return The constraint constant vector
     */
    typename TSparseSpace::VectorType& GetConstraintConstantVector() override
    {
        auto& r_constant_vector = *mpConstantVector;
        return r_constant_vector;
    }

    /**
    * @brief Retrieves the current scale factor.
    * This function returns the current scale factor value.
    * @return Returns the current scale factor.
    */
    double GetScaleFactor()
    {
        return mScaleFactor;
    }

    /**
    * @brief Sets the scale factor.
    * This function sets a new value for the scale factor.
    * @param ScaleFactor The new value for the scale factor.
    */
    void SetScaleFactor(const double ScaleFactor)
    {
        mScaleFactor = ScaleFactor;
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
        return "TrilinosBlockBuilderAndSolver";
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
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    /* Base variables */
    EpetraCommunicatorType& mrComm;                 /// The MPI communicator
    int mGuessRowSize;                              /// The guess row size
    IndexType mLocalSystemSize;                     /// The local system size
    int mFirstMyId;                                 /// Auxiliary Id (the first row of the local system)
    int mLastMyId;                                  /// Auxiliary Id (the last row of the local system) // TODO: This can be removed as can be deduced from mLocalSystemSize
    Kratos::shared_ptr<Epetra_Map> mpMap = nullptr; /// The map considered for the different vectors and matrices
    std::vector<int> mFirstMyIds;                   /// The ids corresponding to each partition (only used with MPC)

    /* MPC variables */
    TSystemMatrixPointerType mpT =  nullptr;              /// This is matrix containing the global relation for the constraints
    TSystemVectorPointerType mpConstantVector =  nullptr; /// This is vector containing the rigid movement of the constraint
    std::vector<IndexType> mSlaveIds;                     /// The equation ids of the slaves
    std::vector<IndexType> mMasterIds;                    /// The equation ids of the master
    std::unordered_set<IndexType> mInactiveSlaveDofs;     /// The set containing the inactive slave dofs
    double mScaleFactor = 1.0;                            /// The manually set scale factor

    /* Flags */
    SCALING_DIAGONAL mScalingDiagonal = SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL; /// We identify the scaling considered for the dirichlet dofs
    Flags mOptions;                                                              /// Some flags used internally

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Constructs the master-slave constraints structure for the given model part.
     * @param rModelPart the model part
     * @throws ErrorType describes an Epetra failure in Graph.InsertGlobalIndices
     */
    virtual void ConstructMasterSlaveConstraintsStructure(ModelPart& rModelPart)
    {
        if (rModelPart.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() > 0) {
            START_TIMER("ConstraintsRelationMatrixStructure", 0)
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // Generate indices database
            const IndexType number_of_local_rows = mLocalSystemSize;

            // Generate map - use the "temp" array here
            const int temp_size = (number_of_local_rows < 1000) ? 1000 : number_of_local_rows;
            std::vector<int> temp_primary(temp_size, 0);
            std::vector<int> temp_secondary(temp_size, 0);
            for (IndexType i = 0; i != number_of_local_rows; i++) {
                temp_primary[i] = mFirstMyId + i;
            }
            Epetra_Map& r_map = GetEpetraMap();
            std::fill(temp_primary.begin(), temp_primary.begin() + number_of_local_rows, 0);

            // The T graph
            Epetra_FECrsGraph Tgraph(Copy, r_map, mGuessRowSize);

            // Adding diagonal values
            int ierr;
            int index_diagonal[1] = {0};
            for (IndexType i = 0; i < number_of_local_rows; ++i) {
                index_diagonal[0] = mFirstMyId + i;
                ierr = Tgraph.InsertGlobalIndices(1, index_diagonal, 1, index_diagonal);
                KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
            }

            // Vector containing indices belonging to slave DoFs, not used for graph, but for master/slave index identification
            std::unordered_set<std::size_t> indices;

            // TODO: Check if these should be local constraints
            auto& r_constraints_array = rModelPart.MasterSlaveConstraints();

            // Assemble all constraints
            Element::EquationIdVectorType slave_equation_ids_vector, master_equation_ids_vector;
            for (auto& r_const : r_constraints_array) {
                r_const.EquationIdVector(slave_equation_ids_vector, master_equation_ids_vector, r_current_process_info);

                // Filling the list of active global indices (non fixed)
                IndexType num_active_slave_indices = 0;
                for (auto& r_slave_id : slave_equation_ids_vector) {
                    temp_primary[num_active_slave_indices] = r_slave_id;
                    ++num_active_slave_indices;
                }
                IndexType num_active_master_indices = 0;
                for (auto& r_master_id : master_equation_ids_vector) {
                    temp_secondary[num_active_master_indices] = r_master_id;
                    ++num_active_master_indices;
                }

                // Adding cross master-slave dofs
                if (num_active_slave_indices > 0 && num_active_master_indices > 0) {
                    int slave_index[1] = {0};
                    for (IndexType i = 0; i < num_active_slave_indices; ++i) {
                        slave_index[0] = temp_primary[i];
                        indices.insert(temp_primary[i]);
                        ierr = Tgraph.InsertGlobalIndices(1, slave_index, num_active_master_indices, temp_secondary.data());
                        KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                    }
                }
                std::fill(temp_primary.begin(), temp_primary.begin() + num_active_slave_indices, 0);
                std::fill(temp_secondary.begin(), temp_secondary.begin() + num_active_master_indices, 0);
            }

            // Finalizing graph construction
            ierr = Tgraph.GlobalAssemble();
            KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Epetra_FECrsGraph.GlobalAssemble. Error code: " << ierr << std::endl;
            ierr = Tgraph.FillComplete();
            KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Epetra_FECrsGraph.FillComplete. Error code: " << ierr << std::endl;
            ierr = Tgraph.OptimizeStorage();
            KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Epetra_FECrsGraph.OptimizeStorage. Error code: " << ierr << std::endl;

            // Generate a new matrix pointer according to this non-zero values
            TSystemMatrixPointerType p_new_T = TSystemMatrixPointerType(new TSystemMatrixType(Copy, Tgraph));

            // Swap matrix
            mpT.swap(p_new_T);

            // Generate the constant vector equivalent
            TSystemVectorPointerType p_new_constant_vector = TSystemVectorPointerType(new TSystemVectorType(r_map));
            mpConstantVector.swap(p_new_constant_vector);

            /* Fill ids for master/slave */

            // First clear the ones considered
            mSlaveIds.clear();
            mMasterIds.clear();

            // Now prepare the auxiliary ones
            auto& r_comm = rModelPart.GetCommunicator();
            const auto& r_data_comm = r_comm.GetDataCommunicator();
            const int current_rank = r_data_comm.Rank();
            const int world_size = r_data_comm.Size();

            // Auxiliary slave ids
            std::vector<std::unordered_set<IndexType>> auxiliary_slave_ids(world_size);
            for (auto index : indices) {
                const IndexType rank = DeterminePartitionIndex(index);
                auxiliary_slave_ids[rank].insert(index);
            }
            const int tag_sync_slave_id = 0;
            for (int i_rank = 0; i_rank < world_size; ++i_rank) {
                if (i_rank != current_rank) {
                    std::vector<IndexType> receive_slave_ids_vector;
                    r_data_comm.Recv(receive_slave_ids_vector, i_rank, tag_sync_slave_id);
                    auxiliary_slave_ids[i_rank].insert(receive_slave_ids_vector.begin(), receive_slave_ids_vector.end());
                } else {
                    for (int j_rank = 0; j_rank < world_size; ++j_rank) {
                        if (j_rank != current_rank) {
                            const auto& r_slave_ids = auxiliary_slave_ids[j_rank];
                            std::vector<IndexType> send_slave_ids_vector(r_slave_ids.begin(), r_slave_ids.end());
                            r_data_comm.Send(send_slave_ids_vector, j_rank, tag_sync_slave_id);
                        }
                    }
                }
            }
            mSlaveIds = std::vector<IndexType>(auxiliary_slave_ids[current_rank].begin(), auxiliary_slave_ids[current_rank].end());

            // Master DoFs are complementary
            std::unordered_set<IndexType> temp_master_ids;
            for (IndexType i = 0; i < number_of_local_rows; ++i) {
                temp_master_ids.insert(mFirstMyId + i);
            }

            // Remove ids
            for (auto id_slave_vector : auxiliary_slave_ids) {
                for (auto id_slave : id_slave_vector) {
                    temp_master_ids.erase(id_slave);
                }
            }
            mMasterIds = std::vector<IndexType>(temp_master_ids.begin(), temp_master_ids.end());

            STOP_TIMER("ConstraintsRelationMatrixStructure", 0)
        }
    }

    /**
     * @brief Builds the master-slave constraints for the given model part.
     * @param rModelPart the model part to build the constraints for
     */
    virtual void BuildMasterSlaveConstraints(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // Reference of the matrix and vector
        auto& r_T = GetConstraintRelationMatrix();
        auto& r_constant_vector = GetConstraintConstantVector();

        TSparseSpace::SetToZero(r_T);
        TSparseSpace::SetToZero(r_constant_vector);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Contributions to the system
        Matrix transformation_matrix = LocalSystemMatrixType(0, 0);
        Vector constant_vector = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType slave_equation_ids, master_equation_ids;

        // Now prepare the auxiliary ones
        auto& r_comm = rModelPart.GetCommunicator();
        const auto& r_data_comm = r_comm.GetDataCommunicator();
        const int current_rank = r_data_comm.Rank();
        const int world_size = r_data_comm.Size();

        // Auxiliary inactive slave ids
        std::vector<std::unordered_set<IndexType>> auxiliary_inactive_slave_ids(world_size);

        // We clear the set
        mInactiveSlaveDofs.clear();
        std::size_t num_inactive_slave_dofs_other_partitions = 0;

        // Iterate over the constraints
        for (auto& r_const : rModelPart.MasterSlaveConstraints()) {
            r_const.EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
            // Detect if the constraint is active or not. If the user did not make any choice the constraint. It is active by default
            if (r_const.IsActive()) {
                r_const.CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);

                TrilinosAssemblingUtilities::AssembleRelationMatrixT(r_T, transformation_matrix, slave_equation_ids, master_equation_ids);
                TrilinosAssemblingUtilities::AssembleConstantVector(r_constant_vector, constant_vector, slave_equation_ids);
            } else { // Taking into account inactive constraints
                // Save the auxiliary ids of the the slave inactive DoFs
                for (auto slave_id : slave_equation_ids) {
                    const int index_rank = DeterminePartitionIndex(slave_id);
                    if (index_rank == current_rank) {
                        mInactiveSlaveDofs.insert(slave_id);
                    } else {
                        auxiliary_inactive_slave_ids[index_rank].insert(slave_id);
                        ++num_inactive_slave_dofs_other_partitions;
                    }
                }
            }
        }

        // Compute total number of inactive slave dofs in other partitions
        num_inactive_slave_dofs_other_partitions = r_data_comm.SumAll(num_inactive_slave_dofs_other_partitions);

        // Now we pass the info between partitions
        if (num_inactive_slave_dofs_other_partitions > 0) {
            const int tag_sync_inactive_slave_id = 0;
            for (int i_rank = 0; i_rank < world_size; ++i_rank) {
                if (i_rank != current_rank) {
                    std::vector<IndexType> receive_inactive_slave_ids_vector;
                    r_data_comm.Recv(receive_inactive_slave_ids_vector, i_rank, tag_sync_inactive_slave_id);
                    mInactiveSlaveDofs.insert(receive_inactive_slave_ids_vector.begin(), receive_inactive_slave_ids_vector.end());
                } else {
                    for (int j_rank = 0; j_rank < world_size; ++j_rank) {
                        if (j_rank != current_rank) {
                            const auto& r_inactive_slave_ids = auxiliary_inactive_slave_ids[j_rank];
                            std::vector<IndexType> send_inactive_slave_ids_vector(r_inactive_slave_ids.begin(), r_inactive_slave_ids.end());
                            r_data_comm.Send(send_inactive_slave_ids_vector, j_rank, tag_sync_inactive_slave_id);
                        }
                    }
                }
            }
        }

        // Setting the master dofs into the T and C system
        for (auto eq_id : mMasterIds) {
            TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(r_constant_vector, eq_id, 0.0);
            TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(r_T, eq_id, eq_id, 1.0);
        }

        // Setting inactive slave dofs in the T and C system
        for (auto eq_id : mInactiveSlaveDofs) {
            TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(r_constant_vector, eq_id, 0.0);
            TrilinosAssemblingUtilities::SetGlobalValueWithoutGlobalAssembly(r_T, eq_id, eq_id, 1.0);
        }

        // Finalizing the assembly
        r_T.GlobalAssemble();
        r_constant_vector.GlobalAssemble();

        KRATOS_CATCH("")
    }

    /**
     * @brief Constructs the matrix structure for the given problem.
     * @param pScheme pointer to the scheme object
     * @param rpA reference to the system matrix pointer
     * @param rpDx reference to the system vector pointer for the solution vector
     * @param rpb reference to the system vector pointer for the right-hand side vector
     * @param rModelPart reference to the model part object
     * @throws Runtime error if there is an error in constructing the matrix structure
     */
    virtual void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& rpA,
        TSystemVectorPointerType& rpDx,
        TSystemVectorPointerType& rpb,
        ModelPart& rModelPart
        )
    {
        // Filling with zero the matrix (creating the structure)
        START_TIMER("MatrixStructure", 0)

        // TODO: Check if these should be local elements, conditions and constraints
        auto& r_elements_array = rModelPart.Elements();
        auto& r_conditions_array = rModelPart.Conditions();
        auto& r_constraints_array = rModelPart.MasterSlaveConstraints();

        // Number of local dofs
        const IndexType number_of_local_rows = mLocalSystemSize;

        // Generate map - use the "temp" array here
        const int temp_size = (number_of_local_rows < 1000) ? 1000 : number_of_local_rows;
        std::vector<int> temp_primary(temp_size, 0);
        std::vector<int> temp_secondary(temp_size, 0);
        for (IndexType i = 0; i != number_of_local_rows; i++) {
            temp_primary[i] = mFirstMyId + i;
        }
        Epetra_Map& r_map = GetEpetraMap();
        std::fill(temp_primary.begin(), temp_primary.begin() + number_of_local_rows, 0);

        // Create and fill the graph of the matrix --> the temp array is
        // reused here with a different meaning
        Epetra_FECrsGraph Agraph(Copy, r_map, mGuessRowSize);
        Element::EquationIdVectorType equation_ids_vector;
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Trilinos error int definition
        int ierr;

        // Assemble all elements
        for (auto& r_elem : r_elements_array) {
            pScheme->EquationId(r_elem, equation_ids_vector, r_current_process_info);

            // Filling the list of active global indices (non fixed)
            IndexType num_active_indices = 0;
            for (auto& r_id : equation_ids_vector) {
                temp_primary[num_active_indices] = r_id;
                ++num_active_indices;
            }

            if (num_active_indices != 0) {
                ierr = Agraph.InsertGlobalIndices(num_active_indices, temp_primary.data(), num_active_indices, temp_primary.data());
                KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
            }
            std::fill(temp_primary.begin(), temp_primary.begin() + num_active_indices, 0);
        }

        // Assemble all conditions
        for (auto& r_cond : r_conditions_array) {
            pScheme->EquationId(r_cond, equation_ids_vector, r_current_process_info);

            // Filling the list of active global indices (non fixed)
            IndexType num_active_indices = 0;
            for (auto& r_id : equation_ids_vector) {
                temp_primary[num_active_indices] = r_id;
                ++num_active_indices;
            }

            if (num_active_indices != 0) {
                ierr = Agraph.InsertGlobalIndices(num_active_indices, temp_primary.data(), num_active_indices, temp_primary.data());
                KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
            }
            std::fill(temp_primary.begin(), temp_primary.begin() + num_active_indices, 0);
        }

        // Assemble all constraints
        Element::EquationIdVectorType slave_equation_ids_vector, master_equation_ids_vector;
        for (auto& r_const : r_constraints_array) {
            r_const.EquationIdVector(slave_equation_ids_vector, master_equation_ids_vector, r_current_process_info);

            // Filling the list of active global indices (non fixed)
            IndexType num_active_slave_indices = 0;
            for (auto& r_slave_id : slave_equation_ids_vector) {
                temp_primary[num_active_slave_indices] = r_slave_id;
                ++num_active_slave_indices;
            }
            IndexType num_active_master_indices = 0;
            for (auto& r_master_id : master_equation_ids_vector) {
                temp_secondary[num_active_master_indices] = r_master_id;
                ++num_active_master_indices;
            }

            // First adding the pure slave dofs
            if (num_active_slave_indices > 0) {
                int index[1] = {0};
                for (IndexType i = 0; i < num_active_slave_indices; ++i) {
                    index[0] = temp_primary[i];
                    ierr = Agraph.InsertGlobalIndices(1, index, 1, index);
                    KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
                // Now adding cross master-slave dofs
                if (num_active_master_indices > 0) {
                    ierr = Agraph.InsertGlobalIndices(num_active_slave_indices, temp_primary.data(), num_active_master_indices, temp_secondary.data());
                    KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
            }
            // Second adding pure master dofs
            if (num_active_master_indices > 0) {
                int index[1] = {0};
                for (IndexType i = 0; i < num_active_master_indices; ++i) {
                    index[0] = temp_secondary[i];
                    ierr = Agraph.InsertGlobalIndices(1, index, 1, index);
                    KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
            }
            std::fill(temp_primary.begin(), temp_primary.begin() + num_active_slave_indices, 0);
            std::fill(temp_secondary.begin(), temp_secondary.begin() + num_active_master_indices, 0);
        }

        // Finalizing graph construction
        ierr = Agraph.GlobalAssemble();
        KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Epetra_FECrsGraph.GlobalAssemble. Error code: " << ierr << std::endl;
        ierr = Agraph.FillComplete();
        KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Epetra_FECrsGraph.FillComplete. Error code: " << ierr << std::endl;
        ierr = Agraph.OptimizeStorage();
        KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Epetra_FECrsGraph.OptimizeStorage. Error code: " << ierr << std::endl;

        // Generate a new matrix pointer according to this graph
        TSystemMatrixPointerType p_new_A = TSystemMatrixPointerType(new TSystemMatrixType(Copy, Agraph));
        rpA.swap(p_new_A);

        // Generate new vector pointers according to the given map
        if (rpb == nullptr || TSparseSpace::Size(*rpb) != BaseType::mEquationSystemSize) {
            TSystemVectorPointerType p_new_b = TSystemVectorPointerType(new TSystemVectorType(r_map));
            rpb.swap(p_new_b);
        }
        if (rpDx == nullptr || TSparseSpace::Size(*rpDx) != BaseType::mEquationSystemSize) {
            TSystemVectorPointerType p_new_Dx = TSystemVectorPointerType(new TSystemVectorType(r_map));
            rpDx.swap(p_new_Dx);
        }
        // If the pointer is not initialized initialize it to an empty matrix
        if (BaseType::mpReactionsVector == nullptr) {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(r_map));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        STOP_TIMER("MatrixStructure", 0)
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Get guess row size
        mGuessRowSize = ThisParameters["guess_row_size"].GetInt();

        // Setting flags
        const std::string& r_diagonal_values_for_dirichlet_dofs = ThisParameters["diagonal_values_for_dirichlet_dofs"].GetString();

        const std::set<std::string> available_options_for_diagonal = {"no_scaling","use_max_diagonal","use_diagonal_norm","defined_in_process_info"};

        if (available_options_for_diagonal.find(r_diagonal_values_for_dirichlet_dofs) == available_options_for_diagonal.end()) {
            std::stringstream msg;
            msg << "Currently prescribed diagonal values for dirichlet dofs : " << r_diagonal_values_for_dirichlet_dofs << "\n";
            msg << "Admissible values for the diagonal scaling are : 'no_scaling', 'use_max_diagonal', 'use_diagonal_norm', or 'defined_in_process_info'" << "\n";
            KRATOS_ERROR << msg.str() << std::endl;
        }

        // The first option will not consider any scaling (the diagonal values will be replaced with 1)
        if (r_diagonal_values_for_dirichlet_dofs == "no_scaling") {
            mScalingDiagonal = SCALING_DIAGONAL::NO_SCALING;
        } else if (r_diagonal_values_for_dirichlet_dofs == "use_max_diagonal") {
            mScalingDiagonal = SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL;
        } else if (r_diagonal_values_for_dirichlet_dofs == "use_diagonal_norm") { // On this case the norm of the diagonal will be considered
            mScalingDiagonal = SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL;
        } else { // Otherwise we will assume we impose a numerical value
            mScalingDiagonal = SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL;
        }
        mOptions.Set(SILENT_WARNINGS, ThisParameters["silent_warnings"].GetBool());
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Generates the EpetraMap used for the vectors and matrices
     * @return Returns the Epetra_Map considered for the graphs
     */
    Epetra_Map& GetEpetraMap()
    {
        if (mpMap == nullptr) {
            // Generate map - use the "temp" array here
            const int temp_size = (mLocalSystemSize < 1000) ? 1000 : mLocalSystemSize;
            std::vector<int> temp_primary(temp_size, 0);
            for (IndexType i = 0; i != mLocalSystemSize; i++) {
                temp_primary[i] = mFirstMyId + i;
            }
            mpMap = Kratos::make_shared<Epetra_Map>(-1, mLocalSystemSize, temp_primary.data(), 0, mrComm);
        }

        return *mpMap;
    }

    /**
     * @brief Determine in which partition the index belongs
     * @param Index The index where determine the partition
     * @return The partition where the index belongs
     */
    IndexType DeterminePartitionIndex(const IndexType Index)
    {
        IndexType index;
        for (index = 0; index < mFirstMyIds.size() - 1; ++index) {
            if (Index < static_cast<unsigned int>(mFirstMyIds[index + 1])) {
                break;
            }
        }
        return index;
    }

    void AssembleLHS_CompleteOnFreeRows(TSystemMatrixType& rA,
                                        LocalSystemMatrixType& rLHS_Contribution,
                                        Element::EquationIdVectorType& rEquationId)
    {
        KRATOS_ERROR << "This method is not implemented for Trilinos";
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; /* Class TrilinosBlockBuilderAndSolver */

///@}

///@name Type Definitions
///@{

// Here one should use the KRATOS_CREATE_LOCAL_FLAG, but it does not play nice with template parameters
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
const Kratos::Flags TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::SILENT_WARNINGS(Kratos::Flags::Create(0));

///@}

} /* namespace Kratos.*/
