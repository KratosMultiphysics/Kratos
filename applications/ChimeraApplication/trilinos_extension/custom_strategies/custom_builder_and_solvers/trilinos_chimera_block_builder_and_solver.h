//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//
#if !defined(KRATOS_CHIMERA_TRILINOS_BLOCK_BUILDER_AND_SOLVER)
#define KRATOS_CHIMERA_TRILINOS_BLOCK_BUILDER_AND_SOLVER

/* System includes */
#include <unordered_set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/timer.h"
#include "includes/master_slave_constraint.h"

/* Application includes */
#include "custom_utilities/helper_classes_for_constraint_builder.h"

/* Trilinos includes */
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"

#define START_TIMER(label, rank) \
    if (mrComm.MyPID() == rank)  \
        Timer::Start(label);
#define STOP_TIMER(label, rank) \
    if (mrComm.MyPID() == rank) \
        Timer::Stop(label);

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
 * @class TrilinosChimeraResidualBasedBlockBuilderAndSolver
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
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosChimeraResidualBasedBlockBuilderAndSolver
    : public BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosChimeraResidualBasedBlockBuilderAndSolver);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    typedef MasterSlaveConstraint MasterSlaveConstraintType;
    typedef typename MasterSlaveConstraint::Pointer MasterSlaveConstraintPointerType;
    typedef Internals::AuxiliaryGlobalMasterSlaveConstraint AuxiliaryGlobalMasterSlaveConstraintType;
    typedef Internals::GlobalMasterSlaveRelationContainerType GlobalMasterSlaveRelationContainerType;
    typedef std::vector<IndexType> EquationIdVectorType;
    typedef std::vector<IndexType> VectorIndexType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef Vector VectorType;
    typedef Internals::ConstraintImposer<TSparseSpace, TDenseSpace, TLinearSolver> ConstraintImposerType;

    /// Epetra definitions
    typedef Epetra_MpiComm EpetraCommunicatorType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef DofType::Pointer DofPointerType;
    /// The DoF pointer vector type definition
    typedef std::vector< DofType::Pointer > DofPointerVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    TrilinosChimeraResidualBasedBlockBuilderAndSolver(EpetraCommunicatorType &rComm,
                                         int GuessRowSize,
                                         typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver),
          mrComm(rComm),
          mGuessRowSize(GuessRowSize)
    {
    }

    /**
     * @brief Default destructor.
     */
    ~TrilinosChimeraResidualBasedBlockBuilderAndSolver() override = default;

    /**
     * Copy constructor
     */
    TrilinosChimeraResidualBasedBlockBuilderAndSolver(const TrilinosChimeraResidualBasedBlockBuilderAndSolver &rOther) = delete;

    /**
     * Assignment operator
     */
    TrilinosChimeraResidualBasedBlockBuilderAndSolver &operator=(const TrilinosChimeraResidualBasedBlockBuilderAndSolver &rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    //**************************************************************************
    //**************************************************************************

    void InitializeSolutionStep(
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo &r_process_info = rModelPart.GetProcessInfo();
        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
#pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k)
        {
            auto it = constraints_begin + k;
            it->InitializeSolutionStep(r_process_info); // Here each constraint constructs and stores its T and C matrices. Also its equation slave_ids.
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep(
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo &r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
#pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k)
        {
            auto it = constraints_begin + k;
            it->FinalizeSolutionStep(r_process_info);
        }
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
                  ModelPart &rModelPart,
                  TSystemMatrixType &rA) override
    {
        KRATOS_TRY
        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_ids_vector;
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();

        // assemble all elements
        for (auto it = rModelPart.Elements().ptr_begin(); it < rModelPart.Elements().ptr_end(); it++)
        {
            pScheme->Calculate_LHS_Contribution(*(it), LHS_Contribution,
                                                equation_ids_vector, r_current_process_info);

            // assemble the elemental contribution
            TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);

            // clean local elemental memory
            pScheme->CleanMemory(*(it));
        }

        LHS_Contribution.resize(0, 0, false);

        // assemble all conditions
        for (auto it = rModelPart.Conditions().ptr_begin(); it < rModelPart.Conditions().ptr_end(); it++)
        {
            // calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(
                *(it), LHS_Contribution, equation_ids_vector, r_current_process_info);

            // assemble the elemental contribution
            TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
        }

        // finalizing the assembly
        rA.GlobalAssemble();
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
                                     ModelPart &rModelPart,
                                     TSystemMatrixType &A) override
    {
        KRATOS_ERROR << "Method BuildLHS_CompleteOnFreeRows not implemented in "
                        "Trilinos Builder And Solver"
                     << std::endl;
    }

    /**
     * @brief This is a call to the linear system solver
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void SystemSolveWithPhysics(TSystemMatrixType &rA,
                                TSystemVectorType &rDx,
                                TSystemVectorType &rb,
                                ModelPart &rModelPart)
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(rb) != 0)
            norm_b = TSparseSpace::TwoNorm(rb);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            if (BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded())
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(
                    rA, rDx, rb, BaseType::mDofSet, rModelPart);

            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        }
        else
        {
            TSparseSpace::SetToZero(rDx);
            KRATOS_WARNING(
                "TrilinosResidualBasedBlockBuilderAndSolver")
                << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // prints informations about the current time
        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", (BaseType::GetEchoLevel() > 1))
            << *(BaseType::mpLinearSystemSolver) << std::endl;

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
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb) override
    {
        KRATOS_TRY

        if (BaseType::GetEchoLevel() > 0)
            START_TIMER("UpdateConstraints", 0)
        this->UpdateConstraintsForBuilding(rModelPart);
        if (BaseType::GetEchoLevel() > 0)
            STOP_TIMER("UpdateConstraints", 0)

        if (BaseType::GetEchoLevel() > 0)
            START_TIMER("Build", 0)
        Build(pScheme, rModelPart, rA, rb);
        if (BaseType::GetEchoLevel() > 0)
            STOP_TIMER("Build", 0)

        this->ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nBefore the solution of the system"
            << "\nSystem Matrix = " << rA << "\nunknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        if (BaseType::GetEchoLevel() > 0)
            START_TIMER("System solve time ", 0);
        this->SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nAfter the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        if (BaseType::GetEchoLevel() > 0)
            START_TIMER("Reconstruct Slaves time ", 0);
        ReconstructSlaveSolutionAfterSolve(rModelPart, rA, rDx, rb);
        if (BaseType::GetEchoLevel() > 0)
            STOP_TIMER("Reconstruct Slaves time ", 0);

        if (BaseType::GetEchoLevel() > 0)
            STOP_TIMER("System solve time ", 0);

        KRATOS_CATCH("")
    }

    /*
    * This function is exactly same as the Build() function in base class except that the function
    * has the call to ApplyConstraints function call once the LHS or RHS are computed by elements and conditions
    */
    void Build(typename TSchemeType::Pointer pScheme,
                              ModelPart &rModelPart,
                              TSystemMatrixType &rA,
                              TSystemVectorType &rb) override
    {
        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector);
        ConstraintImposerType constraint_imposer(mGlobalMasterSlaveConstraints);

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_ids_vector;
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        // assemble all elements
        for (auto it = rModelPart.Elements().ptr_begin(); it < rModelPart.Elements().ptr_end(); it++)
        {
            // detect if the element is active or not. If the user did not make
            // any choice the element is active by default
            const bool element_is_active = !((*it)->IsDefined(ACTIVE)) || (*it)->Is(ACTIVE);

            if (element_is_active)
            {
                // calculate elemental contribution
                pScheme->CalculateSystemContributions(
                    *(it), LHS_Contribution, RHS_Contribution,
                    equation_ids_vector, r_current_process_info);

                constraint_imposer.template ApplyConstraints<Element>(*(*it), LHS_Contribution,
                                                                      RHS_Contribution,
                                                                      equation_ids_vector,
                                                                      r_current_process_info);

                // assemble the elemental contribution
                TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
                TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);

                // clean local elemental memory
                pScheme->CleanMemory(*(it));
            }
        }

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // assemble all conditions
        for (auto it = rModelPart.Conditions().ptr_begin(); it < rModelPart.Conditions().ptr_end(); it++)
        {
            // detect if the element is active or not. If the user did not make
            // any choice the element is active by default
            const bool condition_is_active = !((*it)->IsDefined(ACTIVE)) || (*it)->Is(ACTIVE);
            if (condition_is_active)
            {
                // calculate elemental contribution
                pScheme->Condition_CalculateSystemContributions(
                    *(it), LHS_Contribution, RHS_Contribution,
                    equation_ids_vector, r_current_process_info);

                constraint_imposer.template ApplyConstraints<Condition>(*(*it), LHS_Contribution,
                                                                        RHS_Contribution,
                                                                        equation_ids_vector,
                                                                        r_current_process_info);

                // assemble the condition contribution
                TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
                TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);

                // clean local elemental memory
                pScheme->CleanMemory(*(it));
            }
        }

        // finalizing the assembly
        rA.GlobalAssemble();
        rb.GlobalAssemble();
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
                          ModelPart &rModelPart,
                          TSystemMatrixType &rA,
                          TSystemVectorType &rDx,
                          TSystemVectorType &rb) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, rb);
        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

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
                  ModelPart &rModelPart,
                  TSystemVectorType &rb) override
    {
        KRATOS_TRY
        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_ids_vector;
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();

        // assemble all elements
        for (auto it = rModelPart.Elements().ptr_begin(); it < rModelPart.Elements().ptr_end(); it++)
        {
            // calculate elemental Right Hand Side Contribution
            pScheme->Calculate_RHS_Contribution(*(it), RHS_Contribution,
                                                equation_ids_vector, r_current_process_info);

            // assemble the elemental contribution
            TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);
        }

        RHS_Contribution.resize(0, false);

        // assemble all conditions
        for (auto it = rModelPart.Conditions().ptr_begin(); it < rModelPart.Conditions().ptr_end(); it++)
        {
            // calculate elemental contribution
            pScheme->Condition_Calculate_RHS_Contribution(
                *(it), RHS_Contribution, equation_ids_vector, r_current_process_info);

            // assemble the elemental contribution
            TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);
        }

        // finalizing the assembly
        rb.GlobalAssemble();

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
    void SetUpDofSet(typename TSchemeType::Pointer pScheme, ModelPart &rModelPart) override
    {
        KRATOS_TRY
        typedef Element::DofsVectorType DofsVectorType;
        // Gets the array of elements from the modeler
        ElementsArrayType &r_elements_array =
            rModelPart.GetCommunicator().LocalMesh().Elements();
        DofsVectorType dof_list;
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();

        DofsArrayType temp_dofs_array;
        IndexType guess_num_dofs =
            rModelPart.GetCommunicator().LocalMesh().NumberOfNodes() * 3;
        temp_dofs_array.reserve(guess_num_dofs);
        BaseType::mDofSet = DofsArrayType();

        // Taking dofs of elements
        for (auto it_elem = r_elements_array.ptr_begin(); it_elem != r_elements_array.ptr_end(); ++it_elem)
        {
            const bool element_is_active = !((*it_elem)->IsDefined(ACTIVE)) || (*it_elem)->Is(ACTIVE);

            if (element_is_active)
            {
                pScheme->GetElementalDofList(*(it_elem), dof_list, r_current_process_info);
                for (typename DofsVectorType::iterator i_dof = dof_list.begin();
                    i_dof != dof_list.end(); ++i_dof)
                    temp_dofs_array.push_back(*i_dof);
            }
        }

        // Taking dofs of conditions
        auto &r_conditions_array = rModelPart.Conditions();
        for (auto it_cond = r_conditions_array.ptr_begin(); it_cond != r_conditions_array.ptr_end(); ++it_cond)
        {
            const bool cond_is_active = !((*it_cond)->IsDefined(ACTIVE)) || (*it_cond)->Is(ACTIVE);

            if (cond_is_active)
            {
                pScheme->GetConditionDofList(*(it_cond), dof_list, r_current_process_info);
                for (typename DofsVectorType::iterator i_dof = dof_list.begin();
                    i_dof != dof_list.end(); ++i_dof)
                    temp_dofs_array.push_back(*i_dof);
            }
        }

        temp_dofs_array.Unique();
        BaseType::mDofSet = temp_dofs_array;

        // throws an exception if there are no Degrees of freedom involved in
        // the analysis
        if (BaseType::mDofSet.size() == 0)
            KRATOS_ERROR << "No degrees of freedom!";

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have
        // reactions defined This is to be done only in debug mode
        if (BaseType::GetCalculateReactionsFlag())
        {
            for (auto dof_iterator = BaseType::mDofSet.begin();
                 dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction())
                    << "Reaction variable not set for the following : " << std::endl
                    << "Node : " << dof_iterator->Id() << std::endl
                    << "Dof : " << (*dof_iterator) << std::endl
                    << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif
        BaseType::mDofSetIsInitialized = true;

        KRATOS_CATCH("")
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     *          Sets equation id for degrees of freedom
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(ModelPart &rModelPart) override
    {
        int free_size = 0;
        auto &r_comm = rModelPart.GetCommunicator();
        const auto &r_data_comm = r_comm.GetDataCommunicator();
        int current_rank = r_comm.MyPID();

        // Calculating number of fixed and free dofs
        for (const auto &r_dof : BaseType::mDofSet)
            if (r_dof.GetSolutionStepValue(PARTITION_INDEX) == current_rank)
                free_size++;

        // Calculating the total size and required offset
        int free_offset;
        int global_size;

        // The correspounding offset by the sum of the sizes in thread with
        // inferior current_rank
        free_offset = r_data_comm.ScanSum(free_size);

        // The total size by the sum of all size in all threads
        global_size = r_data_comm.SumAll(free_size);

        // finding the offset for the begining of the partition
        free_offset -= free_size;

        // Now setting the equation id with .
        for (auto &r_dof : BaseType::mDofSet)
            if (r_dof.GetSolutionStepValue(PARTITION_INDEX) == current_rank)
                r_dof.SetEquationId(free_offset++);

        BaseType::mEquationSystemSize = global_size;
        mLocalSystemSize = free_size;
        KRATOS_INFO_IF_ALL_RANKS("TrilinosChimeraResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() > 0)
            << std::endl
            << current_rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize
            << std::endl
            << current_rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl
            << current_rank << " : free_offset = " << free_offset << std::endl;

        // by Riccardo ... it may be wrong!
        mFirstMyId = free_offset - mLocalSystemSize;
        mLastMyId = mFirstMyId + mLocalSystemSize;

        r_comm.SynchronizeDofs();

        // For Constraints
        FormulateGlobalMasterSlaveRelations(rModelPart);

    }

    /**
     * @brief Resizes the system matrix and the vector according to the number
     * of dos in the current rModelPart. This function also decides on the
     * sparsity pattern and the graph of the trilinos csr matrix
     * @param pScheme The integration scheme considered
     * @param rpA The LHS matrix
     * @param rpDx The Unknowns vector
     * @param rpd The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void ResizeAndInitializeVectors(typename TSchemeType::Pointer pScheme,
                                                   TSystemMatrixPointerType &rpA,
                                                   TSystemVectorPointerType &rpDx,
                                                   TSystemVectorPointerType &rpb,
                                                   ModelPart &rModelPart) override
    {
        // resizing the system vectors and matrix
        if (rpA == nullptr || TSparseSpace::Size1(*rpA) == 0 ||
            BaseType::GetReshapeMatrixFlag() == true) // if the matrix is not initialized
        {
            ConstraintImposerType constraint_imposer(mGlobalMasterSlaveConstraints);
            IndexType number_of_local_dofs = mLastMyId - mFirstMyId;
            int temp_size = number_of_local_dofs;
            if (temp_size < 1000)
                temp_size = 1000;
            std::vector<int> temp(temp_size, 0);
            // TODO: Check if these should be local elements and conditions
            auto &r_elements_array = rModelPart.Elements();
            auto &r_conditions_array = rModelPart.Conditions();
            // generate map - use the "temp" array here
            for (IndexType i = 0; i != number_of_local_dofs; i++)
                temp[i] = mFirstMyId + i;
            Epetra_Map my_map(-1, number_of_local_dofs, temp.data(), 0, mrComm);
            // create and fill the graph of the matrix --> the temp array is
            // reused here with a different meaning
            Epetra_FECrsGraph Agraph(Copy, my_map, mGuessRowSize);
            Element::EquationIdVectorType equation_ids_vector;
            ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();

            // assemble all elements
            for (auto it_elem = r_elements_array.ptr_begin(); it_elem != r_elements_array.ptr_end(); ++it_elem)
            {

                const bool element_is_active = !((*it_elem)->IsDefined(ACTIVE)) || (*it_elem)->Is(ACTIVE);

                if (!element_is_active)
                {
                    continue;
                }
                pScheme->EquationId(*(it_elem), equation_ids_vector,
                                    r_current_process_info);

                constraint_imposer.template ApplyConstraints<Element>(*(*(it_elem)), equation_ids_vector,
                                                                      r_current_process_info);
                // filling the list of active global indices (non fixed)
                IndexType num_active_indices = 0;
                for (IndexType i = 0; i < equation_ids_vector.size(); i++)
                {
                    temp[num_active_indices] = equation_ids_vector[i];
                    num_active_indices += 1;
                }

                if (num_active_indices != 0)
                {
                    int ierr = Agraph.InsertGlobalIndices(
                        num_active_indices, temp.data(), num_active_indices, temp.data());
                    KRATOS_ERROR_IF(ierr < 0)
                        << ": Epetra failure in Graph.InsertGlobalIndices. "
                           "Error code: "
                        << ierr << std::endl;
                }
                std::fill(temp.begin(), temp.end(), 0);
            }

            // assemble all conditions
            for (auto it_cond = r_conditions_array.ptr_begin(); it_cond != r_conditions_array.ptr_end(); ++it_cond)
            {
                const bool cond_is_active = !((*it_cond)->IsDefined(ACTIVE)) || (*it_cond)->Is(ACTIVE);

                if (!cond_is_active)
                {
                    continue;
                }
                pScheme->Condition_EquationId(
                    *(it_cond), equation_ids_vector, r_current_process_info);

                constraint_imposer.template ApplyConstraints<Condition>(*(*(it_cond)), equation_ids_vector, r_current_process_info);

                // filling the list of active global indices (non fixed)
                IndexType num_active_indices = 0;
                for (IndexType i = 0; i < equation_ids_vector.size(); i++)
                {
                    temp[num_active_indices] = equation_ids_vector[i];
                    num_active_indices += 1;
                }

                if (num_active_indices != 0)
                {
                    int ierr = Agraph.InsertGlobalIndices(
                        num_active_indices, temp.data(), num_active_indices, temp.data());
                    KRATOS_ERROR_IF(ierr < 0)
                        << ": Epetra failure in Graph.InsertGlobalIndices. "
                           "Error code: "
                        << ierr << std::endl;
                }
                std::fill(temp.begin(), temp.end(), 0);
            }


            // finalizing graph construction
            int ierr = Agraph.GlobalAssemble();
            KRATOS_ERROR_IF(ierr < 0)
                << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr
                << std::endl;
            // generate a new matrix pointer according to this graph
            TSystemMatrixPointerType p_new_A =
                TSystemMatrixPointerType(new TSystemMatrixType(Copy, Agraph));
            rpA.swap(p_new_A);
            // generate new vector pointers according to the given map
            if (rpb == nullptr || TSparseSpace::Size(*rpb) != BaseType::mEquationSystemSize)
            {
                TSystemVectorPointerType p_new_b =
                    TSystemVectorPointerType(new TSystemVectorType(my_map));
                rpb.swap(p_new_b);
            }
            if (rpDx == nullptr || TSparseSpace::Size(*rpDx) != BaseType::mEquationSystemSize)
            {
                TSystemVectorPointerType p_new_Dx =
                    TSystemVectorPointerType(new TSystemVectorType(my_map));
                rpDx.swap(p_new_Dx);
            }
            if (BaseType::mpReactionsVector == nullptr) // if the pointer is not initialized initialize it to an
                                                        // empty matrix
            {
                TSystemVectorPointerType pNewReactionsVector =
                    TSystemVectorPointerType(new TSystemVectorType(my_map));
                BaseType::mpReactionsVector.swap(pNewReactionsVector);
            }
        }
        else if (BaseType::mpReactionsVector == nullptr && this->mCalculateReactionsFlag)
        {
            TSystemVectorPointerType pNewReactionsVector =
                TSystemVectorPointerType(new TSystemVectorType(rpDx->Map()));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }
        else
        {
            if (TSparseSpace::Size1(*rpA) == 0 ||
                TSparseSpace::Size1(*rpA) != BaseType::mEquationSystemSize ||
                TSparseSpace::Size2(*rpA) != BaseType::mEquationSystemSize)
            {
                KRATOS_ERROR
                    << "It should not come here resizing is not allowed this "
                       "way!!!!!!!! ... ";
            }
        }
    }

    //**************************************************************************
    //**************************************************************************
    void CalculateReactions(typename TSchemeType::Pointer pScheme,
                            ModelPart &rModelPart,
                            TSystemMatrixType &rA,
                            TSystemVectorType &rDx,
                            TSystemVectorType &rb) override
    {
        TSparseSpace::SetToZero(rb);

        // refresh RHS to have the correct reactions
        BuildRHS(pScheme, rModelPart, rb);

        // initialize the Epetra importer
        // TODO: this part of the code has been pasted until a better solution
        // is found
        int system_size = TSparseSpace::Size(rb);
        int number_of_dofs = BaseType::mDofSet.size();
        std::vector<int> index_array(number_of_dofs);

        // filling the array with the global ids
        int counter = 0;
        int id = 0;
        for (const auto &dof : BaseType::mDofSet)
        {
            id = dof.EquationId();
            if (id < system_size)
            {
                index_array[counter++] = id;
            }
        }

        std::sort(index_array.begin(), index_array.end());
        std::vector<int>::iterator NewEnd =
            std::unique(index_array.begin(), index_array.end());
        index_array.resize(NewEnd - index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        rb.Comm().SumAll(&tot_update_dofs, &check_size, 1);
        KRATOS_ERROR_IF(check_size < system_size)
            << "Dof count is not correct. There are less dofs than expected.\n"
            << "Expected number of active dofs = " << system_size
            << " dofs found = " << check_size;

        // defining a map as needed
        Epetra_Map dof_update_map(-1, index_array.size(),
                                  &(*(index_array.begin())), 0, rb.Comm());

        // defining the importer class
        Kratos::shared_ptr<Epetra_Import> pDofImporter =
            Kratos::make_shared<Epetra_Import>(dof_update_map, rb.Map());

        // defining a temporary vector to gather all of the values needed
        Epetra_Vector temp_RHS(pDofImporter->TargetMap());

        // importing in the new temp_RHS vector the values
        int ierr = temp_RHS.Import(rb, *pDofImporter, Insert);
        if (ierr != 0)
            KRATOS_ERROR << "Epetra failure found - error code: " << ierr;

        double *temp_RHS_values; // DO NOT make delete of this one!!
        temp_RHS.ExtractView(&temp_RHS_values);

        rb.Comm().Barrier();

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        // store the RHS values in the reaction variable
        // NOTE: dofs are assumed to be numbered consecutively in the
        // BlockBuilderAndSolver
        for (int k = 0; k < ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator =
                BaseType::mDofSet.begin() + k;

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
    void ApplyDirichletConditions(typename TSchemeType::Pointer pScheme,
                                  ModelPart &rModelPart,
                                  TSystemMatrixType &rA,
                                  TSystemVectorType &rDx,
                                  TSystemVectorType &rb) override
    {
        KRATOS_TRY

        // loop over all dofs to find the fixed ones
        std::vector<int> global_ids(BaseType::mDofSet.size());
        std::vector<int> is_dirichlet(BaseType::mDofSet.size());

        IndexType i = 0;
        for (const auto &dof : BaseType::mDofSet)
        {
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

        // importing in the new temp vector the values
        int ierr = fixed.Import(fixed_local, dirichlet_importer, Insert);
        if (ierr != 0)
            KRATOS_ERROR << "Epetra failure found";

        for (int i = 0; i < rA.NumMyRows(); i++)
        {
            int numEntries; // number of non-zero entries
            double *vals;   // row non-zero values
            int *cols;      // column indices of row non-zero values
            rA.ExtractMyRowView(i, numEntries, vals, cols);

            int row_gid = rA.RowMap().GID(i);
            int row_lid = localmap.LID(row_gid);

            if (fixed_local[row_lid] == 0) // not a dirichlet row
            {
                for (int j = 0; j < numEntries; j++)
                {
                    if (fixed[cols[j]] == true)
                        vals[j] = 0.0;
                }
            }
            else // this IS a dirichlet row
            {
                // set to zero the rhs
                rb[0][i] = 0.0; // note that the index of i is expected to be
                                // coherent with the rows of A

                // set to zero the whole row
                for (int j = 0; j < numEntries; j++)
                {
                    int col_gid = rA.ColMap().GID(cols[j]);
                    if (col_gid != row_gid)
                        vals[j] = 0.0;
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    EpetraCommunicatorType &mrComm;
    int mGuessRowSize;
    IndexType mLocalSystemSize;
    int mFirstMyId;
    int mLastMyId;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * @brief   this method condenses the MasterSlaveConstraints which are added on the rModelPart
     *          into objects of AuxilaryGlobalMasterSlaveRelation. One unique object for each unique slave.
     *          these will be used in the ApplyConstraints functions later on.
     * @param   rModelPart The model part of the problem to solve
     */
    void FormulateGlobalMasterSlaveRelations(ModelPart &rModelPart)
    {
        KRATOS_TRY
        MasterSlaveConstraintType::DofPointerVectorType slave_dofs_vector;
        MasterSlaveConstraintType::DofPointerVectorType master_dofs_vector;
        // First delete the existing ones
        mGlobalMasterSlaveConstraints.clear();
        int assembled_count = 0;
        // Getting the array of the conditions
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        // Getting the beginning iterator

        const ModelPart::MasterSlaveConstraintContainerType::iterator constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();

        for (int i_constraints = 0; i_constraints < number_of_constraints; i_constraints++)
        {
            ModelPart::MasterSlaveConstraintContainerType::iterator it = constraints_begin;
            std::advance(it, i_constraints);
            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool constraint_is_active = true;
            if ((it)->IsDefined(ACTIVE))
                constraint_is_active = (it)->Is(ACTIVE);

            if (constraint_is_active)
            {
                //assemble the Constraint contribution
                AssembleConstraint(*it, r_current_process_info);
                ++assembled_count;
            }
        }
        KRATOS_CATCH("TrilinosChimeraResidualBasedBlockBuilderAndSolver::FormulateGlobalMasterSlaveRelations failed ..");
    }

    /**
     * @brief   this method assembles the given master slave constraint to the auxiliary global master slave constraints
     * @param   rMasterSlaveConstraint object of the master slave constraint to be assembled.
     * @param   rCurrentProcessInfo current process info.
     */
    void AssembleConstraint(ModelPart::MasterSlaveConstraintType &rMasterSlaveConstraint, ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY
        int slave_count = 0;
        LocalSystemMatrixType relation_matrix(0, 0);
        LocalSystemVectorType constant_vector(0);
        EquationIdVectorType slave_equation_ids(0);
        EquationIdVectorType master_equation_ids(0);

        //get the equation Ids of the constraint
        rMasterSlaveConstraint.EquationIdVector(slave_equation_ids, master_equation_ids, rCurrentProcessInfo);
        //calculate constraint's T and b matrices
        rMasterSlaveConstraint.CalculateLocalSystem(relation_matrix, constant_vector, rCurrentProcessInfo);

        for (auto slave_equation_id : slave_equation_ids)
        {
            int master_count = 0;
            auto global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            if (global_constraint == mGlobalMasterSlaveConstraints.end())
            {
                mGlobalMasterSlaveConstraints[slave_equation_id] = Kratos::make_unique<AuxiliaryGlobalMasterSlaveConstraintType>(slave_equation_id);
            }

            global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            for (auto master_equation_id : master_equation_ids)
            {
                global_constraint->second->AddMaster(master_equation_id, relation_matrix(slave_count, master_count));
                master_count++;
            }
            slave_count++;
        }
        KRATOS_CATCH("TrilinosChimeraResidualBasedBlockBuilderAndSolver::AssembleSlaves failed ...");
    }

    /**
     * @brief   this method resets the LHS and RHS values of the AuxilaryGlobalMasterSlaveRelation objects
     */
    void ResetConstraintRelations()
    {
        KRATOS_TRY
        const int number_of_constraints = static_cast<int>(mGlobalMasterSlaveConstraints.size());

        // Getting the beginning iterator
        const GlobalMasterSlaveRelationContainerType::iterator constraints_begin = mGlobalMasterSlaveConstraints.begin();

        for (int i_constraints = 0; i_constraints < number_of_constraints; ++i_constraints)
        {
            GlobalMasterSlaveRelationContainerType::iterator it = constraints_begin;
            std::advance(it, i_constraints);

            (it->second)->Reset();
        }
        KRATOS_CATCH("TrilinosChimeraResidualBasedBlockBuilderAndSolver::ResetConstraintRelations failed to reset constraint relations..");
    }

    /**
     * @brief   this method uses the MasterSlaveConstraints objects in rModelPart to reconstruct the LHS and RHS values
     *          of the AuxilaryGlobalMasterSlaveRelation objects. That is the value of Slave as LHS and the T*M+C as RHS value
     * @param   rModelPart The model part of the problem to solve
     */
    void UpdateConstraintsForBuilding(ModelPart &rModelPart)
    {
        KRATOS_TRY
        Communicator& r_comm = rModelPart.GetCommunicator();
        r_comm.SynchronizeNodalSolutionStepsData();
        // Reset the constraint equations
        ResetConstraintRelations();
        // Getting the array of the conditions
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        // Getting the beginning iterator
        const ModelPart::MasterSlaveConstraintContainerType::iterator constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();

#pragma omp parallel for schedule(guided, 512)
        for (int i_constraints = 0; i_constraints < number_of_constraints; i_constraints++)
        {
            ModelPart::MasterSlaveConstraintContainerType::iterator it = constraints_begin;
            std::advance(it, i_constraints);

            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool constraint_is_active = true;
            if ((it)->IsDefined(ACTIVE))
                constraint_is_active = (it)->Is(ACTIVE);

            if (constraint_is_active)
            {
                UpdateMasterSlaveConstraint(*it, r_current_process_info);
            }
        }
        KRATOS_CATCH("TrilinosChimeraResidualBasedBlockBuilderAndSolver::UpdateConstraintsForBuilding failed ..");
    }

    /**
     * @brief   this method uses the MasterSlaveConstraints objects in rModelPart to reconstruct the LHS and RHS values
     *          of the individual AuxilaryGlobalMasterSlaveRelation object. That is the value of Slave as LHS and the T*M+C as RHS value
     * @param   rMasterSlaveConstraint The MasterSlaveConstraint which is to be updated
     */
    void UpdateMasterSlaveConstraint(ModelPart::MasterSlaveConstraintType &rMasterSlaveConstraint, ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY
        //contributions to the system
        LocalSystemMatrixType relation_matrix(0, 0);
        LocalSystemVectorType constant_vector(0);
        EquationIdVectorType slave_equation_ids(0);
        EquationIdVectorType master_equation_ids(0);

        //get the equation Ids of the constraint
        rMasterSlaveConstraint.EquationIdVector(slave_equation_ids, master_equation_ids, rCurrentProcessInfo);
        //calculate constraint's T and b matrices
        rMasterSlaveConstraint.CalculateLocalSystem(relation_matrix, constant_vector, rCurrentProcessInfo);
        // For calculating the constant
        MasterSlaveConstraintType::DofPointerVectorType slave_dofs_vector;
        MasterSlaveConstraintType::DofPointerVectorType master_dofs_vector;
        rMasterSlaveConstraint.GetDofList(slave_dofs_vector, master_dofs_vector, rCurrentProcessInfo);

        int slave_index = 0;
        for (auto &slave_dof : slave_dofs_vector)
        {
            double slave_value_calc = 0.0;
            for (IndexType master_index = 0; master_index < master_dofs_vector.size(); master_index++)
            {
                slave_value_calc += master_dofs_vector[master_index]->GetSolutionStepValue() * relation_matrix(slave_index, master_index);
            }
            slave_value_calc += constant_vector[slave_index];
            auto global_constraint = mGlobalMasterSlaveConstraints.find(slave_dof->EquationId());
            global_constraint->second->SetLeftHandSide(slave_dof->GetSolutionStepValue());
            global_constraint->second->UpdateRightHandSide(slave_value_calc);
            slave_index++;
        }
        KRATOS_CATCH("TrilinosChimeraResidualBasedBlockBuilderAndSolver::UpdateMasterSlaveConstraint failed ..");
    }

    /**
     * @brief This method reconstructs the slave solution after Solving.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void ReconstructSlaveSolutionAfterSolve(
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb)
    {
        KRATOS_TRY
        const int number_of_global_constraints = static_cast<int>(mGlobalMasterSlaveConstraints.size());
        // Getting the beginning iterator

        const GlobalMasterSlaveRelationContainerType::iterator global_constraints_begin = mGlobalMasterSlaveConstraints.begin();
        //contributions to the system
        VectorType master_weights_vector;
        double constant = 0.0;
        const int num_local_elems = rDx.MyLength();
        // std::vector<int> ids_to_gather(num_local_elems, 0);
        std::vector<int> ids_to_gather;
        ids_to_gather.reserve(num_local_elems);
        // rDx.Map().MyGlobalElements(&ids_to_gather[0]);

        IndexType slave_equation_id = 0;
        EquationIdVectorType master_equation_ids = EquationIdVectorType(0);
        DofsVectorType slave_dof_list, master_dof_list;

        for (int i_constraints = 0; i_constraints < number_of_global_constraints; i_constraints++)
        {
            GlobalMasterSlaveRelationContainerType::iterator it = global_constraints_begin;
            std::advance(it, i_constraints);

            //get the equation Ids of the constraint
            (it->second)->EquationIdsVector(slave_equation_id, master_equation_ids);
            for(const auto& master_eq_id : master_equation_ids){
                ids_to_gather.push_back( master_eq_id );
            }
        }
        ids_to_gather.shrink_to_fit();

        // Here we should remove the duplicates in ids_to_gather otherwise it will segfault
        std::sort(ids_to_gather.begin(), ids_to_gather.end());
        ids_to_gather.erase( std::unique( ids_to_gather.begin(), ids_to_gather.end() ), ids_to_gather.end() );
        double gathered_values[ids_to_gather.size()];
        TSparseSpace::GatherValues(rDx, ids_to_gather, gathered_values);

        std::vector<int> slave_eq_ids;
        slave_eq_ids.reserve(number_of_global_constraints);
        std::vector<double> slave_values;
        slave_values.reserve(number_of_global_constraints);

        for (int i_constraints = 0; i_constraints < number_of_global_constraints; i_constraints++)
        {
            GlobalMasterSlaveRelationContainerType::iterator it = global_constraints_begin;
            std::advance(it, i_constraints);
            // (it->second)->PrintInfo(std::cout);

            double slave_dx_value = 0.0;
            //get the equation Ids of the constraint
            (it->second)->EquationIdsVector(slave_equation_id, master_equation_ids);
            //calculate constraint's T and b matrices
            (it->second)->CalculateLocalSystem(master_weights_vector, constant);
            int master_index = 0;
            for (auto &master_equation_id : master_equation_ids)
            {
                auto result = FindInVector(ids_to_gather, static_cast<int>(master_equation_id));
                double master_value = 0.0;
                if(result.first) // This is a remote master
                     master_value = gathered_values[result.second];
                else // This is a local master
                    KRATOS_INFO_ALL_RANKS("3. master not found              : ")<<master_equation_id<<std::endl;

                slave_dx_value += master_value * master_weights_vector(master_index);
                master_index++;
            }
            slave_dx_value += constant;
            const int  slave_eq_id = static_cast<int>(slave_equation_id);
            slave_values.push_back(slave_dx_value);
            slave_eq_ids.push_back(slave_eq_id);
            // KRATOS_INFO_ALL_RANKS("slave dx value : ")<<slave_dx_value<<std::endl;
        }
        int ierr = rDx.ReplaceGlobalValues (number_of_global_constraints, slave_eq_ids.data(), slave_values.data());
        KRATOS_ERROR_IF(ierr < 0) << "Epetra failure when attempting to ReplaceGlobalValues in ReconstructSlave function." << std::endl;

        ierr = rDx.GlobalAssemble(Insert); //Epetra_CombineMode mode=Add);
        KRATOS_ERROR_IF(ierr < 0) << "Epetra failure when attempting to insert value in function SetValue" << std::endl;

        KRATOS_CATCH("TrilinosChimeraResidualBasedBlockBuilderAndSolver::ReconstructSlaveSolutionAfterSolve failed ..");
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

    // This is the set of condenced global constraints.
    GlobalMasterSlaveRelationContainerType mGlobalMasterSlaveConstraints; //This can be changed to more efficient implementation later on.

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    //**************************************************************************
    //**************************************************************************
    void AssembleLHS_CompleteOnFreeRows(TSystemMatrixType &rA,
                                        LocalSystemMatrixType &rLHS_Contribution,
                                        Element::EquationIdVectorType &rEquationId)
    {
        KRATOS_ERROR << "This method is not implemented for Trilinos";
    }

    /*
        Generic function to find an element in vector and also its position.
        It returns a pair of bool & int i.e.

        bool : Represents if element is present in vector or not.
        int : Represents the index of element in vector if its found else -1

    */
    template <typename T>
    std::pair<bool, int> FindInVector(const std::vector<T> &vecOfElements, const T &element)
    {
        std::pair<bool, int> result;
        // Find given element in vector
        auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
        if (it != vecOfElements.end())
        {
            result.second = distance(vecOfElements.begin(), it);
            result.first = true;
        }
        else
        {
            result.first = false;
            result.second = -1;
        }
        return result;
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
}; /* Class TrilinosChimeraResidualBasedBlockBuilderAndSolver */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_CHIMERA_TRILINOS_BLOCK_BUILDER_AND_SOLVER  defined */
