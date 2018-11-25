//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_WITH_CONSTRAINTS)
#define KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_WITH_CONSTRAINTS

/* System includes */
#include <unordered_set>
#include <unordered_map>

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/helper_classes_for_constraint_builder.h"
#include "input_output/logger.h"

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
 * @class ResidualBasedEliminationBuilderAndSolverWithConstraints
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * The system is build in the following manner. A T matrix is assembled and constant vector g is assembled too. The T matrix contains the relations of all the dofs of the system, even the nodes with no master/slave relation. Then the size is n_total x n_red
 *      The relation u = T u_red
 * Then:
 *      A_red = T^t A T
 *      b_red = T^t (b - A g)
 * @todo There is a more efficient way to asemble the system, but more costly, which is the following. In this case T will be only a relation matrix between master and slave dofs. Then n_slave x n_master: us = T um + g
 * Separating into independent dofs, master ans slave dofs:
 *      u = uu
 *          um
 *          us
 *      A = Auu Aum Aus
 *          Amu Amm Ams
 *          Asu Asm Ass
 *      b = bu
 *          bm
 *          bs
 * Finally:
 *  A_red = Auu              Aum + Aus T
 *          Amu + T^t Asu    Amm + T^t Ams^t + Ams T + T^t Ass T
 *  b_red = bu - Aus g
 *          bm - Ams g
 *
 * This system requires extra care and is more complicated and requires to compute the blocks properly
 * @author Vicente Mataix Ferrandiz
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
class ResidualBasedEliminationBuilderAndSolverWithConstraints
    : public ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedEliminationBuilderAndSolverWithConstraints
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedEliminationBuilderAndSolverWithConstraints);

    /// Definition of the base class
    typedef ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    // The size_t types
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
    typedef typename BaseType::NodeType NodeType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Additional definitions
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// DoF types definition
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    /// Set definition
    typedef std::unordered_set<IndexType> IndexSetType;

    /// MPC definitions
    typedef MasterSlaveConstraint MasterSlaveConstraintType;
    typedef typename MasterSlaveConstraint::Pointer MasterSlaveConstraintPointerType;
    typedef Internals::AuxiliaryGlobalMasterSlaveConstraint AuxiliaryGlobalMasterSlaveConstraintType;
    typedef Internals::GlobalMasterSlaveRelationContainerType GlobalMasterSlaveRelationContainerType;
    typedef std::vector<IndexType> VectorIndexType;
    typedef Vector VectorType;
    typedef Internals::ConstraintImposer<TSparseSpace, TDenseSpace, TLinearSolver> ConstraintImposerType;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedEliminationBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "reassemble_lhs"     : false,
            "rebuild_fixed_dofs" : false
        })" );
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mReassembleLHS = ThisParameters["reassemble_lhs"].GetBool();
        mRebuildFixedDoFs = ThisParameters["rebuild_fixed_dofs"].GetBool();
    }

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedEliminationBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        const bool ReassembleLHS = false,
        const bool RebuildFixedDoFs = false
        )
        : BaseType(pNewLinearSystemSolver),
          mReassembleLHS(ReassembleLHS),
          mRebuildFixedDoFs(RebuildFixedDoFs)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedEliminationBuilderAndSolverWithConstraints() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetUpSystem(ModelPart& rModelPart) override
    {
        if(rModelPart.MasterSlaveConstraints().size() > 0)
            SetUpSystemWithConstraints(rModelPart);
        else
            BaseType::SetUpSystem(rModelPart);
    }

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        ) override
    {
        if(rModelPart.MasterSlaveConstraints().size() > 0)
            BuildWithConstraints(pScheme, rModelPart, rA, rb);
        else
            BaseType::Build(pScheme, rModelPart, rA, rb);
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
        if(rModelPart.MasterSlaveConstraints().size() > 0)
            BuildAndSolveWithConstraints(pScheme, rModelPart, A, Dx, b);
        else
            BaseType::BuildAndSolve(pScheme, rModelPart, A, Dx, b);
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
     * and condition its Dofs.
     * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the
     * way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        ) override
    {
        if(rModelPart.MasterSlaveConstraints().size() > 0)
            SetUpDofSetWithConstraints(pScheme, rModelPart);
        else
            BaseType::SetUpDofSet(pScheme, rModelPart);
    }

    /**
     * @brief It applies certain operations at the system of equations at the begining of the solution step
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
     * @param rb The RHS vector of the system of equations
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->InitializeSolutionStep(r_process_info); // Here each constraint constructs and stores its T and C matrices. Also its equation slave_ids.
        }

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints failed to initialize solution step.")
    }

    /**
     * @brief It applies certain operations at the system of equations at the end of the solution step
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
     * @param rb The RHS vector of the system of equations
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->FinalizeSolutionStep(r_process_info);
        }
        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints failed to finalize solution step.")
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedEliminationBuilderAndSolverWithConstraints";
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

    TSystemMatrixPointerType mpTMatrix = NULL;        /// This is matrix containing the global relation for the constraints
    TSystemMatrixPointerType mpOldAMatrix = NULL;     /// This is matrix containing the old LHS structure
    TSystemVectorPointerType mpConstantVector = NULL; /// This is matrix containing the global relation for the constraints
    DofsArrayType mDoFSlaveSet;                       /// The set containing the slave DoF of the system
    SizeType mDoFToSolveSystemSize = 0;               /// Number of degrees of freedom of the problem to actually be solved

    bool mCleared = true;           /// If the system has been reseted
    bool mReassembleLHS = false;    /// If the LHS must be reconstructed after computing
    bool mRebuildFixedDoFs = false; /// If we rebuild the DoF of the slave DoFs assigned to master DoFs

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
    * @brief This method assembles the global relation matrix (T matrix used to impose the MPC)
    * @param rT The global relation matrix
    * @param rTransformationMatrix The local transformation contribution
    * @param rSlaveEquationId The equation id of the slave dofs
    * @param rMasterEquationId The equation id of the master dofs
    * @param rLockArray The lock of the dof
    */
    void AssembleRelationMatrix(
        TSystemMatrixType& rT,
        const LocalSystemMatrixType& rTransformationMatrix,
        const EquationIdVectorType& rSlaveEquationId,
        const EquationIdVectorType& rMasterEquationId
#ifdef USE_LOCKS_IN_ASSEMBLY
        ,std::vector< omp_lock_t >& rLockArray
#endif
        )
    {
        const SizeType local_size_1 = rTransformationMatrix.size1();

        for (IndexType i_local = 0; i_local < local_size_1; ++i_local) {
            IndexType i_global = rSlaveEquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize) {
            #ifdef USE_LOCKS_IN_ASSEMBLY
                omp_set_lock(&rLockArray[i_global]);
            #endif

                BaseType::AssembleRowContributionFreeDofs(rT, rTransformationMatrix, i_global, i_local, rMasterEquationId);

            #ifdef USE_LOCKS_IN_ASSEMBLY
                omp_unset_lock(&rLockArray[i_global]);
            #endif
            }
            //note that computation of reactions is not performed here!
        }
    }

    /**
     * @brief This method construcs the relationship between the DoF
     * @param pScheme The integration scheme
     * @param rA The LHS of the system
     * @param rModelPart The model part which defines the problem
     */
    void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& rA,
        ModelPart& rModelPart
        ) override
    {
        if(rModelPart.MasterSlaveConstraints().size() > 0)
            ConstructMatrixStructureWithConstraints(pScheme, rA, rModelPart);
        else
            BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
    }

    /**
     * @brief The same methods as the base class but with constraints
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
     * @param rb The RHS vector of the system of equations
     */
    void BuildAndSolveWithConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY

        Timer::Start("Build");

        // We do the build (after that we resize the solution vector to avoid problems)
        Build(pScheme, rModelPart, rA, rb);

        Timer::Stop("Build");

        // Now we apply the BC (in the elmination B&S does nothing)
        rDx.resize(mDoFToSolveSystemSize, false);
        this->ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) <<
        "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        // We solve the system of equations
        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");
        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();

        // We reconstruct the Unknowns vector and the residual
        const double start_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        ReconstructSlaveSolutionAfterSolve(pScheme, rModelPart, rA, rDx, rb);

        const double stop_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Reconstruct slaves time: " << stop_reconstruct_slaves - start_reconstruct_slaves << std::endl;

        // Some verbosity
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) <<
        "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element and condition its Dofs.
     * @details Equivalent to the ResidualBasedEliminationBuilderAndSolver but with constraints. The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSetWithConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the dofs" << std::endl;

        DofsVectorType dof_list, second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        typedef std::unordered_set < DofPointerType, DofPointerHasher> set_type;

        // Declaring temporal variables
        DofsArrayType dof_temp_all, dof_temp_solvable, dof_temp_slave;

        // We assign an empty dof array to our dof sets
        BaseType::mDofSet = DofsArrayType(); /// This corresponds with all the DoF of the system
        mDoFSlaveSet = DofsArrayType();    /// This corresponds with the slave (the ones not solved after compacting the system using MPC)

        /**
         * Here we declare three sets.
         * - The global set: Contains all the DoF of the system
         * - The slave set: The DoF that are not going to be solved, due to MPC formulation
         */
        set_type dof_global_set, dof_global_slave_set;

        #pragma omp parallel firstprivate(dof_list, second_dof_list)
        {
            ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // We cleate the temporal set and we reserve some space on them
            set_type dofs_tmp_set, dof_temp_slave_set;
            dofs_tmp_set.reserve(20000);
            dof_temp_slave_set.reserve(200);

            // Gets the array of elements from the modeler
            ElementsArrayType& r_elements_array = rModelPart.Elements();
            const int number_of_elements = static_cast<int>(r_elements_array.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i) {
                auto it_elem = r_elements_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetElementalDofList(*(it_elem.base()), dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of conditions from the modeler
            ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
            const int number_of_conditions = static_cast<int>(r_conditions_array.size());
            #pragma omp for  schedule(guided, 512) nowait
            for (int i = 0; i < number_of_conditions; ++i) {
                auto it_cond = r_conditions_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetConditionDofList(*(it_cond.base()), dof_list, r_current_process_info);
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
                dof_temp_slave_set.insert(dof_list.begin(), dof_list.end());
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
                dof_global_slave_set.insert(dof_temp_slave_set.begin(), dof_temp_slave_set.end());
            }
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        /// We transfer the temporal sets to our DoF set
        dof_temp_all.reserve(dof_global_set.size());
        for (auto& dof : dof_global_set) {
            dof_temp_all.push_back( dof.get() );
        }
        dof_temp_all.Sort();
        BaseType::mDofSet = dof_temp_all;

        dof_temp_slave.reserve(dof_global_slave_set.size());
        for (auto& dof : dof_global_slave_set) {
            dof_temp_slave.push_back( dof.get() );
        }
        dof_temp_slave.Sort();
        mDoFSlaveSet = dof_temp_slave;

        // Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_WARNING_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", mDoFSlaveSet.size() == 0) << "No slave degrees of freedom to solve!" << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

#ifdef USE_LOCKS_IN_ASSEMBLY
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing lock array" << std::endl;

        if (BaseType::mLockArray.size() != 0) {
            for (int i = 0; i < static_cast<int>(BaseType::mLockArray.size()); ++i) {
                omp_destroy_lock(&BaseType::mLockArray[i]);
            }
        }
        BaseType::mLockArray.resize(BaseType::mDofSet.size());

        for (int i = 0; i < static_cast<int>(BaseType::mLockArray.size()); ++i) {
            omp_init_lock(&BaseType::mLockArray[i]);
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "End of setup dof set\n" << std::endl;
#endif

    // If reactions are to be calculated, we check if all the dofs have reactions defined
    // This is tobe done only in debug mode
    #ifdef KRATOS_DEBUG
        if(BaseType::GetCalculateReactionsFlag()) {
            for(auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " << std::endl
                    << "Node : " << dof_iterator->Id()<< std::endl
                    << "Dof : " << (*dof_iterator) << std::endl << "Not possible to calculate reactions." << std::endl;
            }
        }
    #endif

        KRATOS_CATCH("");
    }

    /**
      *@brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void SystemSolveWithPhysics(
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
            norm_b = 0.0;

        if (norm_b != 0.0) {

             // Create the auxiliar dof set
             DofsArrayType aux_dof_set;
             aux_dof_set.reserve(mDoFToSolveSystemSize);
             for (auto& dof : BaseType::mDofSet) {
                 if (dof.EquationId() < BaseType::mEquationSystemSize) {
                    auto it = mDoFSlaveSet.find(dof);
                    if (it == mDoFSlaveSet.end())
                        aux_dof_set.push_back( &dof );
                 }
             }
             aux_dof_set.Sort();

             KRATOS_ERROR_IF_NOT(aux_dof_set.size() == mDoFToSolveSystemSize) << "Inconsistency (I) in system size: " << mDoFToSolveSystemSize << " vs " << aux_dof_set.size() << "\n Size dof set " << BaseType::mDofSet.size() << " vs Size slave dof set " << mDoFSlaveSet.size() << std::endl;
             KRATOS_ERROR_IF_NOT(aux_dof_set.size() == rA.size1()) << "Inconsistency (II) in system size: " << rA.size1() << " vs " << aux_dof_set.size() << "\n Size dof set " << BaseType::mDofSet.size() << " vs Size slave dof set " << mDoFSlaveSet.size() << std::endl;

            // Provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded())
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(rA, rDx, rb, aux_dof_set, rModelPart);

            // Do solve
            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        } else {
            TSparseSpace::SetToZero(rDx);
            KRATOS_WARNING_IF("ResidualBasedEliminationBuilderAndSolver", rModelPart.GetCommunicator().MyPID() == 0) << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // Prints informations about the current time
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")

    }

    /**
     * @brief This function is exactly same as the ConstructMatrixStructure() function in base class except that the function
     * @details Has the call to ApplyConstraints function call once the element and conditions compute their equation ids
     * @todo Move this method to a common class with block builder and solver with constraints
    */
    virtual void ConstructMatrixStructureWithConstraints(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& rA,
        ModelPart& rModelPart
        )
    {
        // Filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        // The total number of dof of the system
        const SizeType equation_size = BaseType::mEquationSystemSize;

        // This vector contains the indexes sets for all rows
        std::vector<IndexSetType> indices(equation_size);

        // We reserve some indexes on each row
        #pragma omp parallel for firstprivate(equation_size)
        for (int index = 0; index < static_cast<int>(equation_size); ++index)
            indices[index].reserve(40);

        /// Definition of the eqautio id vector type
        EquationIdVectorType ids(3, 0);
        EquationIdVectorType second_ids(3, 0); // NOTE: Used only on the constraints to take into account the master dofs

        #pragma omp parallel firstprivate(ids, second_ids)
        {
            // The process info
            ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // We repeat the same declaration for each thead
            std::vector<IndexSetType> temp_indexes(equation_size);

            #pragma omp for
            for (int index = 0; index < static_cast<int>(equation_size); ++index)
                temp_indexes[index].reserve(30);

            // Getting the size of the array of elements from the model
            const int number_of_elements = static_cast<int>(rModelPart.Elements().size());

            // Element initial iterator
            const auto el_begin = rModelPart.ElementsBegin();

            // We iterate over the elements
            #pragma omp for schedule(guided, 512) nowait
            for (int i_elem = 0; i_elem<number_of_elements; ++i_elem) {
                auto it_elem = el_begin + i_elem;
                pScheme->EquationId( *(it_elem.base()), ids, r_current_process_info);

                for (auto& id_i : ids) {
                    if (id_i < BaseType::mEquationSystemSize) {
                        auto& row_indices = temp_indexes[id_i];
                        for (auto& id_j : ids) {
                            if (id_j < BaseType::mEquationSystemSize) {
                                row_indices.insert(id_j);
                            }
                        }
                    }
                }
            }

            // Getting the size of the array of the conditions
            const int number_of_conditions = static_cast<int>(rModelPart.Conditions().size());

            // Condition initial iterator
            const auto cond_begin = rModelPart.ConditionsBegin();

            // We iterate over the conditions
            #pragma omp for schedule(guided, 512) nowait
            for (int i_cond = 0; i_cond<number_of_conditions; ++i_cond) {
                auto it_cond = cond_begin + i_cond;
                pScheme->Condition_EquationId( *(it_cond.base()), ids, r_current_process_info);
                for (auto& id_i : ids) {
                    if (id_i < BaseType::mEquationSystemSize) {
                        auto& row_indices = temp_indexes[id_i];
                        for (auto& id_j : ids) {
                            if (id_j < BaseType::mEquationSystemSize) {
                                row_indices.insert(id_j);
                            }
                        }
                    }
                }
            }

            // Getting the size of the array of the constraints
            const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

            // Constraint initial iterator
            const auto const_begin = rModelPart.MasterSlaveConstraints().begin();

            // We iterate over the constraints
            #pragma omp for schedule(guided, 512) nowait
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = const_begin + i_const;
                it_const->EquationIdVector(ids, second_ids, r_current_process_info);
                // Slave DoFs
                for (auto& id_i : ids) {
                    if (id_i < BaseType::mEquationSystemSize) {
                        auto& row_indices = temp_indexes[id_i];
                        for (auto& id_j : ids) {
                            if (id_j < BaseType::mEquationSystemSize) {
                                row_indices.insert(id_j);
                            }
                        }
                    }
                }
                // Master DoFs
                for (auto& id_i : second_ids) {
                    if (id_i < BaseType::mEquationSystemSize) {
                        auto& row_indices = temp_indexes[id_i];
                        for (auto& id_j : second_ids) {
                            if (id_j < BaseType::mEquationSystemSize) {
                                row_indices.insert(id_j);
                            }
                        }
                    }
                }
            }

            // Merging all the temporal indexes
            #pragma omp critical
            {
                for (int i = 0; i < static_cast<int>(temp_indexes.size()); ++i) {
                    indices[i].insert(temp_indexes[i].begin(), temp_indexes[i].end());
                }
            }
        }

        // Count the row sizes
        SizeType nnz = 0;
        for (IndexType i = 0; i < indices.size(); ++i)
            nnz += indices[i].size();

        rA = CompressedMatrixType(indices.size(), indices.size(), nnz);

        double *Avalues = rA.value_data().begin();
        IndexType *Arow_indices = rA.index1_data().begin();
        IndexType *Acol_indices = rA.index2_data().begin();

        // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(rA.size1()); i++)
            Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rA.size1()); ++i) {
            const IndexType row_begin = Arow_indices[i];
            const IndexType row_end = Arow_indices[i + 1];
            IndexType k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); ++it) {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
        }

        rA.set_filled(indices.size() + 1, nnz);

        Timer::Stop("MatrixStructure");
    }

    /**
     * @brief This function is exactly same as the ConstructMatrixStructure() function in base class except that the function has the call to ApplyConstraints function call once the element and conditions compute their equation slave_ids
     * @param pScheme The pointer to the integration scheme
     * @param rT The global relation matrix
     * @param rModelPart The model part to compute
    */
    virtual void ConstructRelationMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& rT,
        ModelPart& rModelPart
        )
    {
        // Filling with zero the matrix (creating the structure)
        Timer::Start("RelationMatrixStructure");

        std::unordered_map<IndexType, IndexType> solvable_dof_reorder;
        std::unordered_map<IndexType, IndexSetType> master_indices;

        // Filling with "ones"
        typedef std::pair<IndexType, IndexType> IndexIndexPairType;
        typedef std::pair<IndexType, IndexSetType> IndexIndexSetPairType;
        IndexType counter = 0;
        for (auto& dof : BaseType::mDofSet) {
            if (dof.EquationId() < BaseType::mEquationSystemSize) {
                const IndexType equation_id = dof.EquationId();
                auto it = mDoFSlaveSet.find(dof);
                if (it == mDoFSlaveSet.end()) {
                    solvable_dof_reorder.insert(IndexIndexPairType(equation_id, counter));
                    master_indices.insert(IndexIndexSetPairType(equation_id, IndexSetType({counter})));
                    ++counter;
                } else {
                    master_indices.insert(IndexIndexSetPairType(equation_id, IndexSetType({})));
                }
            }
        }

        // The process info
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        /// Definition of the eqautio id vector type
        EquationIdVectorType ids(3, 0);
        EquationIdVectorType second_ids(3, 0); // NOTE: Used only on the constraints to take into account the master dofs

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        // TODO: OMP
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;
            it_const->EquationIdVector(ids, second_ids, r_current_process_info);
            for (auto& slave_id : ids) {
                if (slave_id < BaseType::mEquationSystemSize) {
                    auto it_slave = solvable_dof_reorder.find(slave_id);
                    if (it_slave == solvable_dof_reorder.end()) {
                        for (auto& master_id : second_ids) {
                            if (master_id < BaseType::mEquationSystemSize) {
                                auto& master_row_indices = master_indices[slave_id];
                                master_row_indices.insert(solvable_dof_reorder[master_id]);
                            }
                        }
                    }
                }
            }
        }

        KRATOS_DEBUG_ERROR_IF_NOT(BaseType::mEquationSystemSize == master_indices.size()) << "Inconsistency in the dofs size: " << BaseType::mEquationSystemSize << "\t vs \t" << master_indices.size() << std::endl;

        // Count the row sizes
        SizeType nnz = 0;
        for (IndexType i = 0; i < BaseType::mEquationSystemSize; ++i) {
            nnz += master_indices[i].size();
        }

        rT = CompressedMatrixType(BaseType::mEquationSystemSize, mDoFToSolveSystemSize, nnz);

        double *Tvalues = rT.value_data().begin();
        IndexType *Trow_indices = rT.index1_data().begin();
        IndexType *Tcol_indices = rT.index2_data().begin();

        // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Trow_indices[0] = 0;
        for (IndexType i = 0; i < BaseType::mEquationSystemSize; ++i)
            Trow_indices[i + 1] = Trow_indices[i] + master_indices[i].size();

        KRATOS_DEBUG_ERROR_IF_NOT(Trow_indices[BaseType::mEquationSystemSize] == nnz) << "Nonzero values does not coincide with the row index definition: " << Trow_indices[BaseType::mEquationSystemSize] << " vs " << nnz << std::endl;

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rT.size1()); ++i) {
            const IndexType row_begin = Trow_indices[i];
            const IndexType row_end = Trow_indices[i + 1];
            IndexType k = row_begin;
            for (auto it = master_indices[i].begin(); it != master_indices[i].end(); ++it) {
                Tcol_indices[k] = *it;
                Tvalues[k] = 0.0;
                k++;
            }

            master_indices[i].clear(); //deallocating the memory

            std::sort(&Tcol_indices[row_begin], &Tcol_indices[row_end]);
        }

        rT.set_filled(BaseType::mEquationSystemSize + 1, nnz);

        // Setting ones
        for (auto& solv_dof : solvable_dof_reorder) {
            rT(solv_dof.first, solv_dof.second) = 1.0;
        }

        Timer::Stop("RelationMatrixStructure");
    }

    /**
     * @brief This function is exactly same as the Build() function in base class except that the function
     * @details It has the call to ApplyConstraints function call once the LHS or RHS are computed by elements and conditions
     */
    void BuildWithConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY

        // We build the original system
        BaseType::Build(pScheme, rModelPart, rA, rb);

        // Assemble the constraints
        const double start_build = OpenMPUtils::GetCurrentTime();

        // The current process info
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize the constant vector
        double aux_constant_value = 0.0;
        const IndexType nl_iteration_number = r_current_process_info[NL_ITERATION_NUMBER];
        const IndexType step = r_current_process_info[STEP];
        const bool add_constant_vector = (nl_iteration_number == 1 && step == 1) ? true : false;

        if (add_constant_vector) {
            if (mpConstantVector == NULL) { // If the pointer is not initialized initialize it to an empty vector
                TSystemVectorPointerType pNewConstantVector = TSystemVectorPointerType(new TSystemVectorType(0));
                mpConstantVector.swap(pNewConstantVector);
            }
            if (mpConstantVector->size() != BaseType::mEquationSystemSize) {
                mpConstantVector->resize(BaseType::mEquationSystemSize, false);
            }
        }

        // We build the global T matrix and the g constant vector
        TSystemMatrixType& rTMatrix = *mpTMatrix;
        TSystemVectorType& rConstantVector = *mpConstantVector;

        // Filling constant vector
        if (add_constant_vector) {
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(BaseType::mEquationSystemSize); ++i) {
                rConstantVector[i] = 0.0;
            }
        }

        // We compute only once (or if cleared)
        if (mCleared) {
            mCleared = false;

            // Auxiliar set to reorder master DoFs
            std::unordered_map<IndexType, IndexType> solvable_dof_reorder;

            // Filling with "ones"
            typedef std::pair<IndexType, IndexType> IndexIndexPairType;
            IndexType counter = 0;
            for (auto& dof : BaseType::mDofSet) {
                if (dof.EquationId() < BaseType::mEquationSystemSize) {
                    const IndexType equation_id = dof.EquationId();
                    auto it = mDoFSlaveSet.find(dof);
                    if (it == mDoFSlaveSet.end()) {
                        solvable_dof_reorder.insert(IndexIndexPairType(equation_id, counter));
                        ++counter;
                    }
                }
            }

            // Contributions to the system
            LocalSystemMatrixType transformation_matrix = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType constant_vector = LocalSystemVectorType(0);

            // Vector containing the localization in the system of the different terms
            EquationIdVectorType slave_equation_id, master_equation_id;

            const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

            #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_equation_id, master_equation_id)
            {
                #pragma omp for schedule(guided, 512)
                for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                    auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

                    // Detect if the constraint is active or not. If the user did not make any choice the constraint
                    // It is active by default
                    bool constraint_is_active = true;
                    if (it_const->IsDefined(ACTIVE))
                        constraint_is_active = it_const->Is(ACTIVE);

                    if (constraint_is_active) {
                        it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);

                        it_const->EquationIdVector(slave_equation_id, master_equation_id, r_current_process_info);

                        // Reassign reordered dofs to the master side
                        for (auto& id : master_equation_id) {
                            id = solvable_dof_reorder[id];
                        }

                        if (add_constant_vector) {
                            for (IndexType i = 0; i < slave_equation_id.size(); ++i) {
                                const IndexType i_global = slave_equation_id[i];
                                if (i_global < BaseType::mEquationSystemSize) {
                                    double& value = rConstantVector[i_global];
                                    const double constant_value = constant_vector[i];
                                    #pragma omp atomic
                                    value += constant_value;
                                    #pragma omp atomic
                                    aux_constant_value += std::abs(constant_value);
                                }
                            }
                        }

                        // Assemble the constraint contribution
                    #ifdef USE_LOCKS_IN_ASSEMBLY
                        AssembleRelationMatrix(rTMatrix, transformation_matrix, slave_equation_id, master_equation_id, BaseType::mLockArray);
                    #else
                        AssembleRelationMatrix(rTMatrix, transformation_matrix, slave_equation_id, master_equation_id);
                    #endif
                    }
                }
            }
        }

        // We compute the transposed matrix of the global relation matrix
        TSystemMatrixType T_transpose_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
        SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, rTMatrix, 1.0);

        // The proper way to include the constants is in the RHS as T^t(f - A * g)
        VectorType rb_copy(rb);
        if (add_constant_vector && aux_constant_value > std::numeric_limits<double>::epsilon()) {
            VectorType aux_constant_vector(rConstantVector);
            TSparseSpace::Mult(rA, rConstantVector, aux_constant_vector);
            TSparseSpace::UnaliasedAdd(rb_copy, -1.0, aux_constant_vector);
        }

        // The auxiliar matrix to store the intermediate matrix multiplication
        TSystemMatrixType auxiliar_A_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
        SparseMatrixMultiplicationUtility::MatrixMultiplication(T_transpose_matrix, rA, auxiliar_A_matrix);

        // We do a backup of the matrix before apply the constraints
        if (mpOldAMatrix == NULL) { // If the pointer is not initialized initialize it to an empty matrix
            TSystemMatrixPointerType pNewOldAMatrix = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            mpOldAMatrix.swap(pNewOldAMatrix);
        }
        (*mpOldAMatrix).swap(rA);
        // We resize of system of equations
        rA.resize(mDoFToSolveSystemSize, mDoFToSolveSystemSize, false);
        rb.resize(mDoFToSolveSystemSize, false);

        // Final multiplication
        SparseMatrixMultiplicationUtility::MatrixMultiplication(auxiliar_A_matrix, rTMatrix, rA);
        TSparseSpace::Mult(T_transpose_matrix, rb_copy, rb);

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Constraint relation build time and multiplication: " << stop_build - start_build << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building with constraints" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method resize and initializes the system of euqations
     * @details Additionally what is done in the base class the constraints are initialized
     * @param pA The pointer to the LHS matrix
     * @param pDx The pointer to the vector of Unknowns
     * @param pb The pointer to the RHS vector
     * @param rModelPart The model part to be computed
     */
    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
        ) override
    {
        // We resize the basic system
        BaseType::ResizeAndInitializeVectors(pScheme, pA, pDx, pb, rModelPart);

        // Now we resize the relation matrix used on the MPC solution
        if(rModelPart.MasterSlaveConstraints().size() > 0) {
            if (mpTMatrix == NULL) { // If the pointer is not initialized initialize it to an empty matrix
                TSystemMatrixPointerType pNewT = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
                mpTMatrix.swap(pNewT);
            }

            TSystemMatrixType& rTMatrix = *mpTMatrix;

            // Resizing the system vectors and matrix
            if (rTMatrix.size1() == 0 || BaseType::GetReshapeMatrixFlag() || mCleared) { // If the matrix is not initialized
                rTMatrix.resize(BaseType::mEquationSystemSize, mDoFToSolveSystemSize, false);
                ConstructRelationMatrixStructure(pScheme, rTMatrix, rModelPart);
            } else {
                if (rTMatrix.size1() != BaseType::mEquationSystemSize || rTMatrix.size2() != mDoFToSolveSystemSize) {
                    KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                    rTMatrix.resize(BaseType::mEquationSystemSize, mDoFToSolveSystemSize, true);
                    ConstructRelationMatrixStructure(pScheme, rTMatrix, rModelPart);
                }
            }
        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details In the base ResidualBasedEliminationBuilderAndSolver does nothing, due to the fact that the BC are automatically managed with the elimination. But in the constrints approach the slave DoF depending on fixed DoFs must be reconstructed
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
        // Rebuild the DoFs is expensive, we can choose if we do it or not
        if (mRebuildFixedDoFs) {
            // Vector containing the localization in the system of the different terms
            DofsVectorType slave_dof_list, master_dof_list;

            // The dof set to be reset
            typedef std::unordered_set < DofPointerType, DofPointerHasher> set_type;
            set_type reset_slave_dofs;
            reset_slave_dofs.reserve(mDoFSlaveSet.size());

            const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

            // Set to zero the NODAL_VAUX
            const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

            // Auxiliar zero array
            const array_1d<double, 3> zero_array = ZeroVector(3);

            #pragma omp parallel for
            for(int i = 0;  i< num_nodes; ++i) {
                auto it_node = rModelPart.Nodes().begin() + i;
                it_node->SetValue(NODAL_MAUX, 0.0);
                it_node->SetValue(NODAL_VAUX, zero_array);
            }

            // Contributions to the system
            LocalSystemMatrixType transformation_matrix = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType constant_vector = LocalSystemVectorType(0);

            // First we do a database of the slave dofs to be reset
            #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_dof_list, master_dof_list)
            {
                // The current process info
                ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

                set_type aux_reset_slave_dofs;
                aux_reset_slave_dofs.reserve(mDoFSlaveSet.size());

                #pragma omp for schedule(guided, 512)
                for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                    auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

                    // Detect if the constraint is active or not. If the user did not make any choice the constraint
                    // It is active by default
                    bool constraint_is_active = true;
                    if (it_const->IsDefined(ACTIVE))
                        constraint_is_active = it_const->Is(ACTIVE);

                    if (constraint_is_active) {
                        it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                        it_const->GetDofList(slave_dof_list, master_dof_list, r_current_process_info);

                        // The direct assignment will be done in case all the dofs of the constraint are fixed
                        // Otherwise the behaviour is undetermined
                        // How do you impose directly at the same time you solve the system?, good question
                        SizeType number_dofs = 0;
                        for (auto& master_dof : master_dof_list) {
                            if (master_dof->IsFixed()) {
                                ++number_dofs;
                            }
                        }

                        if (number_dofs == master_dof_list.size()) {
                            aux_reset_slave_dofs.insert(slave_dof_list.begin(), slave_dof_list.end());

                            IndexType slave_index = 0;
                            for (auto& slave_dof : slave_dof_list) {
                                const IndexType node_id = slave_dof->Id();
                                const auto& r_variable = slave_dof->GetVariable();
                                auto pnode = rModelPart.pGetNode(node_id);

                                if (r_variable.IsNotComponent()) {
                                    double& coeff = pnode->GetValue(NODAL_MAUX);

                                    for (IndexType master_index = 0; master_index < master_dof_list.size(); ++master_index) {
                                        #pragma omp atomic
                                        coeff += transformation_matrix(slave_index, master_index);
                                    }
                                } else {
                                    IndexType component_index;
                                    const std::string& variable_name = r_variable.Name();
                                    std::size_t found = variable_name.find("_X");
                                    if (found!=std::string::npos) {
                                        component_index = 0;
                                    } else {
                                        found = variable_name.find("_Y");
                                        if (found!=std::string::npos) {
                                            component_index = 1;
                                        } else {
                                            component_index = 2;
                                        }
                                    }

                                    double& coeff = pnode->GetValue(NODAL_VAUX)[component_index];

                                    for (IndexType master_index = 0; master_index < master_dof_list.size(); ++master_index) {
                                        #pragma omp atomic
                                        coeff += transformation_matrix(slave_index, master_index);
                                    }
                                }
                                ++slave_index;
                            }
                        }
                    }
                }

                // We merge all the sets in one thread
                #pragma omp critical
                {
                    reset_slave_dofs.insert(aux_reset_slave_dofs.begin(), aux_reset_slave_dofs.end());
                }
            }

            // We reset the DoFs
            for (auto& it_dof : reset_slave_dofs) {
                it_dof->GetSolutionStepValue() = 0.0;
            }

            // Auxiliar vectors
            LocalSystemVectorType aux_solution_vector = LocalSystemVectorType(0);
            LocalSystemVectorType aux_master_vector = LocalSystemVectorType(0);

            // Now we add the contributions of the fixed master dofs
            #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_dof_list, master_dof_list, aux_solution_vector, aux_master_vector)
            {
                // The current process info
                ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

                #pragma omp for schedule(guided, 512)
                for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                    auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

                    // Detect if the constraint is active or not. If the user did not make any choice the constraint
                    // It is active by default
                    bool constraint_is_active = true;
                    if (it_const->IsDefined(ACTIVE))
                        constraint_is_active = it_const->Is(ACTIVE);

                    if (constraint_is_active) {
                        it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                        it_const->GetDofList(slave_dof_list, master_dof_list, r_current_process_info);

                        // The direct assignment will be done in case all the dofs of the constraint are fixed
                        // Otherwise the behaviour is undetermined
                        // How do you impose directly at the same time you solve the system?, good question
                        SizeType number_dofs = 0;
                        for (auto& master_dof : master_dof_list) {
                            if (master_dof->IsFixed()) {
                                ++number_dofs;
                            }
                        }

                        // Now we add the contribution
                        if (number_dofs == master_dof_list.size()) {
                            aux_solution_vector.resize(slave_dof_list.size());
                            aux_master_vector.resize(number_dofs);

                            IndexType counter = 0;
                            for (auto& master_dof : master_dof_list) {
                                aux_master_vector[counter] = master_dof->GetSolutionStepValue();
                                ++counter;
                            }
                            noalias(aux_solution_vector) = prod(transformation_matrix, aux_master_vector);

                            counter = 0;
                            for (auto& slave_dof : slave_dof_list) {
                                const IndexType node_id = slave_dof->Id();
                                const auto& r_variable = slave_dof->GetVariable();
                                auto pnode = rModelPart.pGetNode(node_id);
                                double coeff;
                                if (r_variable.IsNotComponent()) {
                                    coeff = pnode->GetValue(NODAL_MAUX);
                                } else {
                                    IndexType component_index;
                                    const std::string& variable_name = r_variable.Name();
                                    std::size_t found = variable_name.find("_X");
                                    if (found!=std::string::npos) {
                                        component_index = 0;
                                    } else {
                                        found = variable_name.find("_Y");
                                        if (found!=std::string::npos) {
                                            component_index = 1;
                                        } else {
                                            component_index = 2;
                                        }
                                    }
                                    coeff = pnode->GetValue(NODAL_VAUX)[component_index];
                                }
                                const double inv_coeff = std::abs(coeff) > std::numeric_limits<double>::epsilon()? 1.0/coeff : 1.0;

                                double& aux_value = slave_dof->GetSolutionStepValue();
                                #pragma omp atomic
                                aux_value += inv_coeff * aux_solution_vector[counter];

                                ++counter;
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();

        mDoFSlaveSet = DofsArrayType();

        // Clear constraint system
        TSystemMatrixType& rTMatrix = *mpTMatrix;
        rTMatrix.resize(0, 0, false);
        TSystemVectorType& rConstantVector = *mpConstantVector;
        rConstantVector.resize(0, false);

        // Set the flag
        mCleared = true;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 1) << "Clear Function called" << std::endl;
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
     * @brief This method computes the equivalent coounter part of the SetUpSystem when using constraints
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystemWithConstraints(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // First we set up the system of equations without constraints
        BaseType::SetUpSystem(rModelPart);

        // Add the computation of the global ids of the solvable dofs
        IndexType counter = 0;
        for (auto& dof : BaseType::mDofSet) {
            if (dof.EquationId() < BaseType::mEquationSystemSize) {
                auto it = mDoFSlaveSet.find(dof);
                if (it == mDoFSlaveSet.end()) {
                    ++counter;
                }
            }
        }

        // The total system of equations to be solved
        mDoFToSolveSystemSize = counter;

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::FormulateGlobalMasterSlaveRelations failed ..");
    }

    /**
     * @brief This method reconstructs the slave solution after Solving.
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rA System matrix
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void ReconstructSlaveSolutionAfterSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY

        // We get the global T matrix and the constant vector
        TSystemMatrixType& rTMatrix = *mpTMatrix;
        TSystemVectorType& rConstantVector = *mpConstantVector;

        // The current process info
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // We reconstruct the complete vector of Unknowns
        VectorType Dx_copy(rDx);
        rDx.resize(BaseType::mEquationSystemSize);
        TSparseSpace::Mult(rTMatrix, Dx_copy, rDx);
        // Add the constant vector
        const IndexType nl_iteration_number = r_current_process_info[NL_ITERATION_NUMBER];
        const IndexType step = r_current_process_info[STEP];
        const bool add_constant_vector = (nl_iteration_number == 1 && step == 1) ? true : false;

        if (add_constant_vector) {
            TSparseSpace::UnaliasedAdd(rDx, 1.0, rConstantVector);
        }

        // We reconstruct the LHS
        if (mReassembleLHS) {
            // We compute the transposed matrix of the global relation matrix
            TSystemMatrixType T_transpose_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, rTMatrix, 1.0);

            // The auxiliar matrix to store the intermediate matrix multiplication
            TSystemMatrixType auxiliar_A_matrix(BaseType::mEquationSystemSize, mDoFToSolveSystemSize);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(rA, T_transpose_matrix, auxiliar_A_matrix);

            // We resize the LHS
            rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);

            // Final multiplication
            SparseMatrixMultiplicationUtility::MatrixMultiplication(rTMatrix, auxiliar_A_matrix, rA);
        } else {
            // Simply restore old LHS
            (rA).swap(*mpOldAMatrix);
            mpOldAMatrix = NULL;
        }

        // Reconstruct the RHS
        VectorType rb_copy(rb);
        rb.resize(BaseType::mEquationSystemSize, false);
        TSparseSpace::Mult(rTMatrix, rb_copy, rb);

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::ReconstructSlaveSolutionAfterSolve failed ..");
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

}; /* Class ResidualBasedEliminationBuilderAndSolverWithConstraints */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
