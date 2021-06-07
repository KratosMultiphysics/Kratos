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
#include "utilities/constraint_utilities.h"
#include "input_output/logger.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"

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
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BuilderAndSolverBaseType;

    /// Definition of the base class
    typedef ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The definition of the current class
    typedef ResidualBasedEliminationBuilderAndSolverWithConstraints<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

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

    /// Map definition
    typedef std::unordered_map<IndexType, IndexType> IndexMapType;

    /// MPC definitions
    typedef MasterSlaveConstraint MasterSlaveConstraintType;
    typedef typename MasterSlaveConstraint::Pointer MasterSlaveConstraintPointerType;
    typedef std::vector<IndexType> VectorIndexType;
    typedef Vector VectorType;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedEliminationBuilderAndSolverWithConstraints() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedEliminationBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedEliminationBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        const bool CheckConstraintRelation = true,
        const bool ResetRelationMatrixEachIteration = false
        )
        : BaseType(pNewLinearSystemSolver),
          mCheckConstraintRelation(CheckConstraintRelation),
          mResetRelationMatrixEachIteration(ResetRelationMatrixEachIteration)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedEliminationBuilderAndSolverWithConstraints() override
    {
    }

    /**
     * @brief Create method
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    typename BuilderAndSolverBaseType::Pointer Create(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(pNewLinearSystemSolver,ThisParameters);
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

        if(rModelPart.MasterSlaveConstraints().size() > 0)
            BuildRHSWithConstraints(pScheme, rModelPart, b);
        else
            BaseType::BuildRHS(pScheme, rModelPart, b);

        KRATOS_CATCH("")

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

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "elimination_builder_and_solver_with_constraints",
            "check_constraint_relation"            : true,
            "reset_relation_matrix_each_iteration" : true
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
        return "elimination_builder_and_solver_with_constraints";
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

    TSystemMatrixPointerType mpTMatrix = NULL;             /// This is matrix containing the global relation for the constraints
    TSystemMatrixPointerType mpOldAMatrix = NULL;          /// This is matrix containing the old LHS structure
    TSystemVectorPointerType mpConstantVector = NULL;      /// This is vector containing the rigid movement of the constraint
    TSystemVectorPointerType mpDeltaConstantVector = NULL; /// This is vector contains the effective constant displacement
    DofsArrayType mDoFMasterFixedSet;                      /// The set containing the fixed master DoF of the system
    DofsArrayType mDoFSlaveSet;                            /// The set containing the slave DoF of the system
    SizeType mDoFToSolveSystemSize = 0;                    /// Number of degrees of freedom of the problem to actually be solved
    IndexMapType mReactionEquationIdMap;                   /// In order to know the corresponding EquaionId for each component of the reaction vector

    bool mCheckConstraintRelation = false;                 /// If we do a constraint check relation
    bool mResetRelationMatrixEachIteration = false;        /// If we reset the relation matrix at each iteration

    bool mComputeConstantContribution = false;             /// If we compute the constant contribution of the MPC
    bool mCleared = true;                                  /// If the system has been reseted

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
    */
    void AssembleRelationMatrix(
        TSystemMatrixType& rT,
        const LocalSystemMatrixType& rTransformationMatrix,
        const EquationIdVectorType& rSlaveEquationId,
        const EquationIdVectorType& rMasterEquationId
        )
    {
        const SizeType local_size_1 = rTransformationMatrix.size1();

        for (IndexType i_local = 0; i_local < local_size_1; ++i_local) {
            IndexType i_global = rSlaveEquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize) {
                BaseType::AssembleRowContributionFreeDofs(rT, rTransformationMatrix, i_global, i_local, rMasterEquationId);
            }
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

        // We apply the master/slave relationship before build
        ApplyMasterSlaveRelation(pScheme, rModelPart, rA, rDx, rb);

        // We compute the effective constant vector
        TSystemVectorType dummy_Dx(mDoFToSolveSystemSize);
        TSparseSpace::SetToZero(dummy_Dx);
        ComputeEffectiveConstant(pScheme, rModelPart, dummy_Dx);

        // We do the build (after that we resize the solution vector to avoid problems)
        BuildWithConstraints(pScheme, rModelPart, rA, rb);

        Timer::Stop("Build");

        // Now we apply the BC
        rDx.resize(mDoFToSolveSystemSize, false);
        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) <<
        "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        // We solve the system of equations
        const auto timer = BuiltinTimer();
        const double start_solve = timer.ElapsedSeconds();
        Timer::Start("Solve");
        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = timer.ElapsedSeconds();

        // We compute the effective constant vector
        ComputeEffectiveConstant(pScheme, rModelPart, rDx);

        // We reconstruct the Unknowns vector and the residual
        const double start_reconstruct_slaves = timer.ElapsedSeconds();
        ReconstructSlaveSolutionAfterSolve(pScheme, rModelPart, rA, rDx, rb);

        const double stop_reconstruct_slaves = timer.ElapsedSeconds();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Reconstruct slaves time: " << stop_reconstruct_slaves - start_reconstruct_slaves << std::endl;

        // Some verbosity
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) <<
        "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief The same methods as the base class but with constraints
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rb The RHS vector of the system of equations
     */
    void BuildRHSWithConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        )
    {
        Timer::Start("Build RHS");

        // Resetting to zero the vector of reactions
        if(BaseType::mCalculateReactionsFlag) {
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));
        }

        // Builing without BC
        BuildRHSNoDirichlet(pScheme,rModelPart,rb);

        Timer::Stop("Build RHS");

        ApplyDirichletConditionsRHS(pScheme, rModelPart, rb);

        // We get the global T matrix
        const TSystemMatrixType& rTMatrix = *mpTMatrix;

        // Reconstruct the RHS
        TSystemVectorType rb_copy = rb;
        rb.resize(BaseType::mEquationSystemSize, false);
        TSparseSpace::Mult(rTMatrix, rb_copy, rb);

        // Adding contribution to reactions
        TSystemVectorType& r_reactions_vector = *BaseType::mpReactionsVector;

        if (BaseType::mCalculateReactionsFlag) {
            for (auto& r_dof : BaseType::mDofSet) {
                const bool is_master_fixed = mDoFMasterFixedSet.find(r_dof) == mDoFMasterFixedSet.end() ? false : true;
                const bool is_slave = mDoFSlaveSet.find(r_dof) == mDoFSlaveSet.end() ? false : true;
                if (is_master_fixed || is_slave) { // Fixed or MPC dof
                    const IndexType equation_id = r_dof.EquationId();
                    r_reactions_vector[mReactionEquationIdMap[equation_id]] += rb[equation_id];
                }
            }
        }

        // Some verbosity
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) <<
        "After the solution of the system" << "\nRHS vector = " << rb << std::endl;
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
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

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
        for (auto p_dof : dof_global_set) {
            dof_temp_all.push_back( p_dof );
        }
        dof_temp_all.Sort();
        BaseType::mDofSet = dof_temp_all;

        dof_temp_slave.reserve(dof_global_slave_set.size());
        for (auto p_dof : dof_global_slave_set) {
            dof_temp_slave.push_back( p_dof );
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
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
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

        double norm_b = 0.0;
        if (TSparseSpace::Size(rb) > 0)
            norm_b = TSparseSpace::TwoNorm(rb);

        if (norm_b > 0.0) {
             // Create the auxiliar dof set
             DofsArrayType aux_dof_set;
             aux_dof_set.reserve(mDoFToSolveSystemSize);
             for (auto& r_dof : BaseType::mDofSet) {
                 if (r_dof.EquationId() < BaseType::mEquationSystemSize) {
                    auto it = mDoFSlaveSet.find(r_dof);
                    if (it == mDoFSlaveSet.end())
                        aux_dof_set.push_back( &r_dof );
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
        block_for_each(indices, [](IndexSetType& rIndices){
            rIndices.reserve(40);
        });

        /// Definition of the eqautio id vector type
        EquationIdVectorType ids(3, 0);
        EquationIdVectorType second_ids(3, 0); // NOTE: Used only on the constraints to take into account the master dofs

        #pragma omp parallel firstprivate(ids, second_ids)
        {
            // The process info
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

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
                pScheme->EquationId(*it_elem, ids, r_current_process_info);

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
                pScheme->EquationId(*it_cond, ids, r_current_process_info);
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

                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = true;
                if( it_const->IsDefined(ACTIVE) ) {
                    constraint_is_active = it_const->Is(ACTIVE);
                }

                if(constraint_is_active) {
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

        IndexPartition<std::size_t>(rA.size1()).for_each([&](std::size_t Index){
            const IndexType row_begin = Arow_indices[Index];
            const IndexType row_end = Arow_indices[Index + 1];
            IndexType k = row_begin;
            for (auto it = indices[Index].begin(); it != indices[Index].end(); ++it) {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

            indices[Index].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
        });

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

        IndexMapType solvable_dof_reorder;
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
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        /// Definition of the eqautio id vector type
        EquationIdVectorType ids(3, 0);
        EquationIdVectorType second_ids(3, 0); // NOTE: Used only on the constraints to take into account the master dofs

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
        // TODO: OMP
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = it_const_begin + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = true;
            if( it_const->IsDefined(ACTIVE) ) {
                constraint_is_active = it_const->Is(ACTIVE);
            }

            if(constraint_is_active) {
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

        IndexPartition<std::size_t>(rT.size1()).for_each([&](std::size_t Index){
            const IndexType row_begin = Trow_indices[Index];
            const IndexType row_end = Trow_indices[Index + 1];
            IndexType k = row_begin;
            for (auto it = master_indices[Index].begin(); it != master_indices[Index].end(); ++it) {
                Tcol_indices[k] = *it;
                Tvalues[k] = 0.0;
                k++;
            }

            master_indices[Index].clear(); //deallocating the memory

            std::sort(&Tcol_indices[row_begin], &Tcol_indices[row_end]);
        });
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
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     * @param UseBaseBuild If the abse Build function will be used
     */
    void BuildWithConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb,
        const bool UseBaseBuild = true
        )
    {
        KRATOS_TRY

        // We build the original system
        if (UseBaseBuild)
            BaseType::Build(pScheme, rModelPart, rA, rb);
        else
            BuildWithoutConstraints(pScheme, rModelPart, rA, rb);

        // Assemble the constraints
        const auto timer = BuiltinTimer();

        // We get the global T matrix
        const TSystemMatrixType& rTMatrix = *mpTMatrix;

        // We compute only once (or if cleared)
        if (mCleared) {
            mCleared = false;
            ComputeConstraintContribution(pScheme, rModelPart, true, mComputeConstantContribution);
        } else if (mResetRelationMatrixEachIteration) {
            ResetConstraintSystem();
            ComputeConstraintContribution(pScheme, rModelPart, mResetRelationMatrixEachIteration, mComputeConstantContribution);
        }

        // We compute the transposed matrix of the global relation matrix
        TSystemMatrixType T_transpose_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
        SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, rTMatrix, 1.0);

        // The proper way to include the constants is in the RHS as T^t(f - A * g)
        TSystemVectorType rb_copy = rb;
        if (mComputeConstantContribution) {
            // We get the g constant vector
            TSystemVectorType& rDeltaConstantVector = *mpDeltaConstantVector;
            TSystemVectorType aux_constant_vector(rDeltaConstantVector);
            TSparseSpace::Mult(rA, rDeltaConstantVector, aux_constant_vector);
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

        // Cleaning up memory
        auxiliar_A_matrix.resize(0, 0, false);
        T_transpose_matrix.resize(0, 0, false);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() >= 1) << "Constraint relation build time and multiplication: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2) << "Finished parallel building with constraints" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS of the system
     */
    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY

        // Assemble the constraints
        const auto timer = BuiltinTimer();

        // We get the global T matrix
        const TSystemMatrixType& rTMatrix = *mpTMatrix;

        // We compute only once (or if cleared)
        if (mCleared) {
            mCleared = false;
            ComputeConstraintContribution(pScheme, rModelPart, true, mComputeConstantContribution);
        } else if (mResetRelationMatrixEachIteration) {
            ResetConstraintSystem();
            ComputeConstraintContribution(pScheme, rModelPart, mResetRelationMatrixEachIteration, mComputeConstantContribution);
        }

        // We compute the transposed matrix of the global relation matrix
        TSystemMatrixType T_transpose_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
        SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, rTMatrix, 1.0);

        // We build the original system
        TSystemMatrixType A; // Dummy auxiliar matrix we ned to build anyway because are needed to impose the rigid displacements
        if (mComputeConstantContribution) {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(pScheme, A, rModelPart);
            BuildWithoutConstraints(pScheme, rModelPart, A, rb);
        } else {
            BuildRHSNoDirichletWithoutConstraints(pScheme, rModelPart, rb);
        }

        // The proper way to include the constants is in the RHS as T^t(f - A * g)
        TSystemVectorType rb_copy = rb;
        if (mComputeConstantContribution) {
            // We get the g constant vector
            TSystemVectorType& rDeltaConstantVector = *mpDeltaConstantVector;
            TSystemVectorType aux_constant_vector(rDeltaConstantVector);
            TSparseSpace::Mult(A, rDeltaConstantVector, aux_constant_vector);
            TSparseSpace::UnaliasedAdd(rb_copy, -1.0, aux_constant_vector);
        }

        rb.resize(mDoFToSolveSystemSize, false);

        // Final multiplication
        TSparseSpace::Mult(T_transpose_matrix, rb_copy, rb);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() >= 1) << "Constraint relation build time and multiplication: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2) << "Finished parallel building with constraints" << std::endl;

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

        // If needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag) {
            const SizeType reactions_vector_size = BaseType::mDofSet.size() - mDoFToSolveSystemSize + mDoFMasterFixedSet.size();
            TSystemVectorType& rReactionsVector = *(BaseType::mpReactionsVector);
            if (rReactionsVector.size() != reactions_vector_size)
                rReactionsVector.resize(reactions_vector_size, false);
        }

        // Now we resize the relation matrix used on the MPC solution
        if(rModelPart.MasterSlaveConstraints().size() > 0) {
            if (mpTMatrix == NULL) { // If the pointer is not initialized initialize it to an empty matrix
                TSystemMatrixPointerType pNewT = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
                mpTMatrix.swap(pNewT);
            }

            // The rigid movement
            if (mpConstantVector == NULL) { // If the pointer is not initialized initialize it to an empty vector
                TSystemVectorPointerType pNewConstantVector = TSystemVectorPointerType(new TSystemVectorType(0));
                mpConstantVector.swap(pNewConstantVector);
            }

            // The effective rigid movement
            if (mpDeltaConstantVector == NULL) { // If the pointer is not initialized initialize it to an empty vector
                TSystemVectorPointerType pNewConstantVector = TSystemVectorPointerType(new TSystemVectorType(0));
                mpDeltaConstantVector.swap(pNewConstantVector);
            }

            // System matrices/vectors
            TSystemMatrixType& rTMatrix = *mpTMatrix;
            TSystemVectorType& rConstantVector = *mpConstantVector;
            TSystemVectorType& rDeltaConstantVector = *mpDeltaConstantVector;

            // Resizing the system matrix
            if (rTMatrix.size1() == 0 || BaseType::GetReshapeMatrixFlag() || mCleared) { // If the matrix is not initialized
                rTMatrix.resize(BaseType::mEquationSystemSize, mDoFToSolveSystemSize, false);
                ConstructRelationMatrixStructure(pScheme, rTMatrix, rModelPart);
            } else {
                if (rTMatrix.size1() != BaseType::mEquationSystemSize || rTMatrix.size2() != mDoFToSolveSystemSize) {
                    KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                    rTMatrix.resize(BaseType::mEquationSystemSize, mDoFToSolveSystemSize, false);
                    ConstructRelationMatrixStructure(pScheme, rTMatrix, rModelPart);
                }
            }

            // Resizing the system vector
            // The rigid movement
            if (rConstantVector.size() != BaseType::mEquationSystemSize || BaseType::GetReshapeMatrixFlag() || mCleared) {
                rConstantVector.resize(BaseType::mEquationSystemSize, false);
                mComputeConstantContribution = ComputeConstraintContribution(pScheme, rModelPart);
            } else {
                if (rConstantVector.size() != BaseType::mEquationSystemSize) {
                    KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                    rConstantVector.resize(BaseType::mEquationSystemSize, false);
                    mComputeConstantContribution = ComputeConstraintContribution(pScheme, rModelPart);
                }
            }
            // The effective rigid movement
            if (mComputeConstantContribution) {
                if (rDeltaConstantVector.size() != BaseType::mEquationSystemSize || BaseType::GetReshapeMatrixFlag() || mCleared) {
                    rDeltaConstantVector.resize(BaseType::mEquationSystemSize, false);
                } else {
                    if (rDeltaConstantVector.size() != BaseType::mEquationSystemSize) {
                        KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                        rDeltaConstantVector.resize(BaseType::mEquationSystemSize, false);
                    }
                }
            }
        }
    }

    /**
     * @brief It computes the reactions of the system
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
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
        KRATOS_TRY

        // Refresh RHS to have the correct reactions
        BuildRHS(pScheme, rModelPart, rb);

        // Adding contribution to reactions
        TSystemVectorType& r_reactions_vector = *BaseType::mpReactionsVector;

        // Updating variables
        for (auto& r_dof : BaseType::mDofSet) {
            if ((r_dof.IsFixed()) || mDoFSlaveSet.find(r_dof) != mDoFSlaveSet.end()) {
                r_dof.GetSolutionStepReactionValue() = -r_reactions_vector[mReactionEquationIdMap[r_dof.EquationId()]];
            }
        }

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::CalculateReactions failed ..");
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
        KRATOS_TRY;

        if (mDoFMasterFixedSet.size() > 0) {
            // We apply the same method as in the block builder and solver but instead of fixing the fixed Dofs, we just fix the master fixed Dofs
            std::vector<double> scaling_factors (mDoFToSolveSystemSize, 0.0);

            // NOTE: Dofs are assumed to be numbered consecutively
            const auto it_dof_begin = BaseType::mDofSet.begin();
            IndexType counter = 0;
            for (IndexType i = 0; i < BaseType::mDofSet.size(); ++i) {
                auto it_dof = it_dof_begin + i;
                const IndexType equation_id = it_dof->EquationId();
                if (equation_id < BaseType::mEquationSystemSize ) {
                    auto it_first_check = mDoFSlaveSet.find(*it_dof);
                    if (it_first_check == mDoFSlaveSet.end()) {
                        auto it_second_check = mDoFSlaveSet.find(*it_dof);
                        if (it_second_check == mDoFSlaveSet.end()) {
                            if(mDoFMasterFixedSet.find(*it_dof) == mDoFMasterFixedSet.end()) {
                                scaling_factors[counter] = 1.0;
                            }
                        }
                        counter += 1;
                    }
                }
            }

            double* Avalues = rA.value_data().begin();
            IndexType* Arow_indices = rA.index1_data().begin();
            IndexType* Acol_indices = rA.index2_data().begin();

            // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
            #pragma omp parallel for
            for(int k = 0; k < static_cast<int>(mDoFToSolveSystemSize); ++k) {
                const IndexType col_begin = Arow_indices[k];
                const IndexType col_end = Arow_indices[k+1];
                bool empty = true;
                for (IndexType j = col_begin; j < col_end; ++j) {
                    if(Avalues[j] != 0.0) {
                        empty = false;
                        break;
                    }
                }

                if(empty) {
                    rA(k,k) = 1.0;
                    rb[k] = 0.0;
                }
            }

            IndexPartition<std::size_t>(mDoFToSolveSystemSize).for_each([&](std::size_t Index){
                const IndexType col_begin = Arow_indices[Index];
                const IndexType col_end = Arow_indices[Index+1];
                const double k_factor = scaling_factors[Index];
                if (k_factor == 0) {
                    // Zero out the whole row, except the diagonal
                    for (IndexType j = col_begin; j < col_end; ++j)
                        if (Acol_indices[j] != Index )
                            Avalues[j] = 0.0;

                    // Zero out the RHS
                    rb[Index] = 0.0;
                } else {
                    // Zero out the column which is associated with the zero'ed row
                    for (IndexType j = col_begin; j < col_end; ++j) {
                        if(scaling_factors[ Acol_indices[j] ] == 0 ) {
                            Avalues[j] = 0.0;
                        }
                    }
                }
            });
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();

        // Reseting auxiliar set of dofs
        mDoFMasterFixedSet = DofsArrayType();
        mDoFSlaveSet = DofsArrayType();

        // Clearing the relation map
        mReactionEquationIdMap.clear();

        // Clear constraint system
        if (mpTMatrix != nullptr)
            TSparseSpace::Clear(mpTMatrix);
        if (mpConstantVector != nullptr)
            TSparseSpace::Clear(mpConstantVector);
        if (mpDeltaConstantVector != nullptr)
            TSparseSpace::Clear(mpDeltaConstantVector);

        // Set the flag
        mCleared = true;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 1) << "Clear Function called" << std::endl;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
        mCheckConstraintRelation = ThisParameters["check_constraint_relation"].GetBool();
        mResetRelationMatrixEachIteration = ThisParameters["reset_relation_matrix_each_iteration"].GetBool();
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
        // Set equation id for degrees of freedom the free degrees of freedom are positioned at the beginning of the system, while the fixed one are at the end (in opposite order).
        //
        // That means that if the EquationId is greater than "mEquationSystemSize" the pointed degree of freedom is restrained
        // This is almost the same SetUpSystem from ResidualBasedEliminationBuilderAndSolver, but we don't discard from the system the fixed dofs that are part of a constraint at the same time

        /// First we detect the master fixed DoFs ///

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Vector containing the localization in the system of the different terms
        DofsVectorType slave_dof_list, master_dof_list;

        // Declaring temporal variables
        DofsArrayType dof_temp_fixed_master;

        typedef std::unordered_set < DofPointerType, DofPointerHasher> set_type;
        set_type dof_global_fixed_master_set;

        // Iterate over constraints
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
        #pragma omp parallel firstprivate(slave_dof_list, master_dof_list)
        {
            // We cleate the temporal set and we reserve some space on them
            set_type dof_temp_fixed_master_set;
            dof_temp_fixed_master_set.reserve(2000);

            #pragma omp for schedule(guided, 512) nowait
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = it_const_begin + i_const;

                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = true;
                if (it_const->IsDefined(ACTIVE))
                    constraint_is_active = it_const->Is(ACTIVE);

                if (constraint_is_active) {
                    it_const->GetDofList(slave_dof_list, master_dof_list, r_current_process_info);

                    // Filling the set of dofs master and fixed at the same time
                    for (auto& master_dof : master_dof_list) {
                        if (master_dof->IsFixed()) {
                            dof_temp_fixed_master_set.insert(master_dof);
                        }
                    }
                }
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                dof_global_fixed_master_set.insert(dof_temp_fixed_master_set.begin(), dof_temp_fixed_master_set.end());
            }
        }

        dof_temp_fixed_master.reserve(dof_global_fixed_master_set.size());
        for (auto p_dof : dof_global_fixed_master_set) {
            dof_temp_fixed_master.push_back( p_dof );
        }
        dof_temp_fixed_master.Sort();
        mDoFMasterFixedSet = dof_temp_fixed_master;

        /// Now we compute as expected ///
        int free_id = 0;
        int fix_id = BaseType::mDofSet.size();

        for (auto& dof : BaseType::mDofSet) {
            if (dof.IsFixed()) {
                auto it = mDoFMasterFixedSet.find(dof);
                if (it == mDoFMasterFixedSet.end()) {
                    dof.SetEquationId(--fix_id);
                } else {
                    dof.SetEquationId(free_id++);
                }
            } else {
                dof.SetEquationId(free_id++);
            }
        }

        BaseType::mEquationSystemSize = fix_id;

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

        // Finally we build the relation between the EquationID and the component of the reaction
        counter = 0;
        for (auto& r_dof : BaseType::mDofSet) {
            const bool is_master_fixed = mDoFMasterFixedSet.find(r_dof) == mDoFMasterFixedSet.end() ? false : true;
            const bool is_slave = mDoFSlaveSet.find(r_dof) == mDoFSlaveSet.end() ? false : true;
            if (is_master_fixed || is_slave) { // Fixed or MPC dof
                mReactionEquationIdMap.insert({r_dof.EquationId(), counter});
                ++counter;
            }
        }

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::SetUpSystemWithConstraints failed ..");
    }

    /**
     * @brief This method initializes the DoF using the master/slave relationship
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
     * @param rb The RHS vector of the system of equations
     */
    void ApplyMasterSlaveRelation(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY

        // First we reset the slave dofs
        ConstraintUtilities::ResetSlaveDofs(rModelPart);

        // Now we apply the constraints
        ConstraintUtilities::ApplyConstraints(rModelPart);

        KRATOS_CATCH("");
    }

    /**
     * @brief This method checks that the master/slave relation is properly set
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rDx The vector of unkowns
     * @param rDxSolved The vector of unkowns actually solved
     */
    bool CheckMasterSlaveRelation(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rDx,
        TSystemVectorType& rDxSolved
        )
    {
        KRATOS_TRY

        // Auxiliar values
        const auto it_dof_begin = BaseType::mDofSet.begin();
        TSystemVectorType current_solution(mDoFToSolveSystemSize);
        TSystemVectorType updated_solution(BaseType::mEquationSystemSize);
        TSystemVectorType residual_solution(BaseType::mEquationSystemSize);

        // Get current values
        IndexType counter = 0;
        for (IndexType i = 0; i < BaseType::mDofSet.size(); ++i) {
            auto it_dof = it_dof_begin + i;
            const IndexType equation_id = it_dof->EquationId();
            if (equation_id < BaseType::mEquationSystemSize ) {
                auto it = mDoFSlaveSet.find(*it_dof);
                if (it == mDoFSlaveSet.end()) {
                    current_solution[counter] = it_dof->GetSolutionStepValue() + rDxSolved[counter];
                    counter += 1;
                }
            }
        }

        block_for_each(BaseType::mDofSet, [&, this](Dof<double>& rDof){
            const IndexType equation_id = rDof.EquationId();
            if (equation_id < this->mEquationSystemSize ) {
                residual_solution[equation_id] = rDof.GetSolutionStepValue() + rDx[equation_id];
            }
        });

        // Apply master slave constraints
        const TSystemMatrixType& rTMatrix = *mpTMatrix;
        TSparseSpace::Mult(rTMatrix, current_solution, updated_solution);

        if (mComputeConstantContribution) {
            ComputeConstraintContribution(pScheme, rModelPart, false, true);
            const TSystemVectorType& rConstantVector = *mpConstantVector;
            TSparseSpace::UnaliasedAdd(updated_solution, 1.0, rConstantVector);
        }

        TSparseSpace::UnaliasedAdd(residual_solution, -1.0, updated_solution);

        // Check database
        for(int k = 0; k < static_cast<int>(BaseType::mEquationSystemSize); ++k) {
            if (std::abs(residual_solution[k]) > std::numeric_limits<double>::epsilon()) return false;
        }

        return true;

        KRATOS_CATCH("");
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
        const TSystemMatrixType& rTMatrix = *mpTMatrix;

        // We reconstruct the complete vector of Unknowns
        TSystemVectorType Dx_copy = rDx;
        rDx.resize(BaseType::mEquationSystemSize);
        TSparseSpace::Mult(rTMatrix, Dx_copy, rDx);

        // Add the constant vector
        if (mComputeConstantContribution) {
            const TSystemVectorType& rDeltaConstantVector = *mpDeltaConstantVector;
            TSparseSpace::UnaliasedAdd(rDx, 1.0, rDeltaConstantVector);
        }

        // We check the solution
        if (mCheckConstraintRelation) {
            KRATOS_ERROR_IF_NOT(CheckMasterSlaveRelation(pScheme, rModelPart, rDx, Dx_copy)) << "The relation between master/slave dofs is not respected" << std::endl;
        }

        // Simply restore old LHS
        (rA).swap(*mpOldAMatrix);
        mpOldAMatrix = NULL;

        // Reconstruct the RHS
        TSystemVectorType rb_copy = rb;
        rb.resize(BaseType::mEquationSystemSize, false);
        TSparseSpace::Mult(rTMatrix, rb_copy, rb);

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::ReconstructSlaveSolutionAfterSolve failed ..");
    }

    /**
     * @brief Function to perform the build the system without constraints
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void BuildWithoutConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        )
    {
        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Getting the array of elements
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditons_array = rModelPart.Conditions();

        // Contributions to the system
        LocalSystemMatrixType lhs_contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType rhs_contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_id;

        // Assemble all elements and conditions
        #pragma omp parallel firstprivate( lhs_contribution, rhs_contribution, equation_id)
        {
            // Elements
            const auto it_elem_begin = r_elements_array.begin();
            const int nelements = static_cast<int>(r_elements_array.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i<nelements; ++i) {
                auto it_elem = it_elem_begin + i;
                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool element_is_active = true;
                if (it_elem->IsDefined(ACTIVE))
                    element_is_active = it_elem->Is(ACTIVE);

                if (element_is_active) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_elem, lhs_contribution, rhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleWithoutConstraints(rA, rb, lhs_contribution, rhs_contribution, equation_id);
                }
            }

            // Conditions
            const auto it_cond_begin = r_conditons_array.begin();
            const int nconditions = static_cast<int>(r_conditons_array.size());
            #pragma omp  for schedule(guided, 512)
            for (int i = 0; i<nconditions; ++i) {
                auto it_cond = it_cond_begin + i;
                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool condition_is_active = true;
                if (it_cond->IsDefined(ACTIVE))
                    condition_is_active = it_cond->Is(ACTIVE);

                if (condition_is_active) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_cond, lhs_contribution, rhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleWithoutConstraints(rA, rb, lhs_contribution, rhs_contribution, equation_id);
                }
            }
        }
    }

    /**
     * @brief Function to perform the build of the RHS without constraints
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS of the system
     */
    void BuildRHSNoDirichletWithoutConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        )
    {
        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Getting the array of elements
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditons_array = rModelPart.Conditions();

        // Contributions to the system
        LocalSystemVectorType rhs_contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_id;

        // Assemble all elements and conditions
        #pragma omp parallel firstprivate( rhs_contribution, equation_id)
        {
            // Elements
            const auto it_elem_begin = r_elements_array.begin();
            const int nelements = static_cast<int>(r_elements_array.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i<nelements; ++i) {
                auto it_elem = it_elem_begin + i;
                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool element_is_active = true;
                if (it_elem->IsDefined(ACTIVE))
                    element_is_active = it_elem->Is(ACTIVE);

                if (element_is_active) {
                    // Calculate elemental Right Hand Side Contribution
                    pScheme->CalculateRHSContribution(*it_elem, rhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleRHSWithoutConstraints(rb, rhs_contribution, equation_id);
                }
            }

            // Conditions
            const auto it_cond_begin = r_conditons_array.begin();
            const int nconditions = static_cast<int>(r_conditons_array.size());
            #pragma omp  for schedule(guided, 512)
            for (int i = 0; i<nconditions; ++i) {
                auto it_cond = it_cond_begin + i;
                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool condition_is_active = true;
                if (it_cond->IsDefined(ACTIVE))
                    condition_is_active = it_cond->Is(ACTIVE);

                if (condition_is_active) {
                    // Calculate elemental contribution
                    pScheme->CalculateRHSContribution(*it_cond, rhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleRHSWithoutConstraints(rb, rhs_contribution, equation_id);
                }
            }
        }
    }

    /**
    * @brief This function does the assembling of the LHS and RHS
    * @note The main difference respect the block builder and solver is the fact that the fixed DoFs are not considered on the assembling
    */
    void AssembleWithoutConstraints(
        TSystemMatrixType& rA,
        TSystemVectorType& rb,
        const LocalSystemMatrixType& rLHSContribution,
        const LocalSystemVectorType& rRHSContribution,
        const Element::EquationIdVectorType& rEquationId
        )
    {
        const SizeType local_size = rLHSContribution.size1();

        // Assemble RHS
        AssembleRHSWithoutConstraints(rb, rRHSContribution, rEquationId);

        // Assemble LHS
        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            const IndexType i_global = rEquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize) {
                BaseType::AssembleRowContributionFreeDofs(rA, rLHSContribution, i_global, i_local, rEquationId);
            }
        }
    }


    /**
     * @brief Assembling local contribution of nodes and elements in the RHS
     * @param rb The RHS vector
     */
    void AssembleRHSWithoutConstraints(
        TSystemVectorType& rb,
        const LocalSystemVectorType& rRHSContribution,
        const Element::EquationIdVectorType& rEquationId
        )
    {
        const SizeType local_size = rRHSContribution.size();

        if (!BaseType::mCalculateReactionsFlag) {
            for (IndexType i_local = 0; i_local < local_size; ++i_local) {
                const IndexType i_global = rEquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) { // free dof
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& r_b_value = rb[i_global];
                    const double rhs_value = rRHSContribution[i_local];

                    AtomicAdd(r_b_value, rhs_value);
                }
            }
        } else {
            TSystemVectorType& r_reactions_vector = *BaseType::mpReactionsVector;
            for (IndexType i_local = 0; i_local < local_size; ++i_local) {
                const IndexType i_global = rEquationId[i_local];
                auto it_dof = BaseType::mDofSet.begin() + i_global;

                const bool is_master_fixed = mDoFMasterFixedSet.find(*it_dof) == mDoFMasterFixedSet.end() ? false : true;
                const bool is_slave = mDoFSlaveSet.find(*it_dof) == mDoFSlaveSet.end() ? false : true;
                if (is_master_fixed || is_slave) { // Fixed or MPC dof
                    double& r_b_value = r_reactions_vector[mReactionEquationIdMap[i_global]];
                    const double rhs_value = rRHSContribution[i_local];

                    AtomicAdd(r_b_value, rhs_value);
                } else if (it_dof->IsFree()) {  // Free dof not in the MPC
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& r_b_value = rb[i_global];
                    const double& rhs_value = rRHSContribution[i_local];

                    AtomicAdd(r_b_value, rhs_value);
                }
            }
        }
    }

    /**
     * @brief This method set to zero the relation matrix
     */
    void ResetConstraintSystem()
    {
        TSystemMatrixType& rTMatrix = *mpTMatrix;
        double *Tvalues = rTMatrix.value_data().begin();

        IndexPartition<std::size_t>(rTMatrix.nnz()).for_each([&Tvalues](std::size_t Index){
            Tvalues[Index] = 0.0;
        });

        IndexMapType solvable_dof_reorder;

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

        // Setting ones
        for (auto& solv_dof : solvable_dof_reorder) {
            rTMatrix(solv_dof.first, solv_dof.second) = 1.0;
        }

        if (mComputeConstantContribution) {
            TSystemVectorType& rConstantVector = *mpConstantVector;
            TSparseSpace::SetToZero(rConstantVector);
        }
    }

    /**
     * @brief This method applies the BC, only in the RHS
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rb The RHS vector of the system of equations
     */
    void ApplyDirichletConditionsRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY;

        if (mDoFMasterFixedSet.size() > 0) {
            // NOTE: dofs are assumed to be numbered consecutively
            const auto it_dof_begin = BaseType::mDofSet.begin();

            IndexPartition<std::size_t>(mDoFToSolveSystemSize).for_each([&, this](std::size_t Index){
                auto it_dof = it_dof_begin + Index;
                if (Index < this->mEquationSystemSize) {
                    auto it = mDoFSlaveSet.find(*it_dof);
                    if (it == mDoFSlaveSet.end()) {
                        if(mDoFMasterFixedSet.find(*it_dof) != mDoFMasterFixedSet.end()) {
                            rb[Index] = 0.0;
                        }
                    }
                }
            });
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief This method computes the absolute constant contribution of the MPC
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param ComputeTranslationMatrix If the translation matrix will be assembled
     * @param ComputeConstantVector If the constant vector will be assembled
     * @return If there are constant constraints
     */
    bool ComputeConstraintContribution(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        const bool ComputeTranslationMatrix = false,
        const bool ComputeConstantVector = false
        )
    {
        KRATOS_TRY;

        // We build the global T matrix and the g constant vector
        TSystemMatrixType& rTMatrix = *mpTMatrix;
        TSystemVectorType& rConstantVector = *mpConstantVector;

        // Filling constant vector
        if (ComputeConstantVector) {
            IndexPartition<std::size_t>(this->mEquationSystemSize).for_each([&rConstantVector](std::size_t Index){
                rConstantVector[Index] = 0.0;
            });
        }

        // Auxiliar set to reorder master DoFs
        IndexMapType solvable_dof_reorder;

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

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize the constant vector
        double aux_constant_value = 0.0;

        // Contributions to the system
        LocalSystemMatrixType transformation_matrix = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType constant_vector = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType slave_equation_id, master_equation_id;

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        std::unordered_set<IndexType> auxiliar_constant_equations_ids;

        #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_equation_id, master_equation_id)
        {
            std::unordered_set<IndexType> auxiliar_temp_constant_equations_ids;
            auxiliar_temp_constant_equations_ids.reserve(2000);

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

                    if (ComputeConstantVector) {
                        for (IndexType i = 0; i < slave_equation_id.size(); ++i) {
                            const IndexType i_global = slave_equation_id[i];
                            if (i_global < BaseType::mEquationSystemSize) {
                                const double constant_value = constant_vector[i];
                                if (std::abs(constant_value) > 0.0) {
                                    auxiliar_temp_constant_equations_ids.insert(i_global);
                                    double& r_value = rConstantVector[i_global];
                                    AtomicAdd(r_value, constant_value);
                                }
                            }
                        }
                    } else {
                        for (IndexType i = 0; i < slave_equation_id.size(); ++i) {
                            const IndexType i_global = slave_equation_id[i];
                            if (i_global < BaseType::mEquationSystemSize) {
                                const double constant_value = constant_vector[i];
                                AtomicAdd(aux_constant_value, std::abs(constant_value));
                            }
                        }
                    }

                    if (ComputeTranslationMatrix) {
                        // Assemble the constraint contribution
                        AssembleRelationMatrix(rTMatrix, transformation_matrix, slave_equation_id, master_equation_id);
                    }
                }
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                auxiliar_constant_equations_ids.insert(auxiliar_temp_constant_equations_ids.begin(), auxiliar_temp_constant_equations_ids.end());
            }
        }

        return aux_constant_value > std::numeric_limits<double>::epsilon() ? true : false;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method computes the efective constant
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rDxSolved The vector of unkowns actually solved
     */
    void ComputeEffectiveConstant(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rDxSolved
        )
    {
        if (mComputeConstantContribution) {
            // We get
            const TSystemMatrixType& rTMatrix = *mpTMatrix;
            TSystemVectorType& rConstantVector = *mpConstantVector;
            TSystemVectorType& rDeltaConstantVector = *mpDeltaConstantVector;
            TSparseSpace::Copy(rConstantVector, rDeltaConstantVector);

            // We reconstruct the complete vector of Unknowns
            TSystemVectorType Dx(BaseType::mEquationSystemSize);
            TSparseSpace::Mult(rTMatrix, rDxSolved, Dx);

            // Compute the effective constant vector
            // Auxiliar initial dof iterator
            const auto it_dof_begin = BaseType::mDofSet.begin();

            TSystemVectorType u(BaseType::mEquationSystemSize);

            block_for_each(BaseType::mDofSet, [&, this](Dof<double>& rDof){
                const IndexType equation_id = rDof.EquationId();
                if (equation_id < this->mEquationSystemSize ) {
                    u[equation_id] = rDof.GetSolutionStepValue() + Dx[equation_id];
                }
            });

            TSystemVectorType u_bar(mDoFToSolveSystemSize);
            IndexType counter = 0;
            for (IndexType i = 0; i < BaseType::mDofSet.size(); ++i) {
                auto it_dof = it_dof_begin + i;
                const IndexType equation_id = it_dof->EquationId();
                if (equation_id < BaseType::mEquationSystemSize ) {

                    auto it = mDoFSlaveSet.find(*it_dof);
                    if (it == mDoFSlaveSet.end()) {
                        u_bar[counter] = it_dof->GetSolutionStepValue() + rDxSolved[counter];
                        counter += 1;
                    }
                }
            }
            TSystemVectorType u_bar_complete(BaseType::mEquationSystemSize);
            TSparseSpace::Mult(rTMatrix, u_bar, u_bar_complete);

            TSparseSpace::UnaliasedAdd(rDeltaConstantVector, 1.0, u_bar_complete);
            TSparseSpace::UnaliasedAdd(rDeltaConstantVector, -1.0, u);
        }
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
