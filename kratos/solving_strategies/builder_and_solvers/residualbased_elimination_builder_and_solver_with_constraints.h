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
// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#endif
// #define USE_LOCKS_IN_ASSEMBLY

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "includes/master_slave_constraint.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/helper_classes_for_constraint_builder.h"

#include "containers/pointer_vector_map.h"
#include "containers/pointer_hash_map_set.h"
#include "containers/data_value_container.h"
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

    /// Definition of the classes from the base class
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
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
    typedef typename BaseType::ElementsContainerType ElementsContainerType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;
    typedef typename BaseType::DofsVectorType DofsVectorType;
    typedef typename BaseType::CompressedMatrixType CompressedMatrixType;

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
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedEliminationBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
//         ) : BaseType(pNewLinearSystemSolver, ThisParameters)
    {
    }

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedEliminationBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
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
        KRATOS_TRY;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the dofs" << std::endl;

        DofsVectorType dof_list, auxiliar_dof_list;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        const SizeType nthreads = OpenMPUtils::GetNumThreads();

    #ifdef USE_GOOGLE_HASH
        typedef google::dense_hash_set < DofPointerType, DofPointerHasher>  set_type;
    #else
        typedef std::unordered_set < DofPointerType, DofPointerHasher>  set_type;
    #endif

        // Declaring temporal variables
        DofsArrayType dof_temp_all, dof_temp_solvable, dof_temp_master;
        BaseType::mDofSet = DofsArrayType();
        mDoFToSolveSet = DofsArrayType();
        mMasterDoFSet = DofsArrayType();

        std::vector<set_type> dofs_aux_list_all(nthreads);
        std::vector<set_type> dofs_aux_list_master(nthreads);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Number of threads" << nthreads << "\n" << std::endl;

        for (int i = 0; i < static_cast<int>(nthreads); i++) {
    #ifdef USE_GOOGLE_HASH
            dofs_aux_list_all[i].set_empty_key(DofPointerType());
            dofs_aux_list_master[i].set_empty_key(DofPointerType());
    #else
            dofs_aux_list_all[i].reserve(rModelPart.Elements().size());
            dofs_aux_list_master[i].reserve(rModelPart.MasterSlaveConstraints().size());
    #endif
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        #pragma omp parallel firstprivate(dof_list)
        {
            // Gets the array of elements from the modeler
            ElementsArrayType& r_elements_array = rModelPart.Elements();
            const int number_of_elements = static_cast<int>(r_elements_array.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i) {
                auto it_elem = r_elements_array.begin() + i;
                const unsigned int this_thread_id = OpenMPUtils::ThisThread();

                // Gets list of Dof involved on every element
                pScheme->GetElementalDofList(*(it_elem.base()), dof_list, r_current_process_info);

                dofs_aux_list_all[this_thread_id].insert(dof_list.begin(), dof_list.end());
            }
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing condition loop" << std::endl;

        #pragma omp parallel firstprivate(dof_list)
        {
            // Gets the array of conditions from the modeler
            ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
            const int number_of_conditions = static_cast<int>(r_conditions_array.size());
            #pragma omp for  schedule(guided, 512)
            for (int i = 0; i < number_of_conditions; ++i) {
                auto it_cond = r_conditions_array.begin() + i;
                const unsigned int this_thread_id = OpenMPUtils::ThisThread();

                // Gets list of Dof involved on every element
                pScheme->GetConditionDofList(*(it_cond.base()), dof_list, r_current_process_info);
                dofs_aux_list_all[this_thread_id].insert(dof_list.begin(), dof_list.end());
            }
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing constraints loop" << std::endl;

        #pragma omp parallel firstprivate(dof_list, auxiliar_dof_list)
        {
            auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
            const int number_of_constraints = static_cast<int>(r_constraints_array.size());
            #pragma omp for  schedule(guided, 512)
            for (int i = 0; i < number_of_constraints; ++i) {
                auto it_const = r_constraints_array.begin() + i;
                const IndexType this_thread_id = OpenMPUtils::ThisThread();

                // Gets list of Dof involved on every element
                it_const->GetDofList(dof_list, auxiliar_dof_list, r_current_process_info);
                dofs_aux_list_all[this_thread_id].insert(dof_list.begin(), dof_list.end());
                dofs_aux_list_all[this_thread_id].insert(auxiliar_dof_list.begin(), auxiliar_dof_list.end());
                dofs_aux_list_master[this_thread_id].insert(auxiliar_dof_list.begin(), auxiliar_dof_list.end()); // We add the master dofs
            }
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing tree reduction\n" << std::endl;

        // Here we do a reduction in a tree so to have everything on thread 0
        // TODO: move to a function
        /* For all DoF */
        IndexType old_max = nthreads;
        IndexType new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max>=1 && new_max != old_max) {
            // Just for debugging
            KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2) << "old_max" << old_max << " new_max:" << new_max << std::endl;
            for (int i = 0; i < static_cast<int>(new_max); ++i) {
                KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2 && (i + new_max < old_max)) << i << " - " << i+new_max << std::endl;
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); ++i) {
                if (i + new_max < old_max) {
                    dofs_aux_list_all[i].insert(dofs_aux_list_all[i+new_max].begin(), dofs_aux_list_all[i+new_max].end());
                    dofs_aux_list_all[i+new_max].clear();
                }
            }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        /* For master DoF */
        old_max = nthreads;
        new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max>=1 && new_max != old_max) {
            // Just for debugging
            KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2) << "old_max" << old_max << " new_max:" << new_max << std::endl;
            for (int i = 0; i < static_cast<int>(new_max); ++i) {
                KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2 && (i + new_max < old_max)) << i << " - " << i+new_max << std::endl;
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); ++i) {
                if (i + new_max < old_max) {
                    dofs_aux_list_master[i].insert(dofs_aux_list_master[i+new_max].begin(), dofs_aux_list_master[i+new_max].end());
                    dofs_aux_list_master[i+new_max].clear();
                }
            }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        dof_temp_all.reserve(dofs_aux_list_all[0].size());
        for (auto& dof : dofs_aux_list_all[0]) {
            dof_temp_all.push_back( dof.get() );
        }
        dof_temp_all.Sort();
        BaseType::mDofSet = dof_temp_all;

        dof_temp_master.reserve(dofs_aux_list_master[0].size());
        for (auto& dof : dofs_aux_list_master[0]) {
            dof_temp_master.push_back( dof.get() );
        }
        dof_temp_master.Sort();
        mMasterDoFSet = dof_temp_master;

        // Repeating the loop for solvable DoF (in serial)
        set_type dofs_aux_list_solvable;
        // Gets the array of elements from the modeler
        ElementsArrayType& r_elements_array = rModelPart.Elements();
        const int number_of_elements = static_cast<int>(r_elements_array.size());
        #pragma omp for schedule(guided, 512) nowait
        for (int i = 0; i < number_of_elements; ++i) {
            auto it_elem = r_elements_array.begin() + i;

            // Gets list of Dof involved on every element
            pScheme->GetElementalDofList(*(it_elem.base()), dof_list, r_current_process_info);
            dofs_aux_list_solvable.insert(dof_list.begin(), dof_list.end());
        }

        // Gets the array of conditions from the modeler
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
        const int number_of_conditions = static_cast<int>(r_conditions_array.size());
        #pragma omp for  schedule(guided, 512)
        for (int i = 0; i < number_of_conditions; ++i) {
            auto it_cond = r_conditions_array.begin() + i;

            // Gets list of Dof involved on every element
            pScheme->GetConditionDofList(*(it_cond.base()), dof_list, r_current_process_info);
            dofs_aux_list_solvable.insert(dof_list.begin(), dof_list.end());
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing constraints loop" << std::endl;

        // We don't do a loop in OMP because erase is not a parallel threadsafe operation
        typedef std::pair<IndexType, IndexSetType> IndexIndexSetPairType;
        auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
        const int number_of_constraints = static_cast<int>(r_constraints_array.size());

        for (int i = 0; i < number_of_constraints; ++i) {
            auto it_const = r_constraints_array.begin() + i;

            // Gets list of Dof involved on every element
            it_const->GetDofList(dof_list, auxiliar_dof_list, r_current_process_info);
            // We remove the slave dofs
            for (auto& dof : dof_list) {
                auto it = dofs_aux_list_solvable.find(dof);
                // Check if Iterator is valid
                if(it != dofs_aux_list_solvable.end())
                    dofs_aux_list_solvable.erase(it);

                // We add the master dofs to the map of slave dofs relation
                const IndexType equation_id_dof = dof->EquationId();
                auto it_relation_dof = mSlaveMasterDoFRelation.find(equation_id_dof);
                if ( it_relation_dof == mSlaveMasterDoFRelation.end()) {
                    IndexSetType dummy_set;
                    dummy_set.reserve(auxiliar_dof_list.size());
                    for (auto& auxiliar_dof : auxiliar_dof_list)
                        dummy_set.insert(auxiliar_dof->EquationId());
                    mSlaveMasterDoFRelation.insert(IndexIndexSetPairType(equation_id_dof, dummy_set));
                } else {
                    IndexSetType& set = it_relation_dof->second;
                    for (auto& auxiliar_dof : auxiliar_dof_list)
                        set.insert(auxiliar_dof->EquationId());
                }
            }
            dofs_aux_list_solvable.insert(auxiliar_dof_list.begin(), auxiliar_dof_list.end()); // We add the master dofs
        }

        dof_temp_solvable.reserve(dofs_aux_list_solvable.size());
        for (auto& dof : dofs_aux_list_solvable) {
            dof_temp_solvable.push_back( dof.get() );
        }
        dof_temp_solvable.Sort();
        mDoFToSolveSet = dof_temp_solvable;

        // Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_ERROR_IF( mDoFToSolveSet.size() == 0) << "No degrees of freedom to solve!" << std::endl;
        KRATOS_WARNING_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", mMasterDoFSet.size() == 0) << "No MPC degrees of freedom to solve!" << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing lock array" << std::endl;

#ifdef _OPENMP
        if (BaseType::mLockArray.size() != 0) {
            for (int i = 0; i < static_cast<int>(BaseType::mLockArray.size()); ++i) {
                omp_destroy_lock(&BaseType::mLockArray[i]);
            }
        }
        BaseType::mLockArray.resize(BaseType::mDofSet.size());

        for (int i = 0; i < static_cast<int>(BaseType::mLockArray.size()); ++i) {
            omp_init_lock(&BaseType::mLockArray[i]);
        }
#endif

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "End of setup dof set\n" << std::endl;

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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
    * @brief This method assembles the global relation matrix
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
        SizeType local_size = rTransformationMatrix.size1();

        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            IndexType i_global = rSlaveEquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize) {
            #ifdef USE_LOCKS_IN_ASSEMBLY
                omp_set_lock(&rLockArray[i_global]);
            #endif
                for (IndexType j_local = 0; j_local < local_size; ++j_local) {
                    IndexType j_global = mSolvableDoFReorder[rMasterEquationId[j_local]];
                    if (j_global < BaseType::mEquationSystemSize) {
                        rT(i_global, j_global) += rTransformationMatrix(i_local, j_local);
                    }
                }
            #ifdef USE_LOCKS_IN_ASSEMBLY
                omp_unset_lock(&rLockArray[i_global]);
            #endif

            }
            //note that assembly on fixed rows is not performed here
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

        Build(pScheme, rModelPart, rA, rb);
        rDx.resize(mDoFToSolveSystemSize, false);

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) <<
        "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");
        this->SystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();

        const double start_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        ReconstructSlaveSolutionAfterSolve(rModelPart, rA, rDx, rb);

        const double stop_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Reconstruct slaves time: " << stop_reconstruct_slaves - start_reconstruct_slaves << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) <<
        "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

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

        // Getting the elements from the model
        const int number_of_elements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int number_of_conditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        auto el_begin = rModelPart.ElementsBegin();
        auto cond_begin = rModelPart.ConditionsBegin();

        const SizeType equation_size = BaseType::mEquationSystemSize;

    #ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<IndexType>> indices(equation_size);
        const SizeType empty_key = 2 * equation_size + 10;
    #else
        std::vector<IndexSetType> indices(equation_size);
    #endif

        #pragma omp parallel for firstprivate(equation_size)
        for (int index = 0; index < static_cast<int>(equation_size); ++index) {
        #ifdef USE_GOOGLE_HASH
            indices[index].set_empty_key(empty_key);
        #else
            indices[index].reserve(40);
        #endif
        }

        /// Definition of the eqautio id vector type
        EquationIdVectorType ids(3, 0);

        #pragma omp parallel for firstprivate(number_of_elements, ids)
        for (int i_elem = 0; i_elem<number_of_elements; i_elem++) {
            auto it_elem = el_begin + i_elem;
            pScheme->EquationId( *(it_elem.base()), ids, r_current_process_info);

            for (IndexType i = 0; i < ids.size(); ++i) {
                if (ids[i] < BaseType::mEquationSystemSize) {
                #ifdef _OPENMP
                    omp_set_lock(&BaseType::mLockArray[ids[i]]);
                #endif
                    auto& row_indices = indices[ids[i]];
                    for (auto it = ids.begin(); it != ids.end(); ++it) {
                        if (*it < BaseType::mEquationSystemSize)
                            row_indices.insert(*it);
                    }
                #ifdef _OPENMP
                    omp_unset_lock(&BaseType::mLockArray[ids[i]]);
                #endif
                }
            }

        }

        #pragma omp parallel for firstprivate(number_of_conditions, ids)
        for (int i_cond = 0; i_cond<number_of_conditions; ++i_cond) {
            auto it_cond = cond_begin + i_cond;
            pScheme->Condition_EquationId( *(it_cond.base()), ids, r_current_process_info);
            for (IndexType i = 0; i < ids.size(); ++i) {
                if (ids[i] < BaseType::mEquationSystemSize) {
                #ifdef _OPENMP
                    omp_set_lock(&BaseType::mLockArray[ids[i]]);
                #endif
                    auto& row_indices = indices[ids[i]];
                    for (auto it = ids.begin(); it != ids.end(); ++it) {
                        if (*it < BaseType::mEquationSystemSize)
                            row_indices.insert(*it);
                    }
                #ifdef _OPENMP
                    omp_unset_lock(&BaseType::mLockArray[ids[i]]);
                #endif
                }
            }
        }

        EquationIdVectorType aux_ids(3, 0);

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        #pragma omp parallel for firstprivate(number_of_constraints, ids, aux_ids)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;
            it_const->EquationIdVector(ids, aux_ids, r_current_process_info);
            for (IndexType i = 0; i < ids.size(); ++i) {
                if (ids[i] < BaseType::mEquationSystemSize) {
                #ifdef _OPENMP
                    omp_set_lock(&BaseType::mLockArray[ids[i]]);
                #endif
                    auto &row_indices = indices[ids[i]];
                    row_indices.insert(ids.begin(), ids.end());
                #ifdef _OPENMP
                    omp_unset_lock(&BaseType::mLockArray[ids[i]]);
                #endif
                }
            }
            for (IndexType i = 0; i < aux_ids.size(); ++i) {
                if (aux_ids[i] < BaseType::mEquationSystemSize) {
                #ifdef _OPENMP
                    omp_set_lock(&BaseType::mLockArray[aux_ids[i]]);
                #endif
                    auto &row_indices = indices[aux_ids[i]];
                    row_indices.insert(aux_ids.begin(), aux_ids.end());
                #ifdef _OPENMP
                    omp_unset_lock(&BaseType::mLockArray[aux_ids[i]]);
                #endif
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

        std::unordered_map<IndexType, bool> row_dof_indices;
        row_dof_indices.reserve(BaseType::mEquationSystemSize);

        std::unordered_map<IndexType, IndexSetType> master_indices;
        master_indices.reserve(mMasterDoFSystemSize);

        // We do a pair to know which DoF are pure master MPC
        typedef std::pair<IndexType, bool> PairIdBoolType;

        // The set containing the solvable DoF that are not a master DoF
        IndexSetType dummy_set;
        typedef std::pair<IndexType, IndexSetType> IndexIndexSetPairType;
        for (auto& dof : BaseType::mDofSet) {
            const IndexType equation_id = dof.EquationId();
            if (equation_id < BaseType::mEquationSystemSize) {
                auto it = mMasterDoFSet.find(dof);
                if (it != mMasterDoFSet.end()) {
                    row_dof_indices.insert(PairIdBoolType(equation_id, false));
                    master_indices.insert(IndexIndexSetPairType(equation_id, dummy_set));
                } else {
                    row_dof_indices.insert(PairIdBoolType(equation_id, true));
                }
            }
        }

        // Reserve on the indexes set
        for (auto& master_set : master_indices) {
            (master_set.second).reserve(3);
        }

        // Defining counter for later use
        IndexType counter = 0;

    #ifdef KRATOS_DEBUG
        // Checking that there is a consistency
        for (auto& to_solve : row_dof_indices) {
            if (!(to_solve.second)) {
                ++counter;
            }
        }
        KRATOS_ERROR_IF_NOT(row_dof_indices.size() == BaseType::mEquationSystemSize) << "Inconsistency in the dofs size: " << row_dof_indices.size() << "\t vs \t" << BaseType::mEquationSystemSize << std::endl;
        KRATOS_ERROR_IF_NOT(counter == mMasterDoFSystemSize) << "Inconsistency in the pure master MPC dofs: " << counter << "\t vs \t" << mMasterDoFSystemSize << std::endl;
    #endif

        // The process info
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        /// Definition of the eqautio id vector type
        EquationIdVectorType ids(3, 0);
        EquationIdVectorType aux_ids(3, 0);

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        // TODO: OMP
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;
            it_const->EquationIdVector(ids, aux_ids, r_current_process_info);
            for (auto& slave_id : ids) {
                if (slave_id < BaseType::mEquationSystemSize) {
                    auto it_slave = mSlaveMasterDoFRelation.find(slave_id);
                    if (it_slave != mSlaveMasterDoFRelation.end()) {
                        auto& master_set = it_slave->second;
                        for (auto& master_id : aux_ids) {
                            auto it_master = master_set.find(master_id);
                            if (it_master != master_set.end()) {
                                auto &master_row_indices = master_indices[master_id];
                                if (master_id < BaseType::mEquationSystemSize) {
                                    master_row_indices.insert(mSolvableDoFReorder[master_id]);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Count the row sizes
        SizeType nnz = 0;
        for (auto& to_solve : row_dof_indices) {
            if (to_solve.second) {
                ++nnz;
            } else {
                nnz += master_indices[to_solve.first].size();
            }
        }

        rT = CompressedMatrixType(BaseType::mEquationSystemSize, mDoFToSolveSystemSize, nnz);

        double *Tvalues = rT.value_data().begin();
        IndexType *Trow_indices = rT.index1_data().begin();
        IndexType *Tcol_indices = rT.index2_data().begin();

        // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Trow_indices[0] = 0;
        counter = 0;
        for (auto& to_solve : row_dof_indices) {
            if (to_solve.second) {
                Trow_indices[counter + 1] = Trow_indices[counter] + 1;
            } else {
                Trow_indices[counter + 1] = Trow_indices[counter] + master_indices[to_solve.first].size();
            }
            ++counter;
        }

        KRATOS_DEBUG_ERROR_IF_NOT(Trow_indices[counter] == nnz) << "Nonzero values does not coincide with the row index definition" << std::endl;

        counter = 0;
        // TODO: OMP
        for (auto& to_solve : row_dof_indices) {
            const IndexType row_begin = Trow_indices[counter];
            const IndexType row_end = Trow_indices[counter + 1];

            IndexType k = row_begin;
            if (to_solve.second) {
                Tcol_indices[k] = mSolvableDoFReorder[to_solve.first];
                Tvalues[k] = 1.0;
            } else {
                for (auto& index : master_indices[to_solve.first]) {
                    Tcol_indices[k] = index;
                    Tvalues[k] = 0.0;
                    k++;
                }
                master_indices[to_solve.first].clear(); //deallocating the memory
            }

            std::sort(&Tcol_indices[row_begin], &Tcol_indices[row_end]);

            ++counter;
        }

        rT.set_filled(BaseType::mEquationSystemSize + 1, nnz);

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

        // We build the global T matrix
        TSystemMatrixType& rTMatrix = *mpTMatrix;

        // We compute only once (or if cleared)
        if (BaseType::GetReshapeMatrixFlag() || mCleared) {
            mCleared = false;

            // Contributions to the system
            LocalSystemMatrixType transformation_matrix = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType constant_vector = LocalSystemVectorType(0);

            // Vector containing the localization in the system of the different terms
            EquationIdVectorType slave_equation_id, master_equation_id;

            // The current process info
            ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
            #pragma omp parallel for firstprivate(number_of_constraints, transformation_matrix, constant_vector, slave_equation_id)
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

                    // Assemble the constraint contribution
                #ifdef USE_LOCKS_IN_ASSEMBLY
                    AssembleRelationMatrix(rTMatrix, transformation_matrix, slave_equation_id, master_equation_id, mLockArray);
                #else
                    AssembleRelationMatrix(rTMatrix, transformation_matrix, slave_equation_id, master_equation_id);
                #endif
                }
            }
        }

        // We compute the transposed matrix of the global relation matrix
        TSystemMatrixType T_transpose_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
        SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, rTMatrix, 1.0);

        // The auxiliar matrix to store the intermediate matrix multiplication
        TSystemMatrixType auxiliar_A_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
        SparseMatrixMultiplicationUtility::MatrixMultiplication(T_transpose_matrix, rA, auxiliar_A_matrix);

        // We resize of system of equations
        rA.resize(mDoFToSolveSystemSize, mDoFToSolveSystemSize, false);
        VectorType rb_copy(rb);
        rb.resize(mDoFToSolveSystemSize, false);

        // Final multiplication
        SparseMatrixMultiplicationUtility::MatrixMultiplication(auxiliar_A_matrix, rTMatrix, rA);
        TSparseSpace::Mult(T_transpose_matrix, rb_copy, rb);

        // TODO: The proper way to include the constants os in the RHS as T^t(f - A * c) -> Include the constants once the problem of the delta is solved

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

        // Now we resize the constraint system
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
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();

        mDoFToSolveSet = DofsArrayType();
        mMasterDoFSet = DofsArrayType();
        mSolvableDoFReorder.clear();
        mSlaveMasterDoFRelation.clear();

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

    TSystemMatrixPointerType mpTMatrix = NULL; /// This is matrix containing the global relation for the constraints

    DofsArrayType mDoFToSolveSet;              /// The set containing the solvable DoF of the system
    DofsArrayType mMasterDoFSet;               /// The set containing the master DoF of the system

    SizeType mDoFToSolveSystemSize = 0;        /// Number of degrees of freedom of the problem to actually be solved
    SizeType mMasterDoFSystemSize = 0;         /// Number of master degrees of freedom of the problem

    bool mCleared = true; /// If the system has been reseted

    std::unordered_map<IndexType, IndexType> mSolvableDoFReorder; /// The correlation between the global total DoF order and the solvable DoF order
    std::unordered_map<IndexType, IndexSetType> mSlaveMasterDoFRelation; /// The relation between the slave and the master DoFs

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

        /// A pair of Ids
        typedef std::pair<IndexType, IndexType> IdPairType;

        // Add the computation of the global ids of the solvable dofs
        mSolvableDoFReorder.reserve(mDoFToSolveSystemSize);
        IndexType counter = 0;
        IndexType master_counter = 0;
        for (auto& dof : BaseType::mDofSet) {
            if (dof.EquationId() < BaseType::mEquationSystemSize) {
                auto it = mDoFToSolveSet.find(dof);
                if (it != mDoFToSolveSet.end()) {
                    mSolvableDoFReorder.insert(IdPairType(dof.EquationId(), counter));
                    ++counter;

                    // We check if this dofs is in the master set too
                    auto it_master = mMasterDoFSet.find(dof);
                    if (it_master != mMasterDoFSet.end()) {
                        ++master_counter;
                    }
                }
            }
        }

        // The total system of equations to be solved
        mDoFToSolveSystemSize = mSolvableDoFReorder.size();
        mMasterDoFSystemSize = master_counter;

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::FormulateGlobalMasterSlaveRelations failed ..");
    }

    /**
     * @brief This method reconstructs the slave solution after Solving.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rA System matrix
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void ReconstructSlaveSolutionAfterSolve(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY

        // We get the global T matrix
        TSystemMatrixType& rTMatrix = *mpTMatrix;

        // We reconstruct the complete vector of Unknowns
        VectorType Dx_copy(rDx);
        rDx.resize(BaseType::mEquationSystemSize);
        TSparseSpace::Mult(rTMatrix, Dx_copy, rDx);

        /* We reconstruct the rest of the system  */
        // We compute the transposed matrix of the global relation matrix
        TSystemMatrixType T_transpose_matrix(mDoFToSolveSystemSize, BaseType::mEquationSystemSize);
        SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, rTMatrix, 1.0);

        // The auxiliar matrix to store the intermediate matrix multiplication
        TSystemMatrixType auxiliar_A_matrix(BaseType::mEquationSystemSize, mDoFToSolveSystemSize);
        SparseMatrixMultiplicationUtility::MatrixMultiplication(rA, T_transpose_matrix, auxiliar_A_matrix);

        // We resize of system of equations
        rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
        VectorType rb_copy(rb);
        rb.resize(BaseType::mEquationSystemSize, false);

        // Final multiplication
        SparseMatrixMultiplicationUtility::MatrixMultiplication(rTMatrix, auxiliar_A_matrix, rA);
        TSparseSpace::Mult(rTMatrix, rb_copy, rb);

        // TODO: Remember that the constants are compued in the building side too
        // TODO: Should be applied only once, because it is increment of the Unknowns (to solve later)
//         const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
//
//         // Getting the beginning iterator
//         const auto it_constraints_begin = rModelPart.MasterSlaveConstraints().begin();
//
//         // Contributions to the system
//         LocalSystemMatrixType transformation_matrix = LocalSystemMatrixType(0, 0);
//         LocalSystemVectorType constant_vector = LocalSystemVectorType(0);
//
//         // Vector containing the localization in the system of the different terms
//         EquationIdVectorType slave_equation_id, master_equation_id;
//
//         // The current process info
//         ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
//
//         #pragma omp parallel for schedule(guided, 512) firstprivate(slave_equation_id, master_equation_id, transformation_matrix, constant_vector)
//         for (int i_constraints = 0; i_constraints < number_of_constraints; ++i_constraints) {
//             auto it_const = it_constraints_begin;
//             std::advance(it_const, i_constraints);
//
//             // Get the equation Ids of the constraint
//             it_const->EquationIdVector(slave_equation_id, master_equation_id, r_current_process_info);
//
//             // Get constraint's T and b matrices
//             it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
//
//             IndexType counter = 0;
//             for (auto id : slave_equation_id) {
//                 rDx[id] += constant_vector[counter];
//                 ++counter;
//             }
//         }

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
