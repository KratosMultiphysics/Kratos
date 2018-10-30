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
        BaseType::SetUpSystem(rModelPart);
        if(rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            FormulateGlobalMasterSlaveRelations(rModelPart);
        }
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
        if(mGlobalMasterSlaveConstraints.size() > 0)
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
        if(mGlobalMasterSlaveConstraints.size() > 0)
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

        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = rModelPart.Elements();
        const int nelements = static_cast<int>(pElements.size());

        Element::DofsVectorType ElementalDofList, AuxiliarDofList;

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        unsigned int nthreads = OpenMPUtils::GetNumThreads();

#ifdef USE_GOOGLE_HASH
        typedef google::dense_hash_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#else
        typedef std::unordered_set < typename NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#endif

        std::vector<set_type> dofs_aux_list(nthreads);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Number of threads" << nthreads << "\n" << std::endl;

        for (int i = 0; i < static_cast<int>(nthreads); i++) {
#ifdef USE_GOOGLE_HASH
            dofs_aux_list[i].set_empty_key(NodeType::DofType::Pointer());
#else
//             dofs_aux_list[i] = set_type( allocators[i]);
            dofs_aux_list[i].reserve(nelements);
#endif
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        #pragma omp parallel firstprivate(nelements, ElementalDofList, AuxiliarDofList)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < nelements; i++)
            {
                auto it = pElements.begin() + i;
                const unsigned int this_thread_id = OpenMPUtils::ThisThread();

                // gets list of Dof involved on every element
                pScheme->GetElementalDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);

                dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
            }

            KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing condition loop" << std::endl;

            ConditionsArrayType& pConditions = rModelPart.Conditions();
            const int nconditions = static_cast<int>(pConditions.size());
            #pragma omp for  schedule(guided, 512)
            for (int i = 0; i < nconditions; i++)
            {
                auto it = pConditions.begin() + i;
                const unsigned int this_thread_id = OpenMPUtils::ThisThread();

                // gets list of Dof involved on every element
                pScheme->GetConditionDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);
                dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());

            }

            auto& pConstraints = rModelPart.MasterSlaveConstraints();
            const int nconstraints = static_cast<int>(pConstraints.size());
            #pragma omp for  schedule(guided, 512)
            for (int i = 0; i < nconstraints; i++)
            {
                auto it = pConstraints.begin() + i;
                const unsigned int this_thread_id = OpenMPUtils::ThisThread();

                // gets list of Dof involved on every element
                it->GetDofList(ElementalDofList, AuxiliarDofList, CurrentProcessInfo);
                dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
                dofs_aux_list[this_thread_id].insert(AuxiliarDofList.begin(), AuxiliarDofList.end());

            }
        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing tree reduction\n" << std::endl;

        // Here we do a reduction in a tree so to have everything on thread 0
        unsigned int old_max = nthreads;
        unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max>=1 && new_max != old_max) {
            // Just for debugging
            KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2) << "old_max" << old_max << " new_max:" << new_max << std::endl;
            for (int i = 0; i < static_cast<int>(new_max); ++i) {
                KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", this->GetEchoLevel() > 2 && (i + new_max < old_max)) << i << " - " << i+new_max << std::endl;
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); ++i) {
                if (i + new_max < old_max) {
                    dofs_aux_list[i].insert(dofs_aux_list[i+new_max].begin(), dofs_aux_list[i+new_max].end());
                    dofs_aux_list[i+new_max].clear();
                }
            }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));

        }

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", ( this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dofs_aux_list[0].size());
        for (auto it= dofs_aux_list[0].begin(); it!= dofs_aux_list[0].end(); ++it) {
            Doftemp.push_back( it->get() );
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

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
        if(mGlobalMasterSlaveConstraints.size() > 0)
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
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        )
    {
        KRATOS_TRY

        const double start_update_constraints = OpenMPUtils::GetCurrentTime();
        this->UpdateConstraintsForBuilding(rModelPart);
        const double stop_update_constraints = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Constraints update time : " << stop_update_constraints - start_update_constraints << std::endl;

        Timer::Start("Build");

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) << "Before the solution of the system"
                                                                                                         << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");
        this->SystemSolveWithPhysics(A, Dx, b, rModelPart);
        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();

        const double start_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        ReconstructSlaveSolutionAfterSolve(rModelPart, A, Dx, b);
        const double stop_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Reconstruct slaves time: " << stop_reconstruct_slaves - start_reconstruct_slaves << std::endl;


        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) << "After the solution of the system"
                                                                                                         << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

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
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        auto el_begin = rModelPart.ElementsBegin();
        auto cond_begin = rModelPart.ConditionsBegin();

        const SizeType equation_size = BaseType::mEquationSystemSize;

    #ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<IndexType>> indices(equation_size);
        const SizeType empty_key = 2 * equation_size + 10;
    #else
        std::vector<std::unordered_set<IndexType>> indices(equation_size);
    #endif

        #pragma omp parallel for firstprivate(equation_size)
        for (int index = 0; index < static_cast<int>(equation_size); ++index) {
        #ifdef USE_GOOGLE_HASH
            indices[index].set_empty_key(empty_key);
        #else
            indices[index].reserve(40);
        #endif
        }

        // The process info
        auto r_process_info = rModelPart.GetProcessInfo();

        /// Definition of the eqautio id vector type
        EquationIdVectorType ids(3, 0);

        #pragma omp parallel for firstprivate(nelements, ids)
        for (int i_elem = 0; i_elem<nelements; i_elem++) {
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

        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int i_cond = 0; i_cond<nconditions; ++i_cond) {
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

        const int nconstraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        #pragma omp parallel for firstprivate(nconstraints, ids, aux_ids)
        for (int i_const = 0; i_const < nconstraints; i_const++) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;
            it_const->EquationIdVector(ids, aux_ids, r_process_info);
            for (IndexType i = 0; i < ids.size(); ++i) {
            #ifdef _OPENMP
                omp_set_lock(&BaseType::mLockArray[ids[i]]);
            #endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
            #ifdef _OPENMP
                omp_unset_lock(&BaseType::mLockArray[ids[i]]);
            #endif
            }
            for (IndexType i = 0; i < aux_ids.size(); ++i) {
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

        const SizeType equation_size = BaseType::mEquationSystemSize;
        const SizeType slave_equation_size = mSlaveDofsSystemSize;

    #ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<IndexType>> master_indices(equation_size);
        std::vector<google::dense_hash_set<IndexType>> slave_indices(slave_equation_size);
        const SizeType empty_key = 2 * equation_size + 10;
    #else
        std::vector<std::unordered_set<IndexType>> master_indices(equation_size);
        std::vector<std::unordered_set<IndexType>> slave_indices(slave_equation_size);
    #endif

        #pragma omp parallel for firstprivate(equation_size)
        for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
        {
        #ifdef USE_GOOGLE_HASH
            master_indices[iii].set_empty_key(empty_key);
        #else
            master_indices[iii].reserve(40);
        #endif
        }

        #pragma omp parallel for firstprivate(slave_equation_size)
        for (int iii = 0; iii < static_cast<int>(slave_equation_size); iii++)
        {
        #ifdef USE_GOOGLE_HASH
            slave_indices[iii].set_empty_key(empty_key);
        #else
            slave_indices[iii].reserve(40);
        #endif
        }

        EquationIdVectorType slave_ids(3, 0), master_ids(3, 0);

        const int nconstraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        #pragma omp parallel for firstprivate(nconstraints, slave_ids, master_ids)
        for (int i_constraint = 0; i_constraint < nconstraints; ++i_constraint) {
            auto it_constraint = rModelPart.MasterSlaveConstraints().begin() + i_constraint;
            it_constraint->EquationIdVector(slave_ids, master_ids, rModelPart.GetProcessInfo());
            for (IndexType i = 0; i < slave_ids.size(); ++i)
            {
            #ifdef _OPENMP
                omp_set_lock(&BaseType::mLockArray[slave_ids[i]]);
            #endif
                auto& row_indices = master_indices[slave_ids[i]];
                row_indices.insert(slave_ids.begin(), slave_ids.end());
            #ifdef _OPENMP
                omp_unset_lock(&BaseType::mLockArray[slave_ids[i]]);
            #endif
            }
            for (IndexType i = 0; i < master_ids.size(); ++i) {
            #ifdef _OPENMP
                omp_set_lock(&BaseType::mLockArray[master_ids[i]]);
            #endif
                auto& row_indices = slave_indices[master_ids[i]];
                row_indices.insert(master_ids.begin(), master_ids.end());
            #ifdef _OPENMP
                omp_unset_lock(&BaseType::mLockArray[master_ids[i]]);
            #endif
            }
        }

        // Count the row sizes
        SizeType nnz = 0;
        for (IndexType i = 0; i < master_indices.size(); ++i)
            nnz += master_indices[i].size();

        rT = CompressedMatrixType(master_indices.size(), slave_indices.size(), nnz);

        double *Avalues = rT.value_data().begin();
        SizeType *Arow_indices = rT.index1_data().begin();
        SizeType *Acol_indices = rT.index2_data().begin();

        // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (IndexType i = 0; i < rT.size1(); i++)
            Arow_indices[i + 1] = Arow_indices[i] + master_indices[i].size();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rT.size1()); i++) {
            const IndexType row_begin = Arow_indices[i];
            const IndexType row_end = Arow_indices[i + 1];
            IndexType k = row_begin;
            for (auto it = master_indices[i].begin(); it != master_indices[i].end(); ++it) {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                ++k;
            }

            master_indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin],& Acol_indices[row_end]);
        }

        rT.set_filled(master_indices.size() + 1, nnz);

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

        // We build the global T matrix

//         KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;
//         ConstraintImposerType constraint_imposer(mGlobalMasterSlaveConstraints);
//
//         // Getting the elements from the model
//         const int nelements = static_cast<int>(rModelPart.Elements().size());
//
//         // Getting the array of the conditions
//         const int nconditions = static_cast<int>(rModelPart.Conditions().size());
//
//         ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
//         const auto el_begin = rModelPart.ElementsBegin();
//         const auto cond_begin = rModelPart.ConditionsBegin();
//
//         //contributions to the system
//         LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
//         LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
//
//         // Vector containing the localization in the system of the different terms
//         EquationIdVectorType equation_id;
//
//         // Assemble all elements
//         const double start_build = OpenMPUtils::GetCurrentTime();
//
//         #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, equation_id, constraint_imposer)
//         {
//             #pragma omp for schedule(guided, 512) nowait
//             for (int k = 0; k < nelements; ++k) {
//                 auto it = el_begin + k;
//
//                 //detect if the element is active or not. If the user did not make any choice the element
//                 //is active by default
//                 bool element_is_active = true;
//                 if ((it)->IsDefined(ACTIVE))
//                     element_is_active = (it)->Is(ACTIVE);
//
//                 if (element_is_active)
//                 {
//                     //calculate elemental contribution
//                     pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);
//                     constraint_imposer.template ApplyConstraints<Element>(*it, LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);
//
//                     // Assemble the elemental contribution
//                 #ifdef USE_LOCKS_IN_ASSEMBLY
//                     this->Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id, BaseType::mLockArray);
//                 #else
//                     this->Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id);
//                 #endif
//                     // Clean local elemental memory
//                     pScheme->CleanMemory(*(it.base()));
//                 }
//             }
//
//
//             #pragma omp for schedule(guided, 512)
//             for (int k = 0; k < nconditions; ++k) {
//                 auto it = cond_begin + k;
//
//                 //detect if the element is active or not. If the user did not make any choice the element
//                 //is active by default
//                 bool condition_is_active = true;
//                 if ((it)->IsDefined(ACTIVE))
//                     condition_is_active = (it)->Is(ACTIVE);
//
//                 if (condition_is_active) {
//                     // Calculate elemental contribution
//                     pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);
//                     constraint_imposer.template ApplyConstraints<Condition>(*it, LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);
//
//                     // Assemble the elemental contribution
//                 #ifdef USE_LOCKS_IN_ASSEMBLY
//                     this->Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id, BaseType::mLockArray);
//                 #else
//                     this->Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id);
//                 #endif
//
//                     // Clean local elemental memory
//                     pScheme->CleanMemory(*(it.base()));
//                 }
//             }
//         }

//         const double stop_build = OpenMPUtils::GetCurrentTime();
//         KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

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
        if(mGlobalMasterSlaveConstraints.size() > 0) {
            if (mpTMatrix == NULL) // If the pointer is not initialized initialize it to an empty matrix
            {
                TSystemMatrixPointerType pNewT = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
                mpTMatrix.swap(pNewT);
            }

            TSystemMatrixType& rTMatrix = *mpTMatrix;

            // Resizing the system vectors and matrix
            if (rTMatrix.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) { // If the matrix is not initialized
                rTMatrix.resize(BaseType::mEquationSystemSize, mSlaveDofsSystemSize, false);
                ConstructRelationMatrixStructure(pScheme, rTMatrix, rModelPart);
            } else {
                if (rTMatrix.size1() != BaseType::mEquationSystemSize || rTMatrix.size2() != mSlaveDofsSystemSize) {
                    KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                    rTMatrix.resize(BaseType::mEquationSystemSize, mSlaveDofsSystemSize, true);
                    ConstructRelationMatrixStructure(pScheme, rTMatrix, rModelPart);
                }
            }
        }
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

    // This is the set of condensed global constraints.
    SizeType mSlaveDofsSystemSize = 0;        /// The number of master DoF defined. These are the total number of degrees of freedom to be used
    TSystemMatrixPointerType mpTMatrix = NULL; /// This is matrix containing the global relation for the constraints
    GlobalMasterSlaveRelationContainerType mGlobalMasterSlaveConstraints; /// This can be changed to more efficient implementation later on. // TODO: Maybe to remove later

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method condenses the MasterSlaveConstraints which are added on the rModelPart into objects of AuxilaryGlobalMasterSlaveRelation. One unique object for each unique slave.
     * @details These will be used in the ApplyConstraints functions later on.
     * @param rModelPart The model part of the problem to solve
     */
    void FormulateGlobalMasterSlaveRelations(ModelPart& rModelPart)
    {
        KRATOS_TRY

        const double start_formulate = OpenMPUtils::GetCurrentTime();
        // First delete the existing ones
        mGlobalMasterSlaveConstraints.clear();
        // Getting the array of the conditions
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        // Getting the beginning iterator
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        #pragma omp parallel for schedule(guided, 512)
        for (int i_constraints = 0; i_constraints < number_of_constraints; ++i_constraints) {
            auto it = constraints_begin;
            std::advance(it, i_constraints);
            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool constraint_is_active = true;
            if ((it)->IsDefined(ACTIVE))
                constraint_is_active = (it)->Is(ACTIVE);

            if (constraint_is_active) {
                // Assemble the Constraint contribution
                #pragma omp critical
                AssembleConstraint(*it, r_current_process_info);
            }
        }
        const double stop_formulate = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Formulate global constraints time: " << stop_formulate - start_formulate << std::endl;

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::FormulateGlobalMasterSlaveRelations failed ..");
    }


    /**
     * @brief This method assembles the given master slave constraint to the auxiliary global master slave constraints
     * @param rMasterSlaveConstraint object of the master slave constraint to be assembled.
     * @param rCurrentProcessInfo current process info.
     */
    void AssembleConstraint(
        ModelPart::MasterSlaveConstraintType& rMasterSlaveConstraint,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY
        int slave_count = 0;
        LocalSystemMatrixType relation_matrix(0,0);
        LocalSystemVectorType constant_vector(0);
        EquationIdVectorType slave_equation_ids(0);
        EquationIdVectorType master_equation_ids(0);

        //get the equation Ids of the constraint
        rMasterSlaveConstraint.EquationIdVector(slave_equation_ids, master_equation_ids, rCurrentProcessInfo);
        //calculate constraint's T and b matrices
        rMasterSlaveConstraint.CalculateLocalSystem(relation_matrix, constant_vector, rCurrentProcessInfo);

        for (auto slave_equation_id : slave_equation_ids) {
            int master_count = 0;
            auto global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            if (global_constraint == mGlobalMasterSlaveConstraints.end()) {
                mGlobalMasterSlaveConstraints[slave_equation_id] = Kratos::make_unique<AuxiliaryGlobalMasterSlaveConstraintType>(slave_equation_id);
            }

            global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            for (auto master_equation_id : master_equation_ids) {
                global_constraint->second->AddMaster(master_equation_id, relation_matrix(slave_count, master_count));
                master_count++;
            }
            slave_count++;
        }
        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::AssembleSlaves failed ...");
    }



    /**
     * @brief This method resets the LHS and RHS values of the AuxilaryGlobalMasterSlaveRelation objects
     */
    void ResetConstraintRelations()
    {
        KRATOS_TRY

        const int number_of_constraints = static_cast<int>(mGlobalMasterSlaveConstraints.size());

        // Getting the beginning iterator
        const auto constraints_begin = mGlobalMasterSlaveConstraints.begin();

        #pragma omp parallel for schedule(guided, 512)
        for (int i_constraints = 0; i_constraints < number_of_constraints; ++i_constraints) {
            auto it = constraints_begin;
            std::advance(it, i_constraints);

            (it->second)->Reset();
        }

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::ResetConstraintRelations failed to reset constraint relations..");
    }

    /**
     * @brief This method uses the MasterSlaveConstraints objects in rModelPart to reconstruct the LHS and RHS values of the AuxilaryGlobalMasterSlaveRelation objects. That is the value of Slave as LHS and the T*M+C as RHS value
     * @param rModelPart The model part of the problem to solve
     */
    void UpdateConstraintsForBuilding(ModelPart& rModelPart)
    {
        KRATOS_TRY
        // Reset the constraint equations
        ResetConstraintRelations();
        // Getting the array of the conditions
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        // Getting the beginning iterator
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        #pragma omp parallel for schedule(guided, 512)
        for (int i_constraints = 0; i_constraints < number_of_constraints; ++i_constraints) {
            auto it = constraints_begin;
            std::advance(it, i_constraints);

            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool constraint_is_active = true;
            if ((it)->IsDefined(ACTIVE))
                constraint_is_active = (it)->Is(ACTIVE);

            if (constraint_is_active) {
                UpdateMasterSlaveConstraint(*it, r_current_process_info);
            }
        }
        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::UpdateConstraintsForBuilding failed ..");
    }


    /**
     * @brief This method uses the MasterSlaveConstraints objects in rModelPart to reconstruct the LHS and RHS values of the individual AuxilaryGlobalMasterSlaveRelation object. That is the value of Slave as LHS and the T*M+C as RHS value
     * @param rMasterSlaveConstraint The MasterSlaveConstraint which is to be updated
     */
    void UpdateMasterSlaveConstraint(
        ModelPart::MasterSlaveConstraintType& rMasterSlaveConstraint,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        // Contributions to the system
        LocalSystemMatrixType relation_matrix(0,0);
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
        for (auto& slave_dof : slave_dofs_vector) {
            double slave_value_calc = 0.0;
            for (IndexType master_index = 0; master_index < master_dofs_vector.size(); ++master_index) {
                slave_value_calc += master_dofs_vector[master_index]->GetSolutionStepValue() * relation_matrix(slave_index, master_index);
            }
            slave_value_calc += constant_vector[slave_index];
            auto global_constraint = mGlobalMasterSlaveConstraints.find(slave_dof->EquationId());
            global_constraint->second->SetLeftHandSide( slave_dof->GetSolutionStepValue() );
            global_constraint->second->UpdateRightHandSide(slave_value_calc);
            slave_index++;
        }

        KRATOS_CATCH("ResidualBasedEliminationBuilderAndSolverWithConstraints::UpdateMasterSlaveConstraint failed ..");
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
        const int number_of_constraints = static_cast<int>(mGlobalMasterSlaveConstraints.size());
        // Getting the beginning iterator

        const auto constraints_begin = mGlobalMasterSlaveConstraints.begin();
        //contributions to the system
        VectorType master_weights_vector;
        double constant = 0.0;

        IndexType slave_equation_id = 0;
        EquationIdVectorType master_equation_ids = EquationIdVectorType(0);

#pragma omp parallel for schedule(guided, 512) firstprivate(slave_equation_id, master_equation_ids, master_weights_vector, constant)
        for (int i_constraints = 0; i_constraints < number_of_constraints; ++i_constraints) {
            //auto it = constraints_begin + i_constraints;
            auto it = constraints_begin;
            std::advance(it, i_constraints);

            double slave_dx_value = 0.0;
            //get the equation Ids of the constraint
            (it->second)->EquationIdsVector(slave_equation_id, master_equation_ids);
            //calculate constraint's T and b matrices
            (it->second)->CalculateLocalSystem(master_weights_vector, constant);
            int master_index = 0;
            for (auto& master_equation_id : master_equation_ids) {
                slave_dx_value += TSparseSpace::GetValue(rDx, master_equation_id) * master_weights_vector(master_index);
                ++master_index;
            }
            slave_dx_value += constant;

            rDx[slave_equation_id] = slave_dx_value; // this access is always unique for an object so no need of special care for openmp
        }

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
