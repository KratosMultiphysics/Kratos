//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:       BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER

/* System includes */
#include <set>
#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#else
#include <unordered_set>
#endif

/* Project includes */
#include "utilities/timer.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"

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
 * @class ResidualBasedEliminationBuilderAndSolver
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedEliminationBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedEliminationBuilderAndSolverWithConstraints
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedEliminationBuilderAndSolver);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
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

    /// Definition of the equation id vector
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    /// Definition of the compressed matrix
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// Node definition
    typedef Node<3> NodeType;

    /// Containers definition
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ResidualBasedEliminationBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedEliminationBuilderAndSolver() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        auto el_begin = rModelPart.ElementsBegin();
        auto cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements, nconditions,  LHS_Contribution, RHS_Contribution, equation_id )
        {
            #pragma omp  for schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++) {
                auto it_elem = el_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool element_is_active = true;
                if ((it_elem)->IsDefined(ACTIVE))
                    element_is_active = (it_elem)->Is(ACTIVE);

                if (element_is_active) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*(it_elem.base()), LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id, mLockArray);
#else
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id);
#endif
                    // Clean local elemental memory
                    pScheme->CleanMemory(*(it_elem.base()));
                }
            }

            #pragma omp  for schedule(guided, 512)
            for (int k = 0; k < nconditions; k++) {
                auto it_cond = cond_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool condition_is_active = true;
                if ((it_cond)->IsDefined(ACTIVE))
                    condition_is_active = (it_cond)->Is(ACTIVE);

                if (condition_is_active) {
                    // Calculate elemental contribution
                    pScheme->Condition_CalculateSystemContributions(*(it_cond.base()), LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);

#ifdef USE_LOCKS_IN_ASSEMBLY
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id, mLockArray);
#else
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id);
#endif

                    // Clean local elemental memory
                    pScheme->CleanMemory(*(it_cond.base()));
                }
            }
        }

        double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0) << "Build time: " << stop_build - start_build << std::endl;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation choosen the size of the matrix could
     * be equal to the total number of Dofs or to the number of unrestrained dofs
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     */
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA
        ) override
    {
        KRATOS_TRY

        // Getting the elements from the model
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();

        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Assemble all elements
        for (auto it_elem = r_elements_array.ptr_begin(); it_elem != r_elements_array.ptr_end(); ++it_elem) {
            // Calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it_elem, LHS_Contribution, equation_id, r_current_process_info);

            // Assemble the elemental contribution
            AssembleLHS(rA, LHS_Contribution, equation_id);

            // Clean local elemental memory
            pScheme->CleanMemory(*it_elem);
        }

        LHS_Contribution.resize(0, 0, false);

        // Assemble all conditions
        for (auto it_cond = r_conditions_array.ptr_begin(); it_cond != r_conditions_array.ptr_end(); ++it_cond) {
            // Calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it_cond, LHS_Contribution, equation_id, r_current_process_info);

            // Assemble the elemental contribution
            AssembleLHS(rA, LHS_Contribution, equation_id);
        }

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
        TSystemMatrixType& rA
        ) override
    {
        KRATOS_TRY

        // Getting the elements from the model
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        // Assemble all elements
        for (auto it_elem = r_elements_array.ptr_begin(); it_elem != r_elements_array.ptr_end(); ++it_elem) {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it_elem, LHS_Contribution, equation_id, r_current_process_info);

            //assemble the elemental contribution
            AssembleLHSCompleteOnFreeRows(rA, LHS_Contribution, equation_id);

            // clean local elemental memory
            pScheme->CleanMemory(*it_elem);
        }

        LHS_Contribution.resize(0, 0, false);
        // Assemble all conditions
        for (auto it_cond = r_conditions_array.ptr_begin(); it_cond != r_conditions_array.ptr_end(); ++it_cond) {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it_cond, LHS_Contribution, equation_id, r_current_process_info);

            //assemble the elemental contribution
            AssembleLHSCompleteOnFreeRows(rA, LHS_Contribution, equation_id);
        }

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
        if (TSparseSpace::Size(rb) != 0)
            norm_b = TSparseSpace::TwoNorm(rb);
        else
            norm_b = 0.0;

        if (norm_b != 0.0) {
            // Do solve
            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        } else
            TSparseSpace::SetToZero(rDx);

        // Prints informations about the current time
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")

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
            // Provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(rA, rDx, rb, BaseType::mDofSet, rModelPart);

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
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, rModelPart, rA, rb);

        Timer::Stop("Build");

//         ApplyPointLoads(pScheme,rModelPart, rb);

        // Does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, b);
        SystemSolve(A, Dx, b);

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
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Resetting to zero the vector of reactions
        if(BaseType::mCalculateReactionsFlag) {
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));
        }

        // Getting the Elements
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        // Assemble all elements
        #pragma omp parallel firstprivate( RHS_Contribution, equation_id)
        {
            const int nelements = static_cast<int>(r_elements_array.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i<nelements; ++i) {
                auto it_elem = r_elements_array.begin() + i;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool element_is_active = true;
                if (it_elem->IsDefined(ACTIVE))
                    element_is_active = it_elem->Is(ACTIVE);

                if (element_is_active) {
                    // Calculate elemental Right Hand Side Contribution
                    pScheme->Calculate_RHS_Contribution(*(it_elem.base()), RHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleRHS(rb, RHS_Contribution, equation_id);
                }
            }

            // Assemble all conditions
            const int nconditions = static_cast<int>(r_conditions_array.size());
            #pragma omp  for schedule(guided, 512)
            for (int i = 0; i<nconditions; ++i) {
                auto it_cond = r_conditions_array.begin() + i;
                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool condition_is_active = true;
                if (it_cond->IsDefined(ACTIVE))
                    condition_is_active = it_cond->Is(ACTIVE);

                if (condition_is_active) {
                    // Calculate elemental contribution
                    pScheme->Condition_Calculate_RHS_Contribution(*(it_cond.base()), RHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleRHS(rb, RHS_Contribution, equation_id);
                }
            }
        }

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
        KRATOS_TRY;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0) << "Setting up the dofs" << std::endl;

        // Gets the array of elements from the modeler
        ElementsArrayType& r_elements_array = rModelPart.Elements();
        const int nelements = static_cast<int>(r_elements_array.size());

        DofsVectorType elemental_dof_list;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        SizeType nthreads = OpenMPUtils::GetNumThreads();

#ifdef USE_GOOGLE_HASH
        typedef google::dense_hash_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#else
        typedef std::unordered_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#endif
        std::vector<set_type> dofs_aux_list(nthreads);

        for (int i = 0; i < static_cast<int>(nthreads); ++i) {
#ifdef USE_GOOGLE_HASH
            dofs_aux_list[i].set_empty_key(NodeType::DofType::Pointer());
#else
            dofs_aux_list[i].reserve(nelements);
#endif
        }

        #pragma omp parallel for firstprivate(nelements, elemental_dof_list)
        for (int i = 0; i < static_cast<int>(nelements); ++i) {
            auto it_elem = r_elements_array.begin() + i;
            const IndexType this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*(it_elem.base()), elemental_dof_list, r_current_process_info);

            dofs_aux_list[this_thread_id].insert(elemental_dof_list.begin(), elemental_dof_list.end());
        }

        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
        const int nconditions = static_cast<int>(r_conditions_array.size());
        #pragma omp parallel for firstprivate(nconditions, elemental_dof_list)
        for (int i = 0; i < nconditions; ++i) {
            auto it_cond = r_conditions_array.begin() + i;
            const IndexType this_thread_id = OpenMPUtils::ThisThread();

            // Gets list of Dof involved on every element
            pScheme->GetConditionDofList(*(it_cond.base()), elemental_dof_list, r_current_process_info);
            dofs_aux_list[this_thread_id].insert(elemental_dof_list.begin(), elemental_dof_list.end());
        }

        // Here we do a reduction in a tree so to have everything on thread 0
        SizeType old_max = nthreads;
        SizeType new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max >= 1 && new_max != old_max) {
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); ++i) {
                if (i + new_max < old_max) {
                    dofs_aux_list[i].insert(dofs_aux_list[i + new_max].begin(), dofs_aux_list[i + new_max].end());
                    dofs_aux_list[i + new_max].clear();
                }
            }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        DofsArrayType dof_temp;
        BaseType::mDofSet = DofsArrayType();

        dof_temp.reserve(dofs_aux_list[0].size());
        for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); ++it) {
            dof_temp.push_back(it->get());
        }
        dof_temp.Sort();

        BaseType::mDofSet = dof_temp;

        // Throws an execption if there are no Degrees of freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the dofs" << std::endl;

#ifdef _OPENMP
        if (mLockArray.size() != 0) {
            for (int i = 0; i < static_cast<int>(mLockArray.size()); i++)
                omp_destroy_lock(&mLockArray[i]);
        }

        mLockArray.resize(BaseType::mDofSet.size());

        for (int i = 0; i < static_cast<int>(mLockArray.size()); i++)
            omp_init_lock(&mLockArray[i]);
#endif

        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
#ifdef KRATOS_DEBUG
        if(BaseType::GetCalculateReactionsFlag()) {
            for(auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                    KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " << std::endl
                        << "Node : " << dof_iterator->Id() << std::endl
                        << "Dof : " << (*dof_iterator) << std::endl << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif

        KRATOS_CATCH("");
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(ModelPart& rModelPart) override
    {
        // Set equation id for degrees of freedom
        // the free degrees of freedom are positioned at the beginning of the system,
        // while the fixed one are at the end (in opposite order).
        //
        // that means that if the EquationId is greater than "mEquationSystemSize"
        // the pointed degree of freedom is restrained
        //
        int free_id = 0;
        int fix_id = BaseType::mDofSet.size();

        for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->IsFixed())
                dof_iterator->SetEquationId(--fix_id);
            else
                dof_iterator->SetEquationId(free_id++);

        BaseType::mEquationSystemSize = fix_id;
    }

    /**
     * @brief This method resize and initializes the system of euqations
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
        KRATOS_TRY

        if (pA == NULL) { // If the pointer is not initialized initialize it to an empty matrix
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) { // If the pointer is not initialized initialize it to an empty matrix
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) { // If the pointer is not initialized initialize it to an empty matrix
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }
        if (BaseType::mpReactionsVector == NULL) { // If the pointer is not initialized initialize it to an empty matrix
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        TSystemMatrixType& rA = *pA;
        TSystemVectorType& rDx = *pDx;
        TSystemVectorType& rb = *pb;

        // Resizing the system vectors and matrix
        if (rA.size1() == 0 || BaseType::GetReshapeMatrixFlag()) { // If the matrix is not initialized
            rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(pScheme, rA, rModelPart);
        } else {
            if (rA.size1() != BaseType::mEquationSystemSize || rA.size2() != BaseType::mEquationSystemSize) {
                KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, rA, rModelPart);
            }
        }
        if (rDx.size() != BaseType::mEquationSystemSize)
            rDx.resize(BaseType::mEquationSystemSize, false);
        if (rb.size() != BaseType::mEquationSystemSize)
            rb.resize(BaseType::mEquationSystemSize, false);

        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag) {
            const SizeType ReactionsVectorSize = BaseType::mDofSet.size();
            if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
        }

        KRATOS_CATCH("")
    }


    /**
     * @brief This method computes the reactions
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part considered
     * @param rA The LHS of the system
     * @param rDx The vector of Unknowns
     * @param rb The RHS vector
     */
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        //refresh RHS to have the correct reactions
        BuildRHS(pScheme, rModelPart, rb);

        int i;
        int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size(*BaseType::mpReactionsVector);

        //updating variables
        TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
        for (auto it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2) {
            i = (*it2)->EquationId();
            i -= systemsize;
            (*it2)->GetSolutionStepReactionValue() = -ReactionsVector[i];

        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver choosen
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
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        this->mDofSet = DofsArrayType();

        if (this->mpReactionsVector != NULL)
            TSparseSpace::Clear((this->mpReactionsVector));
//             this->mReactionsVector = TSystemVectorType();

        this->mpLinearSystemSolver->Clear();

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 1) << "Clear Function called" << std::endl;
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
        return "ResidualBasedEliminationBuilderAndSolver";
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

#ifdef _OPENMP
   std::vector<omp_lock_t> mLockArray;
#endif

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

   /**
    * @brief This method assembles the system
    * @param rA The LHS of the system
    * @param rb The RHS of the system
    * @param rLHSContribution The LHS local contribution
    * @param rRHSContribution The RHS local contribution
    * @param rEquationId The equation id
    * @param rLockArray The lock of the dof
    */
    void Assemble(
        TSystemMatrixType& rA,
        TSystemVectorType& rb,
        const LocalSystemMatrixType& rLHSContribution,
        const LocalSystemVectorType& rRHSContribution,
        const EquationIdVectorType& rEquationId
#ifdef USE_LOCKS_IN_ASSEMBLY
        ,std::vector< omp_lock_t >& rLockArray
#endif
        )
    {
        SizeType local_size = rLHSContribution.size1();

        for (IndexType i_local = 0; i_local < local_size; i_local++) {
            IndexType i_global = rEquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize) {
#ifdef USE_LOCKS_IN_ASSEMBLY
                omp_set_lock(&rLockArray[i_global]);
#endif
                rb[i_global] += rRHSContribution(i_local);
                for (IndexType j_local = 0; j_local < local_size; j_local++) {
                    IndexType j_global = rEquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize) {
                        rA(i_global, j_global) += rLHSContribution(i_local, j_local);
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
    virtual void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& rA,
        ModelPart& rModelPart
        )
    {
        // Filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

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

        const int nelements = static_cast<int>(rModelPart.Elements().size());
        #pragma omp parallel for firstprivate(nelements, ids)
        for (int i_elem = 0; i_elem < nelements; ++i_elem) {
            auto it_elem = rModelPart.Elements().begin() + i_elem;
            pScheme->EquationId(*(it_elem.base()), ids, r_process_info);
            for (IndexType i = 0; i < ids.size(); i++) {
            #ifdef _OPENMP
                omp_set_lock(&mLockArray[ids[i]]);
            #endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

            #ifdef _OPENMP
                omp_unset_lock(&mLockArray[ids[i]]);
            #endif
            }
        }

        const int nconditions = static_cast<int>(rModelPart.Conditions().size());
        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int i_cond = 0; i_cond < nconditions; i_cond++) {
            auto it_cond = rModelPart.Conditions().begin() + i_cond;
            pScheme->Condition_EquationId(*(it_cond.base()), ids, r_process_info);
            for (IndexType i = 0; i < ids.size(); i++) {
            #ifdef _OPENMP
                omp_set_lock(&mLockArray[ids[i]]);
            #endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
            #ifdef _OPENMP
                omp_unset_lock(&mLockArray[ids[i]]);
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
     * @brief This method assembles the LHS of the system
     * @param rA The LHS to assemble
     * @param rLHSContribution The local LHS contribution
     * @param rEquationId The equation id
     */
    void AssembleLHS(
        TSystemMatrixType& rA,
        LocalSystemMatrixType& rLHSContribution,
        EquationIdVectorType& rEquationId
        )
    {
        SizeType local_size = rLHSContribution.size1();

        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            IndexType i_global = rEquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize) {
                for (IndexType j_local = 0; j_local < local_size; ++j_local) {
                    IndexType j_global = rEquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
                        rA(i_global, j_global) += rLHSContribution(i_local, j_local);
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method ensures that the contribution is unique
     */
    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while (i != endit && (*i) != candidate) {
            i++;
        }
        if (i == endit) {
            v.push_back(candidate);
        }

    }

    /**
     * @brief This method assembles the RHS of the system
     * @param rb The RHS to assemble
     * @param rRHSContribution The local RHS contribution
     * @param rEquationId The equation id
     */
    void AssembleRHS(
        TSystemVectorType& rb,
        const LocalSystemVectorType& rRHSContribution,
        const EquationIdVectorType& rEquationId
        )
    {
        SizeType local_size = rRHSContribution.size();

        if (BaseType::mCalculateReactionsFlag == false) {
            for (IndexType i_local = 0; i_local < local_size; ++i_local) {
                const IndexType i_global = rEquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) { // Free dof
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = rb[i_global];
                    const double& rhs_value = rRHSContribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
            }
        } else {
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
            for (IndexType i_local = 0; i_local < local_size; ++i_local) {
                const IndexType i_global = rEquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) { //free dof
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = rb[i_global];
                    const double& rhs_value = rRHSContribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                } else { // Fixed dof
                    double& b_value = ReactionsVector[i_global - BaseType::mEquationSystemSize];
                    const double& rhs_value = rRHSContribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
            }
        }
    }

    /**
     * @brief This method assembles the LHS of the system (on free rows)
     * @param rA The LHS to assemble
     * @param rLHSContribution The local LHS contribution
     * @param rEquationId The equation id
     */
    void AssembleLHSCompleteOnFreeRows(
        TSystemMatrixType& rA,
        LocalSystemMatrixType& rLHSContribution,
        EquationIdVectorType& rEquationId
        )
    {
        SizeType local_size = rLHSContribution.size1();
        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            IndexType i_global = rEquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize) {
                for (IndexType j_local = 0; j_local < local_size; ++j_local) {
                    int j_global = rEquationId[j_local];
                    rA(i_global, j_global) += rLHSContribution(i_local, j_local);
                }
            }
        }
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

    ///@}

}; /* Class ResidualBasedEliminationBuilderAndSolver */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

