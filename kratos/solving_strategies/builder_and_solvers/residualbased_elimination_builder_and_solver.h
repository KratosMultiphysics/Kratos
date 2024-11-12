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
//  Collaborators:   Vicente Mataix
//
//

#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER

/* System includes */
#include <set>
#include <unordered_set>

/* External includes */
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "utilities/timer.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"
#include "utilities/builtin_timer.h"
#include "utilities/atomic_utilities.h"
#include "spaces/ublas_space.h"

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
 * @brief Current class provides an implementation for standard  elimination builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains this information.
 * Calculation of the reactions involves a cost very similar to the calculation of the total residual
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

    /// The definition of the current class
    typedef ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

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

    /// Definition of the equation id vector
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    /// Node definition
    typedef Node NodeType;

    /// Containers definition
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedEliminationBuilderAndSolver() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedEliminationBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Constructor.
     */
    explicit ResidualBasedEliminationBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedEliminationBuilderAndSolver() override
    {
    }

    /**
     * @brief Create method
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(
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
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();

        // Getting the elements from the model
        const int nelements = static_cast<int>(r_elements_array.size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(r_conditions_array.size());

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const auto it_elem_begin = r_elements_array.begin();
        const auto it_cond_begin = r_conditions_array.begin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        // Assemble all elements
        const auto timer = BuiltinTimer();

        #pragma omp parallel firstprivate(LHS_Contribution, RHS_Contribution, equation_id )
        {
            #pragma omp  for schedule(guided, 512) nowait
            for (int k = 0; k < nelements; ++k) {
                auto it_elem = it_elem_begin + k;

                // If the element is active
                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_elem, LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);

                    KRATOS_WATCH("Inside Elimination Builder element loop")
                    KRATOS_WATCH("Before Assemble")
                    KRATOS_WATCH(LHS_Contribution)
                    KRATOS_WATCH(RHS_Contribution)
                    KRATOS_WATCH(equation_id)
                    // Assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id, mLockArray);
#else
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id);

                    KRATOS_WATCH("Inside Elimination Builder element loop")
                    KRATOS_WATCH("After Assemble")
                    KRATOS_WATCH(rA)
                    KRATOS_WATCH(rb)
                    KRATOS_WATCH(LHS_Contribution)
                    KRATOS_WATCH(RHS_Contribution)
                    KRATOS_WATCH(equation_id)

#endif
                }
            }

            #pragma omp  for schedule(guided, 512)
            for (int k = 0; k < nconditions; ++k) {
                auto it_cond = it_cond_begin + k;

                // If the condition is active
                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_cond, LHS_Contribution, RHS_Contribution, equation_id, r_current_process_info);
                    KRATOS_WATCH("Inside Elimination Builder condition loop")
                    KRATOS_WATCH("Before Assemble")
                    KRATOS_WATCH(LHS_Contribution)
                    KRATOS_WATCH(RHS_Contribution)
                    KRATOS_WATCH(equation_id)
#ifdef USE_LOCKS_IN_ASSEMBLY
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id, mLockArray);
#else
                    Assemble(rA, rb, LHS_Contribution, RHS_Contribution, equation_id);
                    KRATOS_WATCH("Inside Elimination Builder condition loop")
                    KRATOS_WATCH("After Assemble")
                    KRATOS_WATCH(rA)
                    KRATOS_WATCH(rb)
                    KRATOS_WATCH(LHS_Contribution)
                    KRATOS_WATCH(RHS_Contribution)
                    KRATOS_WATCH(equation_id)
#endif
                }
            }
        }
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() >=1) << "System build time: " << timer.ElapsedSeconds() << std::endl;


        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 2) << "Finished building" << std::endl;


        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation chosen the size of the matrix could be equal to the total number of Dofs or to the number of unrestrained dofs
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

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();

        // Getting the elements from the model
        const int nelements = static_cast<int>(r_elements_array.size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(r_conditions_array.size());

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const auto it_elem_begin = r_elements_array.begin();
        const auto it_cond_begin = r_conditions_array.begin();

        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        #pragma omp parallel firstprivate(LHS_Contribution, equation_id )
        {
            #pragma omp  for schedule(guided, 512) nowait
            for (int k = 0; k < nelements; ++k) {
                auto it_elem = it_elem_begin + k;

                // If the element is active
                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_elem, LHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHS(rA, LHS_Contribution, equation_id);
                }
            }

            #pragma omp  for schedule(guided, 512)
            for (int k = 0; k < nconditions; ++k) {
                auto it_cond = it_cond_begin + k;

                // If the condition is active
                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_cond, LHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHS(rA, LHS_Contribution, equation_id);
                }
            }
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

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        ElementsArrayType& r_elements_array = rModelPart.Elements();

        // Getting the array of the conditions
        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();

        // Getting the elements from the model
        const int nelements = static_cast<int>(r_elements_array.size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(r_conditions_array.size());

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const auto it_elem_begin = r_elements_array.begin();
        const auto it_cond_begin = r_conditions_array.begin();

        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        #pragma omp parallel firstprivate(LHS_Contribution, equation_id )
        {
            #pragma omp  for schedule(guided, 512) nowait
            for (int k = 0; k < nelements; ++k) {
                auto it_elem = it_elem_begin + k;

                // If the element is active
                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_elem, LHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHSCompleteOnFreeRows(rA, LHS_Contribution, equation_id);
                }
            }

            #pragma omp  for schedule(guided, 512)
            for (int k = 0; k < nconditions; ++k) {
                auto it_cond = it_cond_begin + k;

                // If the condition is active
                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_cond, LHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHSCompleteOnFreeRows(rA, LHS_Contribution, equation_id);
                }
            }
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
        if (TSparseSpace::Size(rb) != 0) {
            norm_b = TSparseSpace::TwoNorm(rb);
        } else {
            norm_b = 0.0;
        }

        if (norm_b != 0.0) {
            // Do solve
            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        } else
            TSparseSpace::SetToZero(rDx);

        // Prints information about the current time
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
        if (TSparseSpace::Size(rb) != 0) {
            norm_b = TSparseSpace::TwoNorm(rb);
        } else {
            norm_b = 0.0;
        }

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

        // Prints information about the current time
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

        KRATOS_WATCH("Inside Elimination Builder")
        KRATOS_WATCH("After Build")
        KRATOS_WATCH(rA)
        KRATOS_WATCH(rb)

        Timer::Stop("Build");

        // Does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() >=1) << "System solve time: " << timer.ElapsedSeconds() << std::endl;


        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, rb);
        SystemSolve(rA, rDx, rb);

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

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;

        // Assemble all elements
        #pragma omp parallel firstprivate( RHS_Contribution, equation_id)
        {
            const auto it_elem_begin = r_elements_array.begin();
            const int nelements = static_cast<int>(r_elements_array.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < nelements; ++i) {
                auto it_elem = it_elem_begin + i;

                // If the element is active
                if (it_elem->IsActive()) {
                    // Calculate elemental Right Hand Side Contribution
                    pScheme->CalculateRHSContribution(*it_elem, RHS_Contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleRHS(rb, RHS_Contribution, equation_id);
                }
            }

            // Assemble all conditions
            const auto it_cond_begin = r_conditions_array.begin();
            const int nconditions = static_cast<int>(r_conditions_array.size());
            #pragma omp  for schedule(guided, 512)
            for (int i = 0; i < nconditions; ++i) {
                auto it_cond = it_cond_begin + i;
                // If the condition is active
                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateRHSContribution(*it_cond, RHS_Contribution, equation_id, r_current_process_info);

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
     * @details The list of dofs is stores inside the BuilderAndSolver as it is closely connected to the
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

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        SizeType nthreads = ParallelUtilities::GetNumThreads();

        typedef std::unordered_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;

        std::vector<set_type> dofs_aux_list(nthreads);

        for (int i = 0; i < static_cast<int>(nthreads); ++i) {
            dofs_aux_list[i].reserve(nelements);
        }

        IndexPartition<std::size_t>(nelements).for_each(elemental_dof_list, [&](std::size_t Index, DofsVectorType& tls_elemental_dof_list){
            auto it_elem = r_elements_array.begin() + Index;
            const IndexType this_thread_id = OpenMPUtils::ThisThread();

            // Gets list of Dof involved on every element
            pScheme->GetDofList(*it_elem, tls_elemental_dof_list, r_current_process_info);

            dofs_aux_list[this_thread_id].insert(tls_elemental_dof_list.begin(), tls_elemental_dof_list.end());
        });

        ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
        const int nconditions = static_cast<int>(r_conditions_array.size());

        IndexPartition<std::size_t>(nconditions).for_each(elemental_dof_list, [&](std::size_t Index, DofsVectorType& tls_elemental_dof_list){
            auto it_cond = r_conditions_array.begin() + Index;
            const IndexType this_thread_id = OpenMPUtils::ThisThread();

            // Gets list of Dof involved on every element
            pScheme->GetDofList(*it_cond, tls_elemental_dof_list, r_current_process_info);
            dofs_aux_list[this_thread_id].insert(tls_elemental_dof_list.begin(), tls_elemental_dof_list.end());
        });

        // Here we do a reduction in a tree so to have everything on thread 0
        SizeType old_max = nthreads;
        SizeType new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max >= 1 && new_max != old_max) {
            IndexPartition<std::size_t>(new_max).for_each([&](std::size_t Index){
                if (Index + new_max < old_max) {
                    dofs_aux_list[Index].insert(dofs_aux_list[Index + new_max].begin(), dofs_aux_list[Index + new_max].end());
                    dofs_aux_list[Index + new_max].clear();
                }
            });

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        DofsArrayType dof_temp;
        BaseType::mDofSet = DofsArrayType();

        dof_temp.reserve(dofs_aux_list[0].size());
        for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); ++it) {
            dof_temp.push_back(*it);
        }
        dof_temp.Sort();

        BaseType::mDofSet = dof_temp;

        // Throws an exception if there are no Degrees of freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the dofs" << std::endl;

#ifdef USE_LOCKS_IN_ASSEMBLY
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

        if (pA == nullptr) { // If the pointer is not initialized initialize it to an empty matrix

            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == nullptr) { // If the pointer is not initialized initialize it to an empty matrix

            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == nullptr) { // If the pointer is not initialized initialize it to an empty matrix

            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }
        if (BaseType::mpReactionsVector == nullptr) { // If the pointer is not initialized initialize it to an empty matrix

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
                KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permitted."<<std::endl;
                rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, rA, rModelPart);
            }
        }
        if (rDx.size() != BaseType::mEquationSystemSize) {
            rDx.resize(BaseType::mEquationSystemSize, false);
        }
        TSparseSpace::SetToZero(rDx);
        if (rb.size() != BaseType::mEquationSystemSize) {
            rb.resize(BaseType::mEquationSystemSize, false);
        }
        TSparseSpace::SetToZero(rb);

        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag == true) {
            const std::size_t reactions_vector_size = BaseType::mDofSet.size() - BaseType::mEquationSystemSize;
            if (BaseType::mpReactionsVector->size() != reactions_vector_size)
                BaseType::mpReactionsVector->resize(reactions_vector_size, false);
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

        // Updating variables
        std::size_t i;
        TSystemVectorType& r_reactions_vector = *BaseType::mpReactionsVector;
        for (auto it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2) {
            i = (*it2)->EquationId();
            if (i >= BaseType::mEquationSystemSize) {
                i -= BaseType::mEquationSystemSize;
                (*it2)->GetSolutionStepReactionValue() = -r_reactions_vector[i];
            }

        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation chosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver chosen
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
        // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal); 
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        this->mDofSet = DofsArrayType();
        this->mpReactionsVector.reset();
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

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "elimination_builder_and_solver",
            "block_builder"                        : false,
            "diagonal_values_for_dirichlet_dofs"   : "use_max_diagonal"
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
        return "elimination_builder_and_solver";
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

#ifdef USE_LOCKS_IN_ASSEMBLY
   std::vector<omp_lock_t> mLockArray; /// TODO: Replace with std::vector<LockObject>
#endif

    double mScaleFactor = 1.0;         /// The manually set scale factor

    SCALING_DIAGONAL mScalingDiagonal = SCALING_DIAGONAL::NO_SCALING;; /// We identify the scaling considered for the dirichlet dofs


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
    * @note The main difference respect the block builder and solver is the fact that the fixed DoFs are not considered on the assembling
    */
    void Assemble(
        TSystemMatrixType& rA,
        TSystemVectorType& rb,
        const LocalSystemMatrixType& rLHSContribution,
        const LocalSystemVectorType& rRHSContribution,
        const Element::EquationIdVectorType& rEquationId
#ifdef USE_LOCKS_IN_ASSEMBLY
        ,std::vector< omp_lock_t >& rLockArray
#endif
        )
    {
        const SizeType local_size = rLHSContribution.size1();

        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            const IndexType i_global = rEquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize) {
#ifdef USE_LOCKS_IN_ASSEMBLY
                omp_set_lock(&rLockArray[i_global]);
                rb[i_global] += rRHSContribution(i_local);
#else
                double& r_a = rb[i_global];
                const double& v_a = rRHSContribution(i_local);
                AtomicAdd(r_a, v_a);
#endif
                AssembleRowContributionFreeDofs(rA, rLHSContribution, i_global, i_local, rEquationId);

#ifdef USE_LOCKS_IN_ASSEMBLY
                omp_unset_lock(&rLockArray[i_global]);
#endif
            }
            //note that computation of reactions is not performed here!
        }
    }

    /**
     * @brief This method constructs the relationship between the DoF
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

        std::vector<std::unordered_set<IndexType> > indices(equation_size);

        block_for_each(indices, [](std::unordered_set<IndexType>& rIndices){
            rIndices.reserve(40);
        });

        Element::EquationIdVectorType ids(3, 0);

        #pragma omp parallel firstprivate(ids)
        {
            // The process info
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // We repeat the same declaration for each thead
            std::vector<std::unordered_set<IndexType> > temp_indexes(equation_size);

            #pragma omp for
            for (int index = 0; index < static_cast<int>(equation_size); ++index)
                temp_indexes[index].reserve(30);

            // Getting the size of the array of elements from the model
            const int number_of_elements = static_cast<int>(rModelPart.Elements().size());

            // Element initial iterator
            const auto it_elem_begin = rModelPart.ElementsBegin();

            // We iterate over the elements
            #pragma omp for schedule(guided, 512) nowait
            for (int i_elem = 0; i_elem<number_of_elements; ++i_elem) {
                auto it_elem = it_elem_begin + i_elem;
                pScheme->EquationId( *it_elem, ids, r_current_process_info);

                for (auto& id_i : ids) {
                    if (id_i < BaseType::mEquationSystemSize) {
                        auto& row_indices = temp_indexes[id_i];
                        for (auto& id_j : ids)
                            if (id_j < BaseType::mEquationSystemSize)
                                row_indices.insert(id_j);
                    }
                }
            }

            // Getting the size of the array of the conditions
            const int number_of_conditions = static_cast<int>(rModelPart.Conditions().size());

            // Condition initial iterator
            const auto it_cond_begin = rModelPart.ConditionsBegin();

            // We iterate over the conditions
            #pragma omp for schedule(guided, 512) nowait
            for (int i_cond = 0; i_cond<number_of_conditions; ++i_cond) {
                auto it_cond = it_cond_begin + i_cond;
                pScheme->EquationId( *it_cond, ids, r_current_process_info);

                for (auto& id_i : ids) {
                    if (id_i < BaseType::mEquationSystemSize) {
                        auto& row_indices = temp_indexes[id_i];
                        for (auto& id_j : ids)
                            if (id_j < BaseType::mEquationSystemSize)
                                row_indices.insert(id_j);
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

        rA = TSystemMatrixType(indices.size(), indices.size(), nnz);

        double* Avalues = rA.value_data().begin();
        std::size_t* Arow_indices = rA.index1_data().begin();
        std::size_t* Acol_indices = rA.index2_data().begin();

        // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (IndexType i = 0; i < rA.size1(); ++i)
            Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();

        IndexPartition<std::size_t>(rA.size1()).for_each([&](std::size_t Index){
            const IndexType row_begin = Arow_indices[Index];
            const IndexType row_end = Arow_indices[Index + 1];
            IndexType k = row_begin;
            for (auto it = indices[Index].begin(); it != indices[Index].end(); ++it) {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                ++k;
            }

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
        });

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
        const SizeType local_size = rLHSContribution.size1();

        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            const IndexType i_global = rEquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize) {
                for (IndexType j_local = 0; j_local < local_size; ++j_local) {
                    const IndexType j_global = rEquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize) {
                        rA(i_global, j_global) += rLHSContribution(i_local, j_local);
                    }
                }
            }
        }
    }

    /**
     * @brief This function is equivalent to the AssembleRowContribution of the block builder and solver
     * @note The main difference respect the block builder and solver is the fact that the fixed DoFs are skipped
     */
    inline void AssembleRowContributionFreeDofs(
        TSystemMatrixType& rA,
        const Matrix& rALocal,
        const IndexType i,
        const IndexType i_local,
        const Element::EquationIdVectorType& EquationId
        )
    {
        double* values_vector = rA.value_data().begin();
        IndexType* index1_vector = rA.index1_data().begin();
        IndexType* index2_vector = rA.index2_data().begin();

        const IndexType left_limit = index1_vector[i];

        // Find the first entry
        // We iterate over the equation ids until we find the first equation id to be considered
        // We count in which component we find an ID
        IndexType last_pos = 0;
        IndexType last_found = 0;
        IndexType counter = 0;
        for(IndexType j=0; j < EquationId.size(); ++j) {
            ++counter;
            const IndexType j_global = EquationId[j];
            if (j_global < BaseType::mEquationSystemSize) {
                last_pos = ForwardFind(j_global,left_limit,index2_vector);
                last_found = j_global;
                break;
            }
        }

        // If the counter is equal to the size of the EquationID vector that means that only one dof will be considered, if the number is greater means that all the dofs are fixed. If the number is below means that at we have several dofs free to be considered
        if (counter <= EquationId.size()) {
#ifndef USE_LOCKS_IN_ASSEMBLY
            double& r_a = values_vector[last_pos];
            const double& v_a = rALocal(i_local,counter - 1);
            AtomicAdd(r_a,  v_a);
#else
            values_vector[last_pos] += rALocal(i_local,counter - 1);
#endif
            // Now find all of the other entries
            IndexType pos = 0;
            for(IndexType j = counter; j < EquationId.size(); ++j) {
                IndexType id_to_find = EquationId[j];
                if (id_to_find < BaseType::mEquationSystemSize) {
                    if(id_to_find > last_found)
                        pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
                    else if(id_to_find < last_found)
                        pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
                    else
                        pos = last_pos;

#ifndef USE_LOCKS_IN_ASSEMBLY
                    double& r = values_vector[pos];
                    const double& v = rALocal(i_local,j);
                    AtomicAdd(r,  v);
#else
                    values_vector[pos] += rALocal(i_local,j);
#endif
                    last_found = id_to_find;
                    last_pos = pos;
                }
            }
        }
    }

    inline IndexType ForwardFind(const IndexType id_to_find,
                                   const IndexType start,
                                   const IndexType* index_vector)
    {
        IndexType pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline IndexType BackwardFind(const IndexType id_to_find,
                                    const IndexType start,
                                    const IndexType* index_vector)
    {
        IndexType pos = start;
        while(id_to_find != index_vector[pos]) pos--;
        return pos;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Setting flags<
        const std::string& r_diagonal_values_for_dirichlet_dofs = ThisParameters["diagonal_values_for_dirichlet_dofs"].GetString();

        std::set<std::string> available_options_for_diagonal = {"no_scaling","use_max_diagonal","use_diagonal_norm","defined_in_process_info"};

        if (available_options_for_diagonal.find(r_diagonal_values_for_dirichlet_dofs) == available_options_for_diagonal.end()) {
            std::stringstream msg;
            msg << "Currently prescribed diagonal values for dirichlet dofs : " << r_diagonal_values_for_dirichlet_dofs << "\n";
            msg << "Admissible values for the diagonal scaling are : no_scaling, use_max_diagonal, use_diagonal_norm, or defined_in_process_info" << "\n";
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

                    AtomicAdd(b_value, rhs_value);
                }
            }
        } else {
            TSystemVectorType& r_reactions_vector = *BaseType::mpReactionsVector;
            for (IndexType i_local = 0; i_local < local_size; ++i_local) {
                const IndexType i_global = rEquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) { //free dof
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = rb[i_global];
                    const double& rhs_value = rRHSContribution[i_local];

                    AtomicAdd(b_value, rhs_value);
                } else { // Fixed dof
                    double& b_value = r_reactions_vector[i_global - BaseType::mEquationSystemSize];
                    const double& rhs_value = rRHSContribution[i_local];

                    AtomicAdd(b_value, rhs_value);
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
        const SizeType local_size = rLHSContribution.size1();
        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            const IndexType i_global = rEquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize) {
                for (IndexType j_local = 0; j_local < local_size; ++j_local) {
                    const IndexType j_global = rEquationId[j_local];
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
