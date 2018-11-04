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
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedEliminationBuilderAndSolver);

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef Node<3> NodeType;

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
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
//         KRATOS_INFO("ResidualBasedEliminationBuilderAndSolver") << "Using the standard builder and solver " << std::endl;
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
     * @param A The LHS matrix
     * @param b The RHS vector
     */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        //getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        //getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements, nconditions,  LHS_Contribution, RHS_Contribution, EquationId )
        {
            #pragma omp  for schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mLockArray);
#else
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));

                }

            }

            #pragma omp  for schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

#ifdef USE_LOCKS_IN_ASSEMBLY
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mLockArray);
#else
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif

                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
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
     * @param A The LHS matrix
     */
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& rElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);
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
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& rElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);
        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);
        }


        KRATOS_CATCH("")

    }

    /**
     * @brief This is a call to the linear system solver
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
            TSparseSpace::SetToZero(Dx);

        // Prints informations about the current time
        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")

    }

    /**
      *@brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void SystemSolveWithPhysics(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        ModelPart& rModelPart
    )
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, rModelPart);

            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
        {
            TSparseSpace::SetToZero(Dx);
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
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

//         ApplyPointLoads(pScheme,rModelPart,b);

        // Does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

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
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        //resetting to zero the vector of reactions

        if(BaseType::mCalculateReactionsFlag)
        {
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));
        }

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& pConditions = rModelPart.Conditions();

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements

        #pragma omp parallel firstprivate( RHS_Contribution, EquationId)
        {
            const int nelements = static_cast<int>(pElements.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i<nelements; i++)
            {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    // Calculate elemental Right Hand Side Contribution
                    pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }

            // assemble all conditions
            const int nconditions = static_cast<int>(pConditions.size());
            #pragma omp  for schedule(guided, 512)
            for (int i = 0; i<nconditions; i++)
            {
                auto it = pConditions.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
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

        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = rModelPart.Elements();
        const int nelements = static_cast<int>(pElements.size());

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        unsigned int nthreads = OpenMPUtils::GetNumThreads();

//         typedef boost::fast_pool_allocator< NodeType::DofType::Pointer > allocator_type;
//         typedef std::unordered_set < NodeType::DofType::Pointer,
//             DofPointerHasher,
//             DofPointerComparor,
//             allocator_type    >  set_type;

#ifdef USE_GOOGLE_HASH
        typedef google::dense_hash_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#else
        typedef std::unordered_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#endif
      //


        std::vector<set_type> dofs_aux_list(nthreads);
//         std::vector<allocator_type> allocators(nthreads);

        for (int i = 0; i < static_cast<int>(nthreads); i++)
        {
#ifdef USE_GOOGLE_HASH
            dofs_aux_list[i].set_empty_key(NodeType::DofType::Pointer());
#else
//             dofs_aux_list[i] = set_type( allocators[i]);
            dofs_aux_list[i].reserve(nelements);
#endif
        }

#pragma omp parallel for firstprivate(nelements, ElementalDofList)
        for (int i = 0; i < static_cast<int>(nelements); i++)
        {
            typename ElementsArrayType::iterator it = pElements.begin() + i;
            const unsigned int this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);

            dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
        }

        ConditionsArrayType& pConditions = rModelPart.Conditions();
        const int nconditions = static_cast<int>(pConditions.size());
#pragma omp parallel for firstprivate(nconditions, ElementalDofList)
        for (int i = 0; i < nconditions; i++)
        {
            typename ConditionsArrayType::iterator it = pConditions.begin() + i;
            const unsigned int this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetConditionDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);
            dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
        }

        //here we do a reduction in a tree so to have everything on thread 0
        unsigned int old_max = nthreads;
        unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max >= 1 && new_max != old_max)
        {
//          //just for debugging
//          std::cout << "old_max" << old_max << " new_max:" << new_max << std::endl;
//          for (int i = 0; i < new_max; i++)
//          {
//             if (i + new_max < old_max)
//             {
//                std::cout << i << " - " << i + new_max << std::endl;
//             }
//          }
//          std::cout << "********************" << std::endl;

#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); i++)
            {
                if (i + new_max < old_max)
                {
                    dofs_aux_list[i].insert(dofs_aux_list[i + new_max].begin(), dofs_aux_list[i + new_max].end());
                    dofs_aux_list[i + new_max].clear();
                }
            }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dofs_aux_list[0].size());
        for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); it++)
        {
            Doftemp.push_back(it->get());
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        // Throws an execption if there are no Degrees of freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the dofs" << std::endl;

#ifdef _OPENMP
        if (mLockArray.size() != 0)
        {
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
        if(BaseType::GetCalculateReactionsFlag())
        {
            for(auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            {
                    KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " <<std::endl
                        << "Node : "<<dof_iterator->Id()<< std::endl
                        << "Dof : "<<(*dof_iterator)<<std::endl<<"Not possible to calculate reactions."<<std::endl;
            }
        }
#endif

        KRATOS_CATCH("");
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(
        ModelPart& rModelPart
        ) override
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

        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->IsFixed())
                dof_iterator->SetEquationId(--fix_id);
            else
                dof_iterator->SetEquationId(free_id++);

        BaseType::mEquationSystemSize = fix_id;

    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
    ) override
    {
        KRATOS_TRY
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }
        if (BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(pScheme, A, rModelPart);
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                //KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, A, rModelPart);
            }
        }
        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        if (b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize, false);

        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag == true)
        {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size();
            if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
        }

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        //refresh RHS to have the correct reactions
        BuildRHS(pScheme, rModelPart, b);

        int i;
        int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size(*BaseType::mpReactionsVector);

        typename DofsArrayType::ptr_iterator it2;

        //updating variables
        TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
        for (it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2)
        {
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
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
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

    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        const Element::EquationIdVectorType& EquationId
#ifdef USE_LOCKS_IN_ASSEMBLY
        ,std::vector< omp_lock_t >& lock_array
#endif
        )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize)
            {
#ifdef USE_LOCKS_IN_ASSEMBLY
                omp_set_lock(&lock_array[i_global]);
                b[i_global] += RHS_Contribution(i_local);
#else
                double& r_a = b[i_global];
                const double& v_a = RHS_Contribution(i_local);
                #pragma omp atomic
                r_a += v_a;
#endif
                AssembleRowContribution(A, LHS_Contribution, i_global, i_local, EquationId);

#ifdef USE_LOCKS_IN_ASSEMBLY
                omp_unset_lock(&lock_array[i_global]);
#endif
            }
            //note that computation of reactions is not performed here!
        }
    }


    //**************************************************************************
   virtual void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& A,
        ModelPart& rModelPart)
   {
      //filling with zero the matrix (creating the structure)
      Timer::Start("MatrixStructure");

      // Getting the elements from the model
      const int nelements = static_cast<int>(rModelPart.Elements().size());

      // Getting the array of the conditions
      const int nconditions = static_cast<int>(rModelPart.Conditions().size());

      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
      ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
      ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

      const std::size_t equation_size = BaseType::mEquationSystemSize;

#ifdef USE_GOOGLE_HASH
      std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
      const std::size_t empty_key = 2 * equation_size + 10;
#else
      std::vector<std::unordered_set<std::size_t> > indices(equation_size);
#endif

#pragma omp parallel for firstprivate(equation_size)
      for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
      {
#ifdef USE_GOOGLE_HASH
         indices[iii].set_empty_key(empty_key);
#else
         indices[iii].reserve(40);
#endif
      }

      Element::EquationIdVectorType ids(3, 0);

#pragma omp parallel for firstprivate(nelements, ids)
      for (int iii = 0; iii<nelements; iii++)
      {
         typename ElementsContainerType::iterator i_element = el_begin + iii;
         pScheme->EquationId( *(i_element.base()), ids, CurrentProcessInfo);

         for (std::size_t i = 0; i < ids.size(); i++)
         {
            if (ids[i] < BaseType::mEquationSystemSize)
            {
#ifdef _OPENMP
                                    omp_set_lock(&mLockArray[ids[i]]);
#endif
               auto& row_indices = indices[ids[i]];
               for (auto it = ids.begin(); it != ids.end(); it++)
               {
                  if (*it < BaseType::mEquationSystemSize)
                     row_indices.insert(*it);
               }
#ifdef _OPENMP
               omp_unset_lock(&mLockArray[ids[i]]);
#endif
            }
         }

      }

#pragma omp parallel for firstprivate(nconditions, ids)
      for (int iii = 0; iii<nconditions; iii++)
      {
         typename ConditionsArrayType::iterator i_condition = cond_begin + iii;
         pScheme->Condition_EquationId( *(i_condition.base()) , ids, CurrentProcessInfo);
         for (std::size_t i = 0; i < ids.size(); i++)
         {
            if (ids[i] < BaseType::mEquationSystemSize)
            {
#ifdef _OPENMP
               omp_set_lock(&mLockArray[ids[i]]);
#endif
               auto& row_indices = indices[ids[i]];
               for (auto it = ids.begin(); it != ids.end(); it++)
               {
                  if (*it < BaseType::mEquationSystemSize)
                     row_indices.insert(*it);
               }
#ifdef _OPENMP
               omp_unset_lock(&mLockArray[ids[i]]);
#endif
            }
         }
      }

      //count the row sizes
      unsigned int nnz = 0;
      for (unsigned int i = 0; i < indices.size(); i++)
         nnz += indices[i].size();

      A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

      double* Avalues = A.value_data().begin();
      std::size_t* Arow_indices = A.index1_data().begin();
      std::size_t* Acol_indices = A.index2_data().begin();

      //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
      Arow_indices[0] = 0;
      for (int i = 0; i < static_cast<int>(A.size1()); i++)
         Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();



#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(A.size1()); i++)
      {
         const unsigned int row_begin = Arow_indices[i];
         const unsigned int row_end = Arow_indices[i + 1];
         unsigned int k = row_begin;
         for (auto it = indices[i].begin(); it != indices[i].end(); it++)
         {
            Acol_indices[k] = *it;
            Avalues[k] = 0.0;
            k++;
         }

         std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

      }

      A.set_filled(indices.size() + 1, nnz);

      Timer::Stop("MatrixStructure");
   }

//    virtual void ConstructMatrixStructure(
//        TSystemMatrixType& A,
//        ElementsContainerType& rElements,
//        ConditionsArrayType& rConditions,
//        ProcessInfo& CurrentProcessInfo)
//    {
//
//        std::size_t equation_size = A.size1();
//        std::vector<std::vector<std::size_t> > indices(equation_size);
//        //            std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));
//
//        Element::EquationIdVectorType ids(3, 0);
//        for (typename ElementsContainerType::iterator i_element = rElements.begin(); i_element != rElements.end(); i_element++)
//        {
//            (i_element)->EquationIdVector(ids, CurrentProcessInfo);
//
//            for (std::size_t i = 0; i < ids.size(); i++)
//                if (ids[i] < equation_size)
//                {
//                    std::vector<std::size_t>& row_indices = indices[ids[i]];
//                    for (std::size_t j = 0; j < ids.size(); j++)
//                        if (ids[j] < equation_size)
//                        {
//                            AddUnique(row_indices, ids[j]);
//                            //indices[ids[i]].push_back(ids[j]);
//                        }
//                }
//
//        }
//
//        for (typename ConditionsArrayType::iterator i_condition = rConditions.begin(); i_condition != rConditions.end(); i_condition++)
//        {
//            (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
//            for (std::size_t i = 0; i < ids.size(); i++)
//                if (ids[i] < equation_size)
//                {
//                    std::vector<std::size_t>& row_indices = indices[ids[i]];
//                    for (std::size_t j = 0; j < ids.size(); j++)
//                        if (ids[j] < equation_size)
//                        {
//                            AddUnique(row_indices, ids[j]);
//                            //   indices[ids[i]].push_back(ids[j]);
//                        }
//                }
//        }
//
//        //allocating the memory needed
//        int data_size = 0;
//        for (std::size_t i = 0; i < indices.size(); i++)
//        {
//            data_size += indices[i].size();
//        }
//        A.reserve(data_size, false);
//
//        //filling with zero the matrix (creating the structure)
//        Timer::Start("MatrixStructure");
//#ifndef _OPENMP
//        for (std::size_t i = 0; i < indices.size(); i++)
//        {
//            std::vector<std::size_t>& row_indices = indices[i];
//            std::sort(row_indices.begin(), row_indices.end());
//
//            for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
//            {
//                A.push_back(i, *it, 0.00);
//            }
//            row_indices.clear();
//        }
//#else
//        int number_of_threads = omp_get_max_threads();
//        vector<unsigned int> matrix_partition;
//        CreatePartition(number_of_threads, indices.size(), matrix_partition);
//        if (this->GetEchoLevel() > 2)
//        {
//            KRATOS_WATCH(matrix_partition);
//        }
//        for (int k = 0; k < number_of_threads; k++)
//        {
//            #pragma omp parallel
//            if (omp_get_thread_num() == k)
//            {
//                for (std::size_t i = matrix_partition[k]; i < matrix_partition[k + 1]; i++)
//                {
//                    std::vector<std::size_t>& row_indices = indices[i];
//                    std::sort(row_indices.begin(), row_indices.end());
//
//                    for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
//                    {
//                        A.push_back(i, *it, 0.00);
//                    }
//                    row_indices.clear();
//                }
//            }
//        }
//#endif
//        Timer::Stop("MatrixStructure");
//    }

    //**************************************************************************

    void AssembleLHS(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                }
            }
        }
    }

    inline void AssembleRowContribution(TSystemMatrixType& A, const Matrix& Alocal, const std::size_t i, const std::size_t i_local, const Element::EquationIdVectorType& EquationId)
    {
        double* values_vector = A.value_data().begin();
        std::size_t* index1_vector = A.index1_data().begin();
        std::size_t* index2_vector = A.index2_data().begin();

        const std::size_t left_limit = index1_vector[i];

        // Find the first entry
        std::size_t last_pos, last_found;
        std::size_t counter = 0;
        for(std::size_t j=0; j < EquationId.size(); ++j) {
            ++counter;
            const std::size_t j_global = EquationId[j];
            if (j_global < BaseType::mEquationSystemSize) {
                last_pos = ForwardFind(j_global,left_limit,index2_vector);
                last_found = j_global;
                break;
            }
        }

        if (counter < EquationId.size()) {
#ifndef USE_LOCKS_IN_ASSEMBLY
            double& r_a = values_vector[last_pos];
            const double& v_a = Alocal(i_local,counter - 1);
            #pragma omp atomic
            r_a +=  v_a;
#else
            values_vector[last_pos] += Alocal(i_local,counter - 1);
#endif
            // Now find all of the other entries
            std::size_t pos = 0;
            for(std::size_t j=counter; j<EquationId.size(); ++j) {
                std::size_t id_to_find = EquationId[j];
                if (id_to_find < BaseType::mEquationSystemSize) {
                    if(id_to_find > last_found)
                        pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
                    else if(id_to_find < last_found)
                        pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
                    else
                        pos = last_pos;

#ifndef USE_LOCKS_IN_ASSEMBLY
                    double& r = values_vector[pos];
                    const double& v = Alocal(i_local,j);
                    #pragma omp atomic
                    r +=  v;
#else
                    values_vector[pos] += Alocal(i_local,j);
#endif
                    last_found = id_to_find;
                    last_pos = pos;
                }
            }
        }
    }

    inline std::size_t ForwardFind(const std::size_t id_to_find,
                                   const std::size_t start,
                                   const std::size_t* index_vector)
    {
        std::size_t pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline std::size_t BackwardFind(const std::size_t id_to_find,
                                    const std::size_t start,
                                    const std::size_t* index_vector)
    {
        std::size_t pos = start;
        while(id_to_find != index_vector[pos]) pos--;
        return pos;
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

    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while (i != endit && (*i) != candidate)
        {
            i++;
        }
        if (i == endit)
        {
            v.push_back(candidate);
        }

    }

    void AssembleRHS(
        TSystemVectorType& b,
        const LocalSystemVectorType& RHS_Contribution,
        const Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = RHS_Contribution.size();

        if (BaseType::mCalculateReactionsFlag == false)
        {
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
            {
                const unsigned int i_global = EquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) //free dof
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = b[i_global];
                    const double& rhs_value = RHS_Contribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
            }
        }
        else
        {
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
            {
                const unsigned int i_global = EquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) //free dof
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = b[i_global];
                    const double& rhs_value = RHS_Contribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
                else //fixed dof
                {
                    double& b_value = ReactionsVector[i_global - BaseType::mEquationSystemSize];
                    const double& rhs_value = RHS_Contribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
            }
        }
    }

    //**************************************************************************

    void AssembleLHS_CompleteOnFreeRows(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    int j_global = EquationId[j_local];
                    A(i_global, j_global) += LHS_Contribution(i_local, j_local);
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
