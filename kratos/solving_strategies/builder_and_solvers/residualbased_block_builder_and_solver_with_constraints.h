//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala / Navaneeth K Narayanan
//
//

#ifndef KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS_H_
#define KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS_H_

/* System includes */

#include "utilities/openmp_utils.h"
#include <unordered_set>
#include <algorithm>
/* External includes */
#include "boost/smart_ptr.hpp"

#include "utilities/timer.h"
#include "utilities/constraint.h"

/* Project includes */
#include "includes/define.h"
#include "residualbased_block_builder_and_solver.h"
#include "includes/model_part.h"

// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#endif
// #define USE_LOCKS_IN_ASSEMBLY
// #include <iostream>
// #include <fstream>

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */

/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations as in 
ResidualBasedBlockBuilderAndSolver 

see 
kratos/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h

In addition this allows enforcing constraints of the form S = T*M + C on the already built 
linear system of equations Ax=b

in the above
S = Slave degrees of freedom
M = Master degrees of freedom
T = Relation between master and slave
C = Constant matrix. 

 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ResidualBasedBlockBuilderAndSolverWithConstraints
    : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithConstraints);

    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

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
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;

    typedef Node<3> NodeType;
    typedef typename ModelPart::NodesContainerType NodesContainerType;
    typedef typename Constraint<TSparseSpace, TDenseSpace>::Pointer ConstraintPointerType;
    typedef boost::shared_ptr<std::vector<ConstraintPointerType>> ConstraintSharedPointerVectorType;

    typedef ProcessInfo ProcessInfoType;

    /*@} */
    /**@name Life Cycle
	 */
    /*@{ */

    /** Constructor.
	 */
    ResidualBasedBlockBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
    }

    void SetUpSystem(
        ModelPart &rModelPart) override
    {
        ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
        mHasConstraints = rCurrentProcessInfo.Has(CONSTRAINTS_CONTAINER);
        BaseType::SetUpSystem(rModelPart);
        if (mHasConstraints)
        {
            mpConstraintVector = rCurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            Constraints_SetUp(rModelPart, mpConstraintVector);
        }
    }

    /** Destructor.
	 */
    virtual ~ResidualBasedBlockBuilderAndSolverWithConstraints()
    {
    }

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        if (mHasConstraints)
        {
            Constraints_ExecuteBeforeBuilding(rModelPart, mpConstraintVector);
            Build(pScheme, rModelPart, A, b);
        }
        else
        {
            BaseType::Build(pScheme, rModelPart, A, b);
        }

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        if (mHasConstraints)
            Constraints_ExecuteAfterBuilding(rModelPart, mpConstraintVector);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");
        if (mHasConstraints)
            Constraints_ExecuteBeforeSolving(A, Dx, b, rModelPart, mpConstraintVector);

        this->SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        double stop_solve = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "system solve time: " << stop_solve - start_solve << std::endl;

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "after the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }
        if (mHasConstraints)
            Constraints_ExecuteAfterSolving(A, Dx, b, rModelPart, mpConstraintVector);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    // This is modified to include the MPC information.
    //
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        if (!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        //getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType rEquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

#pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, rEquationId)
        {
#pragma omp for schedule(guided, 512) nowait
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
                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);
                    // Modifying the local contributions for MPC
                    this->ApplyConstraints(*(*(it.base())), LHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo, mpConstraintVector);
//assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, rEquationId, BaseType::mlock_array);
#else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, rEquationId);
#endif
                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }

//#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, rEquationId ) schedule(dynamic, 1024)
#pragma omp for schedule(guided, 512)
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
                    pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);

                    // Modifying the local contributions for MPC
                    this->ApplyConstraints(*(*(it.base())), LHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo, mpConstraintVector);

//assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, rEquationId, BaseType::mlock_array);
#else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, rEquationId);
#endif

                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }
        }

        // equation_ids.close();
        // KRATOS_THROW_ERROR(std::logic_error,"i want to stop here :-D","")

        double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        if (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH("finished parallel building");
        }

        KRATOS_CATCH("")
    }

  protected:
    void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType &A,
        ElementsContainerType &rElements,
        ConditionsArrayType &rConditions,
        ProcessInfo &rCurrentProcessInfo) override
    {
        if (mHasConstraints)
            Constraints_ConstructMatrixStructure(pScheme, A, rElements, rConditions, rCurrentProcessInfo);
        else
            BaseType::ConstructMatrixStructure(pScheme, A, rElements, rConditions, rCurrentProcessInfo);
    }

    void Constraints_ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType &A,
        ElementsContainerType &rElements,
        ConditionsArrayType &rConditions,
        ProcessInfo &rCurrentProcessInfo)
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        const std::size_t equation_size = BaseType::mDofSet.size();

#ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<std::size_t>> indices(equation_size);
        const std::size_t empty_key = 2 * equation_size + 10;
#else
        std::vector<std::unordered_set<std::size_t>> indices(equation_size);
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

        const int nelements = static_cast<int>(rElements.size());
#pragma omp parallel for firstprivate(nelements, ids)
        for (int iii = 0; iii < nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
            (i_element)->EquationIdVector(ids, rCurrentProcessInfo);

            // Modifying the equation IDs of this element to suit MPCs
            this->ModifyEquationIdsForConstraints(*(*(i_element.base())), ids, rCurrentProcessInfo, mpConstraintVector);

            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&(BaseType::mlock_array[ids[i]]));
#endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
                omp_unset_lock(&(BaseType::mlock_array[ids[i]]));
#endif
            }
        }
        const int nconditions = static_cast<int>(rConditions.size());
#pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii < nconditions; iii++)
        {
            typename ConditionsArrayType::iterator i_condition = rConditions.begin() + iii;
            (i_condition)->EquationIdVector(ids, rCurrentProcessInfo);
            // Modifying the equation IDs of this element to suit MPCs
            this->ModifyEquationIdsForConstraints(*(*(i_condition.base())), ids, rCurrentProcessInfo, mpConstraintVector);

            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&(BaseType::mlock_array[ids[i]]));
#endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
                omp_unset_lock(&(BaseType::mlock_array[ids[i]]));
#endif
            }
        }
        //count the row sizes
        unsigned int nnz = 0;
        for (unsigned int i = 0; i < indices.size(); i++)
            nnz += indices[i].size();

        A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

        double *Avalues = A.value_data().begin();
        std::size_t *Arow_indices = A.index1_data().begin();
        std::size_t *Acol_indices = A.index2_data().begin();

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

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
        }

        A.set_filled(indices.size() + 1, nnz);
        Timer::Stop("MatrixStructure");
    }

  private:
    /*@} */
    /**@name Private member Variables */
    /*@{ */
    bool mHasConstraints;
    ConstraintSharedPointerVectorType mpConstraintVector;

    /*
     * This function changes/extends the element LHS and RHS to apply constraints
     */

    void ApplyConstraints(Element &rCurrentElement,
                                  LocalSystemMatrixType &LHS_Contribution,
                                  LocalSystemVectorType &RHS_Contribution,
                                  Element::EquationIdVectorType &rEquationId,
                                  ProcessInfo &rCurrentProcessInfo,
                                  ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ApplyConstraints(rCurrentElement, LHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);
            }
        }

    } // End of function

    /*
     * This function changes/extends the condition LHS and RHS to apply constraints
     */
    void ApplyConstraints(Condition &rCurrentCondition,
                                    LocalSystemMatrixType &LHS_Contribution,
                                    LocalSystemVectorType &RHS_Contribution,
                                    Condition::EquationIdVectorType &rEquationId,
                                    ProcessInfo &rCurrentProcessInfo,
                                    ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ApplyConstraints(rCurrentCondition, LHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);
            }
        }

    } // End of the function

    /*
     * This function Formulates the constraint data in equation ID terms for formulate the sparsity pattern of the sparse matrix
     */
    void ModifyEquationIdsForConstraints(Element &rCurrentElement,
                                                 Element::EquationIdVectorType &rEquationId,
                                                 ProcessInfo &rCurrentProcessInfo,
                                                 ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ModifyEquationIdsForConstraints(rCurrentElement, rEquationId, rCurrentProcessInfo);
            }
        }
    }

    /*
     * This function Formulates the constraint data in equation ID terms for formulate the sparsity pattern of the sparse matrix
     */
    void ModifyEquationIdsForConstraints(Condition &rCurrentCondition,
                                                   Condition::EquationIdVectorType &rEquationId,
                                                   ProcessInfo &rCurrentProcessInfo,
                                                   ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ModifyEquationIdsForConstraints(rCurrentCondition, rEquationId, rCurrentProcessInfo);
            }
        }
    }

    void Constraints_ExecuteBeforeBuilding(ModelPart &rModelPart, ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ExecuteBeforeBuilding(rModelPart.Nodes());
            }
        }
    }

    void Constraints_ExecuteAfterBuilding(ModelPart &rModelPart, ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ExecuteAfterBuilding(rModelPart.Nodes());
            }
        }
    }

    void Constraints_SetUp(ModelPart &rModelPart, ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->SetUp(rModelPart.Nodes());
            }
        }
    }

    void Constraints_ExecuteBeforeSolving(TSystemMatrixType &A,
                                          TSystemVectorType &Dx,
                                          TSystemVectorType &b,
                                          ModelPart &rModelPart,
                                          ConstraintSharedPointerVectorType &rConstraintVector)

    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ExecuteBeforeSolving(A, Dx, b);
            }
        }
    }

    void Constraints_ExecuteAfterSolving(TSystemMatrixType &A,
                                         TSystemVectorType &Dx,
                                         TSystemVectorType &b,
                                         ModelPart &rModelPart,
                                         ConstraintSharedPointerVectorType &rConstraintVector)
    {
        for (auto &constraint : (*rConstraintVector))
        {
            if (constraint->IsActive())
            {
                constraint->ExecuteAfterSolving(A, Dx, b);
            }
        }
    }
};
}

#endif /* KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_ */
