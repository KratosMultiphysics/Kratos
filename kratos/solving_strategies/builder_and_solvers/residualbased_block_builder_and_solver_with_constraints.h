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
#include "utilities/multipoint_constraint.h"

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

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


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
    ModelPart& r_model_part
    ) override
    {
        BaseType::SetUpSystem(r_model_part);
        ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();
        Constraints_ExecuteBeforeBuilding(r_model_part, CurrentProcessInfo);
    }

    /** Destructor.
	 */
    virtual ~ResidualBasedBlockBuilderAndSolverWithConstraints()
    {
    }

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();

        Constraints_ExecuteBeforeBuilding(r_model_part, CurrentProcessInfo);

        Build(pScheme, r_model_part, A, b);

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

        Constraints_ExecuteAfterBuilding(r_model_part, CurrentProcessInfo);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        Constraints_ExecuteBeforeSolving(A, Dx, b, r_model_part, CurrentProcessInfo);

        this->SystemSolveWithPhysics(A, Dx, b, r_model_part);

        Timer::Stop("Solve");
        double stop_solve = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "system solve time: " << stop_solve - start_solve << std::endl;

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "after the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        Constraints_ExecuteAfterSolving(A, Dx, b, r_model_part, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    // This is modified to include the MPC information.
    //
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        if (!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        const int nelements = static_cast<int>(r_model_part.Elements().size());

        //getting the array of the conditions
        const int nconditions = static_cast<int>(r_model_part.Conditions().size());

        ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = r_model_part.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

#pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId)
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
                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    // Modifying the local contributions for MPC
                    this->Element_ApplyConstraints(*(*(it.base())), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
//assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
#else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }

//#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
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
                    pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Modifying the local contributions for MPC
                    this->Condition_ApplyConstraints(*(*(it.base())), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

//assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
#else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif

                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }
        }

        // equation_ids.close();
        // KRATOS_THROW_ERROR(std::logic_error,"i want to stop here :-D","")

        double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
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
        ProcessInfo &CurrentProcessInfo) override
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
            (i_element)->EquationIdVector(ids, CurrentProcessInfo);

            // Modifying the equation IDs of this element to suit MPCs
            this->Element_ModifyEquationIdsForConstraints(*(*(i_element.base())), ids, CurrentProcessInfo);

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
            (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
            // Modifying the equation IDs of this element to suit MPCs
            this->Condition_ModifyEquationIdsForConstraints(*(*(i_condition.base())), ids, CurrentProcessInfo);

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

  protected:
    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*
     * This function changes/extends the element LHS and RHS to apply constraints
     */

    void Element_ApplyConstraints(Element &rCurrentElement,
                                  LocalSystemMatrixType &LHS_Contribution,
                                  LocalSystemVectorType &RHS_Contribution,
                                  Element::EquationIdVectorType &EquationId,
                                  ProcessInfo &CurrentProcessInfo)
    {

        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->Element_ApplyConstraints(rCurrentElement, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                }
            }
        }

    } // End of function

    /*
     * This function changes/extends the condition LHS and RHS to apply constraints
     */
    void Condition_ApplyConstraints(Condition &rCurrentElement,
                                    LocalSystemMatrixType &LHS_Contribution,
                                    LocalSystemVectorType &RHS_Contribution,
                                    Element::EquationIdVectorType &EquationId,
                                    ProcessInfo &CurrentProcessInfo)
    {
        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->Condition_ApplyConstraints(rCurrentElement, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                }
            }
        }

    } // End of the function

    /*
     * This function Formulates the constraint data in equation ID terms for formulate the sparsity pattern of the sparse matrix
     */
    void Element_ModifyEquationIdsForConstraints(Element &rCurrentElement,
                                                 Element::EquationIdVectorType &EquationId,
                                                 ProcessInfo &CurrentProcessInfo)
    {
        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->Element_ModifyEquationIdsForConstraints(rCurrentElement, EquationId, CurrentProcessInfo);
                }
            }
        }
    }

    /*
     * This function Formulates the constraint data in equation ID terms for formulate the sparsity pattern of the sparse matrix
     */
    void Condition_ModifyEquationIdsForConstraints(Condition &rCurrentCondition,
                                                   Condition::EquationIdVectorType &EquationId,
                                                   ProcessInfo &CurrentProcessInfo)
    {
        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->Condition_ModifyEquationIdsForConstraints(rCurrentCondition, EquationId, CurrentProcessInfo);
                }
            }
        }
    }

    void Constraints_ExecuteBeforeBuilding(ModelPart &model_part, ProcessInfo &CurrentProcessInfo)
    {

        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->ExecuteBeforeBuilding(model_part.Nodes());
                }
            }
        }
    }

    void Constraints_ExecuteAfterBuilding(ModelPart &model_part, ProcessInfo &CurrentProcessInfo)
    {

        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->ExecuteAfterBuilding(model_part.Nodes());
                }
            }
        }
    }

    void Constraints_ExecuteBeforeSolving(TSystemMatrixType &A,
                                          TSystemVectorType &Dx,
                                          TSystemVectorType &b,
                                          ModelPart &r_model_part,
                                          ProcessInfo &CurrentProcessInfo)
    {

        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->ExecuteBeforeSolving(A, Dx, b);
                }
            }
        }
    }

    void Constraints_ExecuteAfterSolving(TSystemMatrixType &A,
                                         TSystemVectorType &Dx,
                                         TSystemVectorType &b,
                                         ModelPart &r_model_part,
                                         ProcessInfo &CurrentProcessInfo)
    {

        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto &constraint : (*constraintVector))
            {
                if (constraint->IsActive())
                {
                    constraint->ExecuteAfterSolving(A, Dx, b);
                }
            }
        }
    }
};
}

#endif /* KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_ */
