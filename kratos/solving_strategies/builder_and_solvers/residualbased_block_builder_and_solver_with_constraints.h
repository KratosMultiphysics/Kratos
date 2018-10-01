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
#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS)
#define KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS

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
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/master_slave_constraint.h"
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
 * @class ResidualBasedBlockBuilderAndSolverWithConstraints
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similar to the calculation of the total residual
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
class ResidualBasedBlockBuilderAndSolverWithConstraints
    : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithConstraints);

    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
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
    typedef MasterSlaveConstraint MasterSlaveConstraintType;
    typedef typename MasterSlaveConstraint::Pointer MasterSlaveConstraintPointerType;
    typedef Internals::AuxiliaryGlobalMasterSlaveConstraint AuxiliaryGlobalMasterSlaveConstraintType;
    typedef Internals::GlobalMasterSlaveRelationContainerType GlobalMasterSlaveRelationContainerType;
    typedef std::vector<IndexType> EquationIdVectorType;
    typedef std::vector<IndexType> VectorIndexType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef Vector VectorType;
    typedef Internals::ConstraintImposer<TSparseSpace, TDenseSpace, TLinearSolver> ConstraintImposerType;
    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverWithConstraints() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetUpSystem(
        ModelPart &rModelPart) override
    {
        BaseType::SetUpSystem(rModelPart);
        if(rModelPart.NumberOfMasterSlaveConstraints() > 0)
        {
            FormulateGlobalMasterSlaveRelations(rModelPart);
        }
    }

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
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &b) override
    {
        if(mGlobalMasterSlaveConstraints.size() > 0)
            BuildWithConstraints(pScheme, rModelPart, A, b);
        else
            BaseType::Build(pScheme, rModelPart, A, b);
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
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        if(mGlobalMasterSlaveConstraints.size() > 0)
            BuildAndSolveWithConstraints(pScheme, rModelPart, A, Dx, b);
        else
            BaseType::BuildAndSolve(pScheme, rModelPart, A, Dx, b);
    }

    void InitializeSolutionStep(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, A, Dx, b);
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
#pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; k++)
        {
            auto it = constraints_begin + k;
            it->InitializeSolutionStep(rModelPart.GetProcessInfo()); // Here each constraint constructs and stores its T and C matrices. Also its equation ids.
        }

        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints failed to initialize solution step.")
    }


    void FinalizeSolutionStep(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        BaseType::FinalizeSolutionStep(rModelPart, A, Dx, b);
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
#pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; k++)
        {
            auto it = constraints_begin + k;
            it->FinalizeSolutionStep(rModelPart.GetProcessInfo());
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints failed to finalize solution step.")
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

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
    virtual void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType &A,
        ModelPart &rModelPart) override
    {
        if(mGlobalMasterSlaveConstraints.size() > 0)
            ConstructMatrixStructureWithConstraints(pScheme, A, rModelPart);
        else
            BaseType::ConstructMatrixStructure(pScheme, A, rModelPart);
    }

    void BuildAndSolveWithConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b)
    {
        KRATOS_TRY

        const double start_update_constraints = OpenMPUtils::GetCurrentTime();
        this->UpdateConstraintsForBuilding(rModelPart);
        const double stop_update_constraints = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Constraints update time : " << stop_update_constraints - start_update_constraints << std::endl;


        Timer::Start("Build");

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) << "Before the solution of the system"
                                                                                                         << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");
        this->SystemSolveWithPhysics(A, Dx, b, rModelPart);
        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();

        const double start_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        ReconstructSlaveSolutionAfterSolve(rModelPart, A, Dx, b);
        const double stop_reconstruct_slaves = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Reconstruct slaves time: " << stop_reconstruct_slaves - start_reconstruct_slaves << std::endl;


        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() == 3)) << "After the solution of the system"
                                                                                                         << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        KRATOS_CATCH("")
    }

    /*
    * This function is exactly same as the ConstructMatrixStructure() function in base class except that the function
    * has the call to ApplyConstraints function call once the element and conditions compute their equation ids
    */
    virtual void ConstructMatrixStructureWithConstraints(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType &A,
        ModelPart &rModelPart)
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");
        ConstraintImposerType constraint_imposer(mGlobalMasterSlaveConstraints);

        const std::size_t equation_size = BaseType::mEquationSystemSize;

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

        const int nelements = static_cast<int>(rModelPart.Elements().size());
#pragma omp parallel for firstprivate(nelements, ids, constraint_imposer)
        for (int iii = 0; iii < nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = rModelPart.Elements().begin() + iii;
            pScheme->EquationId(*(i_element.base()), ids, rModelPart.GetProcessInfo());
            constraint_imposer.template ApplyConstraints<Element>(*i_element, ids, rModelPart.GetProcessInfo());
            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&BaseType::mlock_array[ids[i]]);
#endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
                omp_unset_lock(&BaseType::mlock_array[ids[i]]);
#endif
            }
        }

        const int nconditions = static_cast<int>(rModelPart.Conditions().size());
#pragma omp parallel for firstprivate(nconditions, ids, constraint_imposer)
        for (int iii = 0; iii < nconditions; iii++)
        {
            typename ConditionsArrayType::iterator i_condition = rModelPart.Conditions().begin() + iii;
            pScheme->Condition_EquationId(*(i_condition.base()), ids, rModelPart.GetProcessInfo());
            constraint_imposer.template ApplyConstraints<Condition>(*i_condition, ids, rModelPart.GetProcessInfo());
            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&BaseType::mlock_array[ids[i]]);
#endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
                omp_unset_lock(&BaseType::mlock_array[ids[i]]);
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

    /*
    * This function is exactly same as the Build() function in base class except that the function
    * has the call to ApplyConstraints function call once the LHS or RHS are computed by elements and conditions
    */
    void BuildWithConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &b)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;
        ConstraintImposerType constraint_imposer(mGlobalMasterSlaveConstraints);

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        const ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        const ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

#pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, constraint_imposer)
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
                    constraint_imposer.template ApplyConstraints<Element>(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

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
                    constraint_imposer.template ApplyConstraints<Condition>(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

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

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
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

    /**
     * @brief   this method condenses the MasterSlaveConstraints which are added on the rModelPart
     *          into objects of AuxilaryGlobalMasterSlaveRelation. One unique object for each unique slave.
     *          these will be used in the ApplyConstraints functions later on.
     * @param   rModelPart The model part of the problem to solve
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

        const ModelPart::MasterSlaveConstraintContainerType::iterator constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        //contributions to the system
        LocalSystemMatrixType relation_matrix;
        LocalSystemVectorType constant_vector;
        EquationIdVectorType slave_equation_ids;
        EquationIdVectorType master_equation_ids;

#pragma omp parallel for schedule(guided, 512) private (relation_matrix, constant_vector, slave_equation_ids, master_equation_ids)
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
                it->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
                //calculate constraint's T and b matrices
                it->CalculateLocalSystem(relation_matrix, constant_vector, r_current_process_info);
                //assemble the Constraint contribution
                AssembleSlaves(slave_equation_ids, master_equation_ids, relation_matrix);
            }
        }
        const double stop_formulate = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Formulate global constraints time: " << stop_formulate - start_formulate << std::endl;

        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::FormulateGlobalMasterSlaveRelations failed ..");
    }


    /**
     * @brief   this method assembles the given list of slave equation ids and corresponding master equation ids
     * @param   rSlaveEquationIdVector list of slave equation ids to be assembled.
     * @param   rMasterEquationIdVector list of corresponding master equation ids.
     * @param   rRelationMatrix the matrix relating the rSlaveEquationIdVector with rMasterEquationIdVector
     */
    void AssembleSlaves(EquationIdVectorType& rSlaveEquationIdVector, EquationIdVectorType& rMasterEquationIdVector, LocalSystemMatrixType& rRelationMatrix)
    {
        KRATOS_TRY
        int slave_count = 0;
        for (auto slave_equation_id : rSlaveEquationIdVector)
        {
            /*int master_count = 0;
            auto global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            if (global_constraint == mGlobalMasterSlaveConstraints.end())
            {
                mGlobalMasterSlaveConstraints.insert(mGlobalMasterSlaveConstraints.begin(), Kratos::make_shared<AuxiliaryGlobalMasterSlaveConstraintType>(slave_equation_id));
                global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            }*/
            int master_count = 0;
            auto global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            if (global_constraint == mGlobalMasterSlaveConstraints.end())
            {
                mGlobalMasterSlaveConstraints[slave_equation_id] = Kratos::make_unique<AuxiliaryGlobalMasterSlaveConstraintType>(slave_equation_id);
                global_constraint = mGlobalMasterSlaveConstraints.find(slave_equation_id);
            }

            for (auto master_equation_id : rMasterEquationIdVector)
            {
                global_constraint->second->AddMaster(master_equation_id, rRelationMatrix(slave_count, master_count));
                master_count++;
            }
            slave_count++;
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::AssembleSlaves failed ...");
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

		#pragma omp parallel for schedule(guided, 512)
        for (int i_constraints = 0; i_constraints < number_of_constraints; ++i_constraints)
        {
            //GlobalMasterSlaveRelationContainerType::iterator it = constraints_begin + i_constraints;
            GlobalMasterSlaveRelationContainerType::iterator it = constraints_begin;
            std::advance(it, i_constraints);

            (it->second)->Reset();
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::ResetConstraintRelations failed to reset constraint relations..");
    }

    /**
     * @brief   this method uses the MasterSlaveConstraints objects in rModelPart to reconstruct the LHS and RHS values
     *          of the AuxilaryGlobalMasterSlaveRelation objects. That is the value of Slave as LHS and the T*M+C as RHS value
     * @param   rModelPart The model part of the problem to solve
     */
    void UpdateConstraintsForBuilding(ModelPart& rModelPart)
    {
        KRATOS_TRY
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
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::UpdateConstraintsForBuilding failed ..");
    }


    /**
     * @brief   this method uses the MasterSlaveConstraints objects in rModelPart to reconstruct the LHS and RHS values
     *          of the individual AuxilaryGlobalMasterSlaveRelation object. That is the value of Slave as LHS and the T*M+C as RHS value
     * @param   rMasterSlaveConstraint The MasterSlaveConstraint which is to be updated
     */
    void UpdateMasterSlaveConstraint(ModelPart::MasterSlaveConstraintType& rMasterSlaveConstraint, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        //contributions to the system
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
        for (auto &slave_dof : slave_dofs_vector)
        {
            double slave_value_calc = 0.0;
            for (IndexType master_index = 0; master_index < master_dofs_vector.size(); master_index++)
            {
                slave_value_calc += master_dofs_vector[master_index]->GetSolutionStepValue() * relation_matrix(slave_index, master_index);
            }
            slave_value_calc += constant_vector[slave_index];
            auto global_constraint = mGlobalMasterSlaveConstraints.find(slave_dof->EquationId());
            global_constraint->second->SetLeftHandSide( slave_dof->GetSolutionStepValue() );
            global_constraint->second->UpdateRightHandSide(slave_value_calc);
            slave_index++;
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::UpdateMasterSlaveConstraint failed ..");
    }

    /**
     * @brief This method reconstructs the slave solution after Solving.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void ReconstructSlaveSolutionAfterSolve(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb)
    {
        KRATOS_TRY
        const int number_of_constraints = static_cast<int>(mGlobalMasterSlaveConstraints.size());
        // Getting the beginning iterator

        const GlobalMasterSlaveRelationContainerType::iterator constraints_begin = mGlobalMasterSlaveConstraints.begin();
        //contributions to the system
        VectorType master_weights_vector;
        double constant = 0.0;

        IndexType slave_equation_id = 0;
        EquationIdVectorType master_equation_ids = EquationIdVectorType(0);

#pragma omp parallel for schedule(guided, 512) firstprivate(slave_equation_id, master_equation_ids, master_weights_vector, constant)
        for (int i_constraints = 0; i_constraints < number_of_constraints; i_constraints++)
        {
            //GlobalMasterSlaveRelationContainerType::iterator it = constraints_begin + i_constraints;
            GlobalMasterSlaveRelationContainerType::iterator it = constraints_begin;
            std::advance(it, i_constraints);

            double slave_dx_value = 0.0;
            //get the equation Ids of the constraint
            (it->second)->EquationIdsVector(slave_equation_id, master_equation_ids);
            //calculate constraint's T and b matrices
            (it->second)->CalculateLocalSystem(master_weights_vector, constant);
            int master_index = 0;
            for (auto &master_equation_id : master_equation_ids)
            {
                slave_dx_value += TSparseSpace::GetValue(rDx, master_equation_id) * master_weights_vector(master_index);
                master_index++;
            }
            slave_dx_value += constant;

            rDx[slave_equation_id] = slave_dx_value; // this access is always unique for an object so no need of special care for openmp
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::ReconstructSlaveSolutionAfterSolve failed ..");
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

}; /* Class ResidualBasedBlockBuilderAndSolverWithConstraints */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
