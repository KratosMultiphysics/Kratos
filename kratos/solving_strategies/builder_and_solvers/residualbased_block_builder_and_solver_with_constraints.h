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
#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS )
#define  KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS


/* System includes */
#include <unordered_set>
// #include <iostream>
// #include <fstream>

/* External includes */
// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
    #include "sparsehash/dense_hash_set" //included in external libraries
#endif
// #define USE_LOCKS_IN_ASSEMBLY

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_flags.h"
#include "utilities/master_slave_constraint.h"
#include "containers/constraint_equation.h"

#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "containers/variable_data.h"
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
 * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolverWithConstraints
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
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

    typedef Node<3> NodeType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    typedef MasterSlaveConstraint MasterSlaveConstraintType;
    typedef typename MasterSlaveConstraint::Pointer MasterSlaveConstraintPointerType;

    typedef MasterSlaveRelation MasterSlaveRelationType;

    typedef PointerVectorSet<MasterSlaveRelationType, IndexedObject> MasterSlaveRelationContainerType;
    typedef std::size_t SizeType;
    typedef std::vector<SizeType> EquationIdVectorType;
    typedef std::vector<SizeType> VectorIndexType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ResidualBasedBlockBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
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
        ModelPart &rModelPart
        ) override
    {
        BaseType::SetUpSystem(rModelPart);
        FormulateGlobalMasterSlaveRelations(rModelPart);
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
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
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

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
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
                    ApplyConstraints<Element>(rModelPart, *(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }

            }


            //#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
            #pragma omp for  schedule(guided, 512)
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
                    ApplyConstraints<Condition>(rModelPart, *(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
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

        this->ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        this->SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithConstraints", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        KRATOS_CATCH("")
    }



    //**************************************************************************
    //**************************************************************************

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        ModelPart::MasterSlaveConstraintContainerType::iterator constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp for  schedule(guided, 512)
        for (int k = 0; k < n_constraints; k++)
        {
            ModelPart::MasterSlaveConstraintContainerType::iterator it = constraints_begin + k;
            it->InitializeSolutionStep(); //TODO: Here each constraint constructs and stores its T and C matrices. Also its equation ids.
        }


        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        ModelPart::MasterSlaveConstraintContainerType::iterator constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp for  schedule(guided, 512)
        for (int k = 0; k < n_constraints; k++)
        {
            ModelPart::MasterSlaveConstraintContainerType::iterator it = constraints_begin + k;
            it->FinalizeSolutionStep();
        }

        // TODO: Here the solution to the slaves is reconstructed.
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
#ifdef _OPENMP
    std::vector< omp_lock_t > mlock_array;
#endif
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& A,
        ModelPart& rModelPart) override
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        const std::size_t equation_size = BaseType::mEquationSystemSize;

#ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
        const std::size_t empty_key = 2*equation_size + 10;
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

        const int nelements = static_cast<int>(rModelPart.Elements().size());
        #pragma omp parallel for firstprivate(nelements, ids)
        for(int iii=0; iii<nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = rModelPart.Elements().begin() + iii;
            pScheme->EquationId( *(i_element.base()) , ids, rModelPart.GetProcessInfo());
            ApplyConstraints<Element>(rModelPart, *(i_element.base()) , ids, rModelPart.GetProcessInfo());
            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&mlock_array[ids[i]]);
#endif
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
                omp_unset_lock(&mlock_array[ids[i]]);
#endif
            }

        }

        const int nconditions = static_cast<int>(rModelPart.Conditions().size());
        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii<nconditions; iii++)
        {
            typename ConditionsArrayType::iterator i_condition = rModelPart.Conditions().begin() + iii;
            pScheme->Condition_EquationId( *(i_condition.base()), ids, rModelPart.GetProcessInfo());
            ApplyConstraints<Condition>(rModelPart, *(i_condition.base()), ids, rModelPart.GetProcessInfo());
            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&mlock_array[ids[i]]);
#endif
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
                omp_unset_lock(&mlock_array[ids[i]]);
#endif
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
            Arow_indices[i+1] = Arow_indices[i] + indices[i].size();



        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
        {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i+1];
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

        A.set_filled(indices.size()+1, nnz);

        Timer::Stop("MatrixStructure");
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
    MasterSlaveRelationContainerType mGlobalMasterSlaveRelations; //This can be changed to more efficient implementation later on.

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void FormulateGlobalMasterSlaveRelations(ModelPart& rModelPart)
    {

        // Getting the array of the conditions
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        // Getting the beginning iterator

        ModelPart::MasterSlaveConstraintContainerType::iterator constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        //contributions to the system
        LocalSystemMatrixType relation_matrix = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType constraint_vector = LocalSystemVectorType(0);
        EquationIdVectorType  slave_equation_ids = EquationIdVectorType(0);
        EquationIdVectorType  master_equation_ids = EquationIdVectorType(0);

        for (int i_constraints = 0; i_constraints < number_of_constraints; i_constraints++)
        {
            ModelPart::MasterSlaveConstraintContainerType::iterator it = constraints_begin + i_constraints;

            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool constraint_is_active = true;
            if ((it)->IsDefined(ACTIVE))
                constraint_is_active = (it)->Is(ACTIVE);

            if (constraint_is_active)
            {
                //get the equation Ids of the constraint
                (*it).EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
                //calculate constraint's T and b matrices
                (*it).CalculateLocalSystem(relation_matrix,constraint_vector,r_current_process_info);

                //assemble the Constraint contribution
                int slave_count = 0;
                for (auto& slave_equation_id : slave_equation_ids)
                {
                    int master_count = 0;
                    auto global_constraint = mGlobalMasterSlaveRelations.find(slave_equation_id);
                    if(global_constraint == mGlobalMasterSlaveRelations.end())
                    {
                        mGlobalMasterSlaveRelations.insert(mGlobalMasterSlaveRelations.begin(), Kratos::make_shared<MasterSlaveRelationType>(slave_equation_id));
                        global_constraint = mGlobalMasterSlaveRelations.find(slave_equation_id);
                    }
                    for(auto& master_equation_id : master_equation_ids)
                    {
                        global_constraint->AddMaster(master_equation_id, relation_matrix(slave_count, master_count));
                        master_count++;
                    }
                    slave_count++;
                }
            }
        }
    }


    // This adds the equation IDs of masters of all the slaves correspoining to pCurrentElement to EquationIds
    template<class TContainerType>
    void ApplyConstraints(  ModelPart& rModelPart,
                            typename TContainerType::Pointer pCurrentContainer,
                            typename TContainerType::EquationIdVectorType& rEquationIds,
                            ProcessInfo& rCurrentProcessInfo
                        )
    {
        const SizeType number_of_nodes = pCurrentContainer->GetGeometry().PointsNumber();
        // For each node check if it is a slave or not If it is .. we change the Transformation matrix
        for (IndexType j = 0; j < number_of_nodes; j++) {
            DofsVectorType element_dofs;
            pCurrentContainer->GetDofList(element_dofs, rCurrentProcessInfo);
            SizeType number_dofs_per_node = element_dofs.size() / number_of_nodes;
            if (pCurrentContainer->GetGeometry()[j].Is(SLAVE)) {

                IndexType slave_equation_id;
                IndexType start_position_node_dofs = number_dofs_per_node * (j);

                for (IndexType i = 0; i < number_dofs_per_node; i++) {
                    slave_equation_id = element_dofs[start_position_node_dofs + i]->EquationId(); // consider everything as a slave.
                    // Get the global constraint equation for this slave.
                    auto global_master_slave_constraint = mGlobalMasterSlaveRelations.find(slave_equation_id);
                    if (global_master_slave_constraint != mGlobalMasterSlaveRelations.end()) { // if a equation exists for this slave

                        typename TContainerType::EquationIdVectorType slave_equation_ids;
                        typename TContainerType::EquationIdVectorType master_equation_ids;
                        global_master_slave_constraint->EquationIdVector(slave_equation_ids, master_equation_ids, rCurrentProcessInfo); // get the slave and master equation ids for this slave.
                        for (auto& master_eq_id : master_equation_ids) {
                            // Add the current slaves master eq ids to the equation ids
                            rEquationIds.push_back(master_eq_id);
                        }
                    }
                }
            }
        }
    }

    template<class TContainerType>
    void ApplyConstraints( ModelPart& rModelPart,
                           typename TContainerType::Pointer pCurrentContainer,
                           LocalSystemMatrixType& rLHSContribution,
                           LocalSystemVectorType& rRHSContribution,
                           typename TContainerType::EquationIdVectorType& rEquationIds,
                           ProcessInfo& rCurrentProcessInfo
                        )
    {
        // This function modifies LHS and RHS contributions with T and C matrices of the constraints present in this current pElement

        KRATOS_TRY
        bool slaveFound = false;
        auto& geometry = pCurrentContainer->GetGeometry();
        const SizeType number_of_nodes = geometry.PointsNumber();
        for (IndexType j = 0; j < number_of_nodes; j++) {
            bool node_is_slave = true;
            if (geometry[j].IsDefined(SLAVE))
                node_is_slave = geometry[j].Is(SLAVE);
            if (node_is_slave) { // temporary, will be checked once at the beginning only
                slaveFound = true;
                break;
            }
        }
        // If no slave is found for this container , no need of going on
        if (!slaveFound) {
            return;
        }
        // TODO: make shure that the order of the masters is same in the loops, in that case, the order of all the vectors should be same

        // Saving th original system size
        const SizeType initial_sys_size = rLHSContribution.size1();

        // Formulating which are the internal, slave indices locally.
        VectorIndexType::iterator it;
        VectorIndexType local_slave_index_vector;
        VectorIndexType local_index_vector;
        VectorIndexType local_internal_index_vector;
        VectorIndexType local_master_index_vector;
        int index = 0;
        for (auto& eq_id : rEquationIds)
        {
            auto global_master_slave_constraint = mGlobalMasterSlaveRelations.find(eq_id);
            if (global_master_slave_constraint != mGlobalMasterSlaveRelations.end()) {
                local_slave_index_vector.push_back(index);
            }
            local_index_vector.push_back(index);
            index++;
        }
        std::sort(local_index_vector.begin(), local_index_vector.end());
        std::sort(local_slave_index_vector.begin(), local_slave_index_vector.end());
        std::set_difference(local_index_vector.begin(), local_index_vector.end(), local_slave_index_vector.begin(), local_slave_index_vector.end(), std::back_inserter(local_internal_index_vector));

        // first fill in the rEquationIds using the above function (overloaded one)
        ApplyConstraints<TContainerType>(rModelPart, pCurrentContainer, rEquationIds, rCurrentProcessInfo); // now rEquationIds has all the slave equation ids appended to it.

        // Get number of master indices for this current container
        const SizeType total_number_of_masters = rEquationIds.size() - initial_sys_size;
        for (unsigned int i=initial_sys_size; i<rEquationIds.size(); i++){
            local_master_index_vector.push_back(i);
        }

        // We resize the LHS and RHS contribution with the master sizes
        const SizeType lhs_size_1 = initial_sys_size + total_number_of_masters;
        const SizeType lhs_size_2 = initial_sys_size + total_number_of_masters;
        rLHSContribution.resize(lhs_size_1, lhs_size_2, true); //true for Preserving the data and resizing the matrix
        rRHSContribution.resize(lhs_size_1, true);
        // Making the extra part of matrx
        for (IndexType m = initial_sys_size; m < lhs_size_1; m++) {
            for (IndexType n = 0; n < lhs_size_1; n++) {
                rLHSContribution(m, n) = 0.0;
                rLHSContribution(n, m) = 0.0;
            }
            rRHSContribution(m) = 0.0;
        }


        int slave_count = 0;
        int local_master_index = initial_sys_size;
        MatrixType transformation_matrix;
        VectorType constant_vector;
        std::vector<double> weights_corresponding_to_master;
        std::vector<IndexType> slaves_corresponding_to_masters;

        for (auto& local_slave_index : local_slave_index_vector) { // Loop over all the slaves for this container

            // Get the golobal equation for this constraint
            auto& global_master_slave_constraint = mGlobalMasterSlaveRelations(rEquationIds[local_slave_index]);
            // Get the tranformation matrix and constant_vector from the current slave
            global_master_slave_constraint->CalculateLocalSystem(transformation_matrix, constant_vector, rCurrentProcessInfo);

            typename TContainerType::EquationIdVectorType slave_equation_ids;
            typename TContainerType::EquationIdVectorType master_equation_ids;
            global_master_slave_constraint->EquationIdVector(slave_equation_ids, master_equation_ids, rCurrentProcessInfo);

            KRATOS_ERROR_IF( master_equation_ids.size() != transformation_matrix.size2() )<< " dimensions does not match "<<std::endl;
            KRATOS_ERROR_IF( slave_equation_ids.size() != 1 )<< " slave equation id dimensions does not match to 1 "<<std::endl;

            int slave_master_index = 0;
            for (auto& master_eq_id : master_equation_ids) { // Loop over all the masters the slave has

                const double weight = transformation_matrix(0,slave_master_index);
                const double constant = constant_vector(0);
                for (auto local_internal_index : local_internal_index_vector) {
                    rRHSContribution(local_internal_index) += -rLHSContribution(local_internal_index, local_slave_index) * constant;
                }

                // For K(m,u) and K(u,m)
                for (auto local_internal_index : local_internal_index_vector) { // Loop over all the local equation ids
                    rLHSContribution(local_internal_index, local_master_index) += rLHSContribution(local_internal_index, local_slave_index) * weight;
                    rLHSContribution(local_master_index, local_internal_index) += rLHSContribution(local_slave_index, local_internal_index) * weight;
                } // Loop over all the local equation ids

                // For RHS(m) += A'*LHS(s,s)*B
                for (auto local_slave_index_vector_other : local_slave_index_vector) {
                    MatrixType transformation_matrix_other;
                    VectorType constant_vector_other;
                    auto& global_master_slave_constraint_other = mGlobalMasterSlaveRelations(rEquationIds[local_slave_index_vector_other]);
                    global_master_slave_constraint_other->CalculateLocalSystem(transformation_matrix_other, constant_vector_other, rCurrentProcessInfo);

                    double constant_other = constant_vector_other(0);
                    rRHSContribution(local_master_index) += rLHSContribution(local_slave_index, local_slave_index_vector_other) * weight * constant_other;
                }

                // Changing the RHS side of the equation
                rRHSContribution(local_master_index) += weight * rRHSContribution(local_slave_index);

                weights_corresponding_to_master.push_back(weight);
                slaves_corresponding_to_masters.push_back(local_slave_index);
                local_master_index ++ ;
                slave_master_index++;
            } // Loop over all the masters the slave has

            rRHSContribution(local_slave_index) = 0.0;

            //Adding contribution from slave to Kmm
            for (IndexType local_master_index = initial_sys_size; local_master_index < rEquationIds.size(); local_master_index++) {
                for (IndexType local_master_index_other = initial_sys_size; local_master_index_other < rEquationIds.size(); local_master_index_other++) {


                    rLHSContribution(local_master_index, local_master_index_other) += weights_corresponding_to_master[local_master_index-initial_sys_size] *
                                                                                        rLHSContribution(slaves_corresponding_to_masters[local_master_index], slaves_corresponding_to_masters[local_master_index_other]) * weights_corresponding_to_master[local_master_index-initial_sys_size];

                }
            }

            slave_count ++ ;
        } // Loop over all the slaves for this node

        // For K(u,s) and K(s,u)
        for (auto local_slave_index : local_slave_index_vector) { // Loop over all the slaves for this node
            for (auto local_internal_index : local_internal_index_vector) { // Loop over all the local equation ids
                rLHSContribution(local_slave_index, local_internal_index) = 0.0;
                rLHSContribution(local_internal_index, local_slave_index) = 0.0;
            }
        }

        KRATOS_CATCH("Applying Multipoint constraints failed ..");

    }

    /**
     * @brief This method reconstructs the slave solution for iteration step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void ReconstructSlaveSolutionAfterSolve(
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb
        )
    {

        const int number_of_constraints = static_cast<int>(mGlobalMasterSlaveRelations.size());
        // Getting the beginning iterator

        MasterSlaveRelationContainerType::iterator constraints_begin = mGlobalMasterSlaveRelations.begin();
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        //contributions to the system
        LocalSystemMatrixType relation_matrix = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType constraint_vector = LocalSystemVectorType(0);
        EquationIdVectorType  slave_equation_ids = EquationIdVectorType(0);
        EquationIdVectorType  master_equation_ids = EquationIdVectorType(0);

        for (int i_constraints = 0; i_constraints < number_of_constraints; i_constraints++)
        {
            MasterSlaveRelationContainerType::iterator it = constraints_begin + i_constraints;

            double slave_dx_value = 0.0;
            //get the equation Ids of the constraint
            (*it).EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
            //calculate constraint's T and b matrices
            (*it).CalculateLocalSystem(relation_matrix,constraint_vector,r_current_process_info);
            int master_index = 0;
            for(auto& master_equation_id : master_equation_ids)
            {
                slave_dx_value = slave_dx_value + TSparseSpace::GetValue(rDx, master_equation_id) * relation_matrix(0,master_index);
                master_index++;
            }

            rDx[slave_equation_ids[0]] = slave_dx_value;
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

}; /* Class ResidualBasedBlockBuilderAndSolverWithConstraints */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
