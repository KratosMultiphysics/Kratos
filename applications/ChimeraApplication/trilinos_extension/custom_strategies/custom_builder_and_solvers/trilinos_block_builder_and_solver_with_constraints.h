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
#if !defined(KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS)
#define KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS

/* System includes */
#include <unordered_set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "utilities/timer.h"

/* Trilinos includes */
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#define START_TIMER(label, rank) \
    if (BaseType::mrComm.MyPID() == rank)  \
        Timer::Start(label);
#define STOP_TIMER(label, rank) \
    if (BaseType::mrComm.MyPID() == rank) \
        Timer::Stop(label);

namespace Kratos {

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
 * @class TrilinosBlockBuilderAndSolverWithConstraints
 * @ingroup TrilinosApplication
 * @brief Current class provides an extension to the trilinos b&s with constraints
 * @details
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosBlockBuilderAndSolverWithConstraints
    : public TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosBlockBuilderAndSolverWithConstraints);

    /// Definition of the base class
    typedef TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The size_t types
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
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /// Epetra definitions
    typedef Epetra_MpiComm EpetraCommunicatorType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef DofType::Pointer DofPointerType;

    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    TrilinosBlockBuilderAndSolverWithConstraints(EpetraCommunicatorType& rComm,
                                  int GuessRowSize,
                                  typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(rComm, GuessRowSize, pNewLinearSystemSolver)
    {
    }

    /**
     * @brief Default destructor.
     */
    ~TrilinosBlockBuilderAndSolverWithConstraints() override = default;

    /**
     * Copy constructor
     */
    TrilinosBlockBuilderAndSolverWithConstraints(const TrilinosBlockBuilderAndSolverWithConstraints& rOther) = delete;

    /**
     * Assignment operator
     */
    TrilinosBlockBuilderAndSolverWithConstraints& operator=(const TrilinosBlockBuilderAndSolverWithConstraints& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to perform the build the system matrix and the residual
     * vector
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void Build(typename TSchemeType::Pointer pScheme,
               ModelPart& rModelPart,
               TSystemMatrixType& rA,
               TSystemVectorType& rb) override
    {
        KRATOS_TRY
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;
        BaseType::Build(pScheme, rModelPart, rA, rb);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation choosen the size of the matrix
     * could be equal to the total number of Dofs or to the number of
     * unrestrained dofs
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     */
    void BuildLHS(typename TSchemeType::Pointer pScheme,
                  ModelPart& rModelPart,
                  TSystemMatrixType& rA) override
    {
        KRATOS_TRY
        BaseType::BuildLHS(pScheme, rModelPart, rA);
        KRATOS_CATCH("")
    }

    /**
     * @brief Build a rectangular matrix of size n*N where "n" is the number of
     * unrestrained degrees of freedom and "N" is the total number of degrees of
     * freedom involved.
     * @details This matrix is obtained by building the total matrix without the
     * lines corresponding to the fixed degrees of freedom (but keeping the
     * columns!!)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS_CompleteOnFreeRows(typename TSchemeType::Pointer pScheme,
                                     ModelPart& rModelPart,
                                     TSystemMatrixType& A) override
    {
        KRATOS_ERROR << "Method BuildLHS_CompleteOnFreeRows not implemented in "
                        "Trilinos Builder And Solver"
                     << std::endl;
    }


    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is
     * possible to solve just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildAndSolve(typename TSchemeType::Pointer pScheme,
                       ModelPart& rModelPart,
                       TSystemMatrixType& rA,
                       TSystemVectorType& rDx,
                       TSystemVectorType& rb) override
    {
        KRATOS_TRY

        if (BaseType::GetEchoLevel() > 0)
            START_TIMER("Build", 0)

        Build(pScheme, rModelPart, rA, rb);

        if (BaseType::GetEchoLevel() > 0)
            STOP_TIMER("Build", 0)

        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        if(global_num_constraints > 0) {
            if (BaseType::GetEchoLevel() > 0) START_TIMER("ApplyConstraints", 0);
            ApplyConstraints(pScheme, rModelPart, rA, rb);
            if (BaseType::GetEchoLevel() > 0) STOP_TIMER("ApplyConstraints", 0);
        }

        // apply dirichlet conditions
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nBefore the solution of the system"
            << "\nSystem Matrix = " << rA << "\nunknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        if (BaseType::GetEchoLevel() > 0)
            START_TIMER("System solve time ", 0)

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        if (BaseType::GetEchoLevel() > 0)
            STOP_TIMER("System solve time ", 0)

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
            << "\nAfter the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;
        KRATOS_CATCH("")
    }

    void SystemSolveWithPhysics(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        ModelPart& rModelPart
    )
    {
        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        if(global_num_constraints > 0) {
            TSystemVectorType dx_mod(rb.Map());

            InternalSystemSolveWithPhysics(rA, dx_mod, rb, rModelPart);
            //recover solution of the original problem
            // TODO: Sparse matrix vector multiplication to get Dx
            TSparseSpace::Mult(*mpT, dx_mod, rDx);
        } else {
            InternalSystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        }
    }

    /**
     * @brief This is a call to the linear system solver
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void InternalSystemSolveWithPhysics(TSystemMatrixType &rA,
                                TSystemVectorType &rDx,
                                TSystemVectorType &rb,
                                ModelPart &rModelPart)
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(rb) != 0)
            norm_b = TSparseSpace::TwoNorm(rb);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            if (BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded())
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(
                    rA, rDx, rb, BaseType::mDofSet, rModelPart);

            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        }
        else
        {
            TSparseSpace::SetToZero(rDx);
            KRATOS_WARNING(
                "TrilinosResidualBasedBlockBuilderAndSolver")
                << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // prints informations about the current time
        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", (BaseType::GetEchoLevel() > 1))
            << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the
     * number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(typename TSchemeType::Pointer pScheme,
                  ModelPart& rModelPart,
                  TSystemVectorType& rb) override
    {
        KRATOS_TRY
        BaseType::BuildRHS(pScheme, rModelPart, rb);

        KRATOS_CATCH("")
    }

    /**
     * @brief Resizes the system matrix and the vector according to the number
     * of dos in the current rModelPart. This function also decides on the
     * sparsity pattern and the graph of the trilinos csr matrix
     * @param pScheme The integration scheme considered
     * @param rpA The LHS matrix
     * @param rpDx The Unknowns vector
     * @param rpd The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void ResizeAndInitializeVectors(typename TSchemeType::Pointer pScheme,
                                    TSystemMatrixPointerType& rpA,
                                    TSystemVectorPointerType& rpDx,
                                    TSystemVectorPointerType& rpb,
                                    ModelPart& rModelPart) override
    {
        KRATOS_TRY
        BaseType::ResizeAndInitializeVectors(pScheme, rpA, rpDx, rpb, rModelPart);

        ConstructMasterSlaveConstraintsStructure(rModelPart);
        KRATOS_CATCH("")
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

protected:
    ///@name Protected static Member Variables
    ///@{

    TSystemMatrixPointerType mpT;
    TSystemMatrixPointerType mpL;
    TSystemVectorPointerType mpConstantVector; /// This is vector containing the rigid movement of the constraint
    TSystemVectorPointerType mpModDx; /// This is vector containing the rigid movement of the constraint
    std::vector<IndexType> mSlaveIds;  /// The equation ids of the slaves
    std::vector<IndexType> mMasterIds; /// The equation ids of the master
    std::unordered_set<IndexType> mInactiveSlaveDofs; /// The set containing the inactive slave dofs

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb) override
    {
        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        if (global_num_constraints > 0) {
            // We make T and L matrices and C vector
            BuildMasterSlaveConstraints(rModelPart);
            // TSparseSpace::WriteMatrixMarketMatrix("T_matrix.mm", *mpT, false);
            // TSparseSpace::WriteMatrixMarketMatrix("A_matrix.mm", rA, false);
            // TSparseSpace::WriteMatrixMarketVector("b_vec.mm", rb);
            KRATOS_INFO("1 Two norm b : ")<<TSparseSpace::TwoNorm(rb)<<std::endl;
            {
                TSystemVectorType res_b(rb.Map());
                TSparseSpace::TransposeMult(*mpT, rb, res_b);
                // TSparseSpace::WriteMatrixMarketVector("res_b_vec.mm", res_b);
                res_b.GlobalAssemble();
                KRATOS_INFO("2 Two norm res : ")<<TSparseSpace::TwoNorm(res_b)<<std::endl;
                TSparseSpace::Copy(res_b, rb);
            }
            KRATOS_INFO("3. Two norm b : ")<<TSparseSpace::TwoNorm(rb)<<std::endl;;
            // First we do aux = T'A
            TSystemMatrixType aux_mat(Copy, mpT->RowMap(), 50);
            EpetraExt::MatrixMatrix::Multiply(*mpT, true, rA, false, aux_mat, true);
            // Now we do A = aux*T
            EpetraExt::MatrixMatrix::Multiply(aux_mat, false, *mpT, false, rA, true);

            double inf_norm = rA.NormInf(); // Direct call from epetra
            KRATOS_INFO("4. Inf norm rA : ")<<inf_norm<<std::endl;
            KRATOS_INFO("5. One norm  T : ")<<mpT->NormOne()<<std::endl;
            std::vector<int> active_slave_eq_ids;
            active_slave_eq_ids.reserve(BaseType::mLocalSystemSize);

            // Apply diagonal values on slaves
            for (int i = 0; i < static_cast<int>(mSlaveIds.size()); ++i) {
                const IndexType slave_equation_id = mSlaveIds[i];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    active_slave_eq_ids.push_back(slave_equation_id);
                }
            }
            double zero = 0.0;
            for(const int& inactive_eq_id : active_slave_eq_ids){
                rA.ReplaceGlobalValues(inactive_eq_id, 1, &inf_norm, &inactive_eq_id);
                rb.ReplaceGlobalValues(1, &inactive_eq_id, &zero);
            }
        }
    }

    virtual void BuildMasterSlaveConstraints(ModelPart& rModelPart)
    {
        KRATOS_TRY

        TSparseSpace::SetToZero(*mpT);
        // TSparseSpace::SetToZero(mConstantVector);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Vector containing the localization in the system of the different terms
        DofsVectorType slave_dof_list, master_dof_list;

        // Contributions to the system
        Matrix transformation_matrix = LocalSystemMatrixType(0, 0);
        Vector constant_vector = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        EquationIdVectorType slave_equation_ids, master_equation_ids;

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        // We clear the set
        mInactiveSlaveDofs.clear();

        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = it_const->IsDefined(ACTIVE) ? constraint_is_active = it_const->Is(ACTIVE) : true;

            if (constraint_is_active) {
                it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);

                // Assemble transformation matrix
                AssembleTMatrixContribution(*mpT, transformation_matrix, slave_equation_ids, master_equation_ids);
                // Assemble the constant vector
                TSparseSpace::AssembleRHS(*mpConstantVector, constant_vector, slave_equation_ids);
            } else { // Taking into account inactive constraints
                it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
                mInactiveSlaveDofs.insert(slave_equation_ids.begin(), slave_equation_ids.end());
            }
        }


        // Setting the master dofs into the T and C system
        {
            // MasterIds
            {
                std::vector<int> master_dof_ids(mMasterIds.begin(), mMasterIds.end());
                double value = 1.0;
                for(const int& master_eq_id : master_dof_ids){
                    mpT->ReplaceGlobalValues(master_eq_id, 1, &value, &master_eq_id);
                }
            }
            // InactiveSlaveDofs
            {
                std::vector<int> inactive_dof_ids(mInactiveSlaveDofs.begin(), mInactiveSlaveDofs.end());
                double value = 1.0;
                for(const int& inactive_eq_id : inactive_dof_ids){
                    mpT->ReplaceGlobalValues(inactive_eq_id, 1, &value, &inactive_eq_id);
                }

            }
        }

        mpT->GlobalAssemble();
        // mpL->GlobalAssemble();

        KRATOS_CATCH("")
    }

    virtual void ConstructMasterSlaveConstraintsStructure(ModelPart& rModelPart)
    {
        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        if (global_num_constraints > 0) {
            IndexType number_of_local_dofs = BaseType::mLastMyId - BaseType::mFirstMyId;
            std::vector<int> local_eq_ids;
            local_eq_ids.reserve(number_of_local_dofs);
            // generate map - use the "temp" array here
            for (IndexType i = 0; i != number_of_local_dofs; i++)
                local_eq_ids.push_back( BaseType::mFirstMyId + i );

            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // Vector containing the localization in the system of the different terms
            DofsVectorType slave_dof_list, master_dof_list;

            // Constraint initial iterator
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            std::map<IndexType, std::unordered_set<IndexType>> indices;

            Element::EquationIdVectorType slave_ids(3);
            Element::EquationIdVectorType master_ids(3);

            for (int i_const = 0; i_const < static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i_const) {
                auto it_const = it_const_begin + i_const;

                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = true;
                if( it_const->IsDefined(ACTIVE) ) {
                    constraint_is_active = it_const->Is(ACTIVE);
                }

                if(constraint_is_active) {
                    it_const->EquationIdVector(slave_ids, master_ids, r_current_process_info);

                    // Slave DoFs
                    for (auto &id_i : slave_ids) {
                        indices[id_i].insert(master_ids.begin(), master_ids.end());
                    }
                }
            }
            mSlaveIds.clear();
            mMasterIds.clear();
            for(const auto& slave_masters_pair : indices){
                mSlaveIds.push_back(slave_masters_pair.first);
            }
            for(const auto& local_eq_id:local_eq_ids)
                indices[local_eq_id].insert(local_eq_id);

            // Count the row sizes
            std::size_t nnz = 10; // This is a guess. This means each slave has some number of masters
            const Epetra_Map row_map(-1, local_eq_ids.size(), local_eq_ids.data(), 0, BaseType::mrComm);

            std::unordered_set<IndexType> column_eq_ids_set;
            for(auto& slave_masters_pair : indices){
                for (auto it = slave_masters_pair.second.begin(); it != slave_masters_pair.second.end(); ++it) {
                    column_eq_ids_set.insert(*it);
                }
            }
            std::vector<int> col_eq_ids_vector(column_eq_ids_set.begin(), column_eq_ids_set.end());
            column_eq_ids_set.clear();
            const Epetra_Map col_map(-1, col_eq_ids_vector.size(), col_eq_ids_vector.data(), 0, BaseType::mrComm);
            Epetra_FECrsGraph t_graph(Copy, row_map, col_map, nnz);

            // Actually inserting indices into the graph
            int slave_eq_id = 0;
            for(auto& slave_masters_pair : indices){
                slave_eq_id = slave_masters_pair.first;
                std::vector<int> master_eq_ids(slave_masters_pair.second.begin(), slave_masters_pair.second.end());

                int ierr = t_graph.InsertGlobalIndices(1, &slave_eq_id, master_eq_ids.size(), master_eq_ids.data());
                KRATOS_ERROR_IF(ierr < 0)
                    << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr
                    << std::endl;

                slave_masters_pair.second.clear(); //deallocating the memory
            }

            int ierr = t_graph.GlobalAssemble();
            KRATOS_ERROR_IF(ierr < 0)
                << ": Epetra failure in Graph.GlobalAssemble. Error code: " << ierr
                << std::endl;

            // generate a new matrix pointer according to this graph
            TSystemMatrixPointerType p_new_t =
                TSystemMatrixPointerType(new TSystemMatrixType(Copy, t_graph));
            mpT.swap(p_new_t);

            TSystemMatrixPointerType p_new_l =
                TSystemMatrixPointerType(new TSystemMatrixType(Copy, t_graph));
            mpL.swap(p_new_l);

            TSystemVectorPointerType p_new_Dx =
                TSystemVectorPointerType(new TSystemVectorType(row_map));
            mpModDx.swap(p_new_Dx);

            TSystemVectorPointerType p_new_Cvec =
                TSystemVectorPointerType(new TSystemVectorType(row_map));
            mpConstantVector.swap(p_new_Cvec);

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
    int GetGlobalNumberOfConstraints(ModelPart& rModelPart){
        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int local_num_constraints = rModelPart.MasterSlaveConstraints().size();
        const int global_num_constraints = r_data_comm.SumAll(local_num_constraints);
        return global_num_constraints;
    }


    inline static void AssembleTMatrixContribution(
        TSystemMatrixType& rMatrix,
        LocalSystemMatrixType& rContribution,
        EquationIdVectorType& rRowEquationIds,
        EquationIdVectorType& rColEquationIds
    )
    {
        //size Epetra vectors
        Epetra_IntSerialDenseVector row_indices(rRowEquationIds.size());
        Epetra_IntSerialDenseVector col_indices(rColEquationIds.size());
        Epetra_SerialDenseMatrix values(rRowEquationIds.size(), rColEquationIds.size());

        //fill epetra vectors
        int loc_i = 0;
        for (unsigned int i = 0; i < rRowEquationIds.size(); ++i)
        {
            row_indices[loc_i] = rRowEquationIds[i];

            int loc_j = 0;
            for (unsigned int j = 0; j < rColEquationIds.size(); ++j)
            {
                col_indices[loc_j] = rColEquationIds[j];
                values(loc_i, loc_j) = rContribution(i, j);
                ++loc_j;
            }
            ++loc_i;
        }

        int ierr = rMatrix.SumIntoGlobalValues(row_indices, col_indices, values);
        if(ierr != 0) KRATOS_ERROR<<"Epetra failure found"<<std::endl;
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
}; /* Class TrilinosBlockBuilderAndSolverWithConstraints */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS  defined */
