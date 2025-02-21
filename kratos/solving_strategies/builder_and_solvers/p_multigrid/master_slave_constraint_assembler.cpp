//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.hpp" // MasterSlaveConstraintAssembler
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // GetDiagonalScaleFactor
#include "solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp" // MakeSparseTopology
#include "includes/define.h" // KRATOS_TRY, KRATOS_CATCH
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility
#include "utilities/atomic_utilities.h" // AtomicAdd


namespace Kratos {


template <class TSparse, class TDense>
MasterSlaveConstraintAssembler<TSparse,TDense>::MasterSlaveConstraintAssembler() noexcept
    : MasterSlaveConstraintAssembler(Parameters())
{}


template <class TSparse, class TDense>
MasterSlaveConstraintAssembler<TSparse,TDense>::MasterSlaveConstraintAssembler(Parameters Settings)
    : Base(ConstraintImposition::MasterSlave),
      mSlaveIds(),
      mMasterIds(),
      mInactiveSlaveIds(),
      mDiagonalScaling(DiagonalScaling::None)
{
    KRATOS_TRY
    Settings.ValidateAndAssignDefaults(MasterSlaveConstraintAssembler::GetDefaultParameters());
    mDiagonalScaling = ParseDiagonalScaling(Settings);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Allocate(const typename Base::ConstraintArray& rConstraints,
                                                              const ProcessInfo& rProcessInfo,
                                                              typename TSparse::MatrixType& rLhs,
                                                              typename TSparse::VectorType& rSolution,
                                                              typename TSparse::VectorType& rRhs,
                                                              typename Base::DofSet& rDofSet)
{
    KRATOS_TRY
    std::vector<std::unordered_set<IndexType>> indices(rDofSet.size());
    std::vector<LockObject> mutexes(indices.size());
    const auto it_const_begin = rConstraints.begin();

    #pragma omp parallel
    {
        Element::EquationIdVectorType slave_ids;
        Element::EquationIdVectorType master_ids;
        std::unordered_map<IndexType, std::unordered_set<IndexType>> temp_indices;

        #pragma omp for nowait
        for (int i_const = 0; i_const < static_cast<int>(rConstraints.size()); ++i_const) {
            auto it_const = it_const_begin + i_const;
            it_const->EquationIdVector(slave_ids, master_ids, rProcessInfo);

            // Slave DoFs
            for (auto &id_i : slave_ids) {
                temp_indices[id_i].insert(master_ids.begin(), master_ids.end());
            }
        }

        // Merging temporary indices into global rows.
        for (auto& pair_temp_indices : temp_indices) {
            std::scoped_lock<LockObject> lock(mutexes[pair_temp_indices.first]);
            indices[pair_temp_indices.first].insert(pair_temp_indices.second.begin(), pair_temp_indices.second.end());
        }
    }

    mSlaveIds.clear();
    mMasterIds.clear();
    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        if (indices[i].size() == 0) // Master dof!
            mMasterIds.push_back(i);
        else // Slave dof
            mSlaveIds.push_back(i);
        indices[i].insert(i); // Ensure that the diagonal is there in T
    }

    MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                         indices.size(),
                                                         this->GetRelationMatrix(),
                                                         /*EnsureDiagonal=*/false);
    this->GetConstraintGapVector().resize(indices.size(), false);

    KRATOS_CATCH("")

    KRATOS_TRY
    MergeMatrices<typename TSparse::DataType>(rLhs, this->GetRelationMatrix());
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Assemble(const typename Base::ConstraintArray& rConstraints,
                                                              const ProcessInfo& rProcessInfo,
                                                              typename Base::DofSet& rDofSet,
                                                              const bool AssembleLhs,
                                                              const bool AssembleRhs)
{
    KRATOS_TRY

    /// @todo @p AssembleLhs and @p AssembleRhs are currently ignored,
    ///       and everything gets assembled. Think it through what needs
    ///       to be done and what can be omitted for each case. - matekelemen

    // Init.
    mInactiveSlaveIds.clear();
    TSparse::SetToZero(this->GetRelationMatrix());
    TSparse::SetToZero(this->GetConstraintGapVector());
    std::vector<LockObject> mutexes(rDofSet.size());

    // Declare local containers.
    typename TDense::MatrixType local_relation_matrix;
    typename TDense::VectorType local_constraint_gap_vector;
    Element::EquationIdVectorType slave_equation_ids, master_equation_ids;
    const int number_of_constraints = static_cast<int>(rConstraints.size());

    #pragma omp parallel firstprivate(local_relation_matrix, local_constraint_gap_vector, slave_equation_ids, master_equation_ids)
    {
        std::unordered_set<IndexType> tls_inactive_slave_dofs;

        #pragma omp for
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rConstraints.begin() + i_const;

            // If the constraint is active
            if (it_const->IsActive()) {
                it_const->EquationIdVector(slave_equation_ids,
                                           master_equation_ids,
                                           rProcessInfo);
                it_const->CalculateLocalSystem(local_relation_matrix,
                                               local_constraint_gap_vector,
                                               rProcessInfo);

                for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                    const IndexType i_global = slave_equation_ids[i];

                    // Assemble matrix row.
                    {
                        std::scoped_lock<LockObject> lock(mutexes[i_global]);
                        MapRowContribution<TSparse,TDense>(this->GetRelationMatrix(),
                                                           local_relation_matrix,
                                                           i_global,
                                                           i,
                                                           master_equation_ids);
                    }

                    // Assemble constant vector
                    KRATOS_ERROR_IF_NOT(i_global < this->GetConstraintGapVector().size())
                        << "constraint gap vector size " << this->GetConstraintGapVector().size()
                        << " <= insertion index " << i_global << "\n";
                    const double constant_value = local_constraint_gap_vector[i];
                    AtomicAdd(this->GetConstraintGapVector()[i_global], typename TSparse::DataType(constant_value));
                }
            } else { // Taking into account inactive constraints
                tls_inactive_slave_dofs.insert(slave_equation_ids.begin(), slave_equation_ids.end());
            }
        }

        // We merge all the sets in one thread
        #pragma omp critical
        {
            mInactiveSlaveIds.insert(tls_inactive_slave_dofs.begin(), tls_inactive_slave_dofs.end());
        }
    }

    // Setting the master dofs into the T and C system
    for (auto eq_id : mMasterIds) {
        std::scoped_lock<LockObject> lock(mutexes[eq_id]);
        this->GetConstraintGapVector()[eq_id] = 0.0;
        this->GetRelationMatrix()(eq_id, eq_id) = 1.0;
    }

    // Setting inactive slave dofs in the T and C system
    for (auto eq_id : mInactiveSlaveIds) {
        std::scoped_lock<LockObject> lock(mutexes[eq_id]);
        this->GetConstraintGapVector()[eq_id] = 0.0;
        this->GetRelationMatrix()(eq_id, eq_id) = 1.0;
    }

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Initialize(typename TSparse::MatrixType& rLhs,
                                                                typename TSparse::VectorType& rRhs)
{
    KRATOS_TRY
    // Compute the transposed matrix of the global relation matrix
    {
        // Storage for an intermediate matrix transpose(relation_matrix) * lhs_matrix
        typename TSparse::MatrixType left_multiplied_lhs(this->GetRelationMatrix().size2(), rLhs.size2());
        {
            // Transpose the relation matrix
            typename TSparse::MatrixType transposed_relation_matrix(this->GetRelationMatrix().size2(), this->GetRelationMatrix().size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix(transposed_relation_matrix, this->GetRelationMatrix(), 1.0);

            typename TSparse::VectorType b_modified(rRhs.size());
            TSparse::Mult(transposed_relation_matrix, rRhs, b_modified);
            TSparse::Copy(b_modified, rRhs);

            SparseMatrixMultiplicationUtility::MatrixMultiplication(transposed_relation_matrix, rLhs, left_multiplied_lhs);
        } // deallocate transposed_relation_matrix

        SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, this->GetRelationMatrix(), rLhs);
    } // deallocate left_multiplied_lhs

    // Compute the scale factor for slave DoFs.
    typename TSparse::DataType diagonal_scale_factor = GetDiagonalScaleFactor<TSparse>(rLhs, this->mDiagonalScaling);

    // Apply diagonal values on slaves.
    block_for_each(this->mSlaveIds, [this, &rLhs, &rRhs, diagonal_scale_factor](const auto iSlave){
        if (this->mInactiveSlaveIds.find(iSlave) == this->mInactiveSlaveIds.end()) {
            rLhs(iSlave, iSlave) = diagonal_scale_factor;
            rRhs[iSlave] = 0.0;
        }
    });
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
typename MasterSlaveConstraintAssembler<TSparse,TDense>::Base::Status
MasterSlaveConstraintAssembler<TSparse,TDense>::FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                                                     typename TSparse::VectorType& rSolution,
                                                                     typename TSparse::VectorType& rRhs,
                                                                     const std::size_t iIteration)
{
    return typename Base::Status {/*finished=*/true, /*converged=*/true};
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Finalize(typename TSparse::MatrixType& rLhs,
                                                              typename TSparse::VectorType& rSolution,
                                                              typename TSparse::VectorType& rRhs,
                                                              typename Base::DofSet& rDofSet)
{
    KRATOS_TRY
    typename TSparse::VectorType original_solution(rSolution.size());
    TSparse::SetToZero(original_solution);
    TSparse::Mult(this->GetRelationMatrix(),
                  rSolution,
                  original_solution);
    rSolution.swap(original_solution);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Clear()
{
    Base::Clear();
    mSlaveIds = decltype(mSlaveIds)();
    mMasterIds = decltype(mMasterIds)();
    mInactiveSlaveIds = decltype(mInactiveSlaveIds)();
}


template <class TSparse, class TDense>
Parameters MasterSlaveConstraintAssembler<TSparse,TDense>::GetDefaultParameters()
{
    return Parameters(R"({
"method" : "master_slave",
"diagonal_scaling" : "none"
})");
}


template class MasterSlaveConstraintAssembler<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;

template class MasterSlaveConstraintAssembler<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos
