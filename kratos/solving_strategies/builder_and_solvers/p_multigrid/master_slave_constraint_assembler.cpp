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
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // ParseDiagonalScaling, GetDiagonalScaleFactor
#include "includes/define.h" // KRATOS_TRY, KRATOS_CATCH
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility
#include "utilities/atomic_utilities.h" // AtomicAdd


namespace Kratos {


template <class TSparse, class TDense>
struct MasterSlaveConstraintAssembler<TSparse,TDense>::Impl
{
    std::vector<std::size_t> mSlaveIds;

    std::vector<std::size_t> mMasterIds;

    std::unordered_set<std::size_t> mInactiveSlaveIds;

    std::unique_ptr<Scaling> mpDiagonalScaling;

    int mVerbosity;
};


template <class TSparse, class TDense>
MasterSlaveConstraintAssembler<TSparse,TDense>::MasterSlaveConstraintAssembler() noexcept
    : MasterSlaveConstraintAssembler(Parameters())
{}


template <class TSparse, class TDense>
MasterSlaveConstraintAssembler<TSparse,TDense>::MasterSlaveConstraintAssembler(Parameters Settings)
    : MasterSlaveConstraintAssembler(Settings, "unnamed")
{}


template <class TSparse, class TDense>
MasterSlaveConstraintAssembler<TSparse,TDense>::MasterSlaveConstraintAssembler(Parameters Settings,
                                                                               std::string&& rInstanceName)
    : Base(ConstraintImposition::MasterSlave, std::move(rInstanceName)),
      mpImpl(new Impl)
{
    KRATOS_TRY
    // Parse diagonal scaling and validate other settings.
    Parameters default_parameters = this->GetDefaultParameters();
    Parameters default_diagonal_scaling = default_parameters["diagonal_scaling"].Clone();
    default_parameters.RemoveValue("diagonal_scaling");
    std::optional<Parameters> maybe_diagonal_scaling = Settings.Has("diagonal_scaling") ? Settings["diagonal_scaling"].Clone() : std::optional<Parameters>();
    if (maybe_diagonal_scaling.has_value()) Settings.RemoveValue("diagonal_scaling");
    Settings.ValidateAndAssignDefaults(default_parameters);
    Settings.AddValue("diagonal_scaling", maybe_diagonal_scaling.has_value() ? maybe_diagonal_scaling.value() : default_diagonal_scaling);
    mpImpl->mpDiagonalScaling = std::make_unique<Scaling>(Settings["diagonal_scaling"]);
    mpImpl->mVerbosity = Settings["verbosity"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
MasterSlaveConstraintAssembler<TSparse,TDense>::~MasterSlaveConstraintAssembler() = default;


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

    mpImpl->mSlaveIds.clear();
    mpImpl->mMasterIds.clear();
    for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        if (indices[i].size() == 0) // Master dof!
            mpImpl->mMasterIds.push_back(i);
        else // Slave dof
            mpImpl->mSlaveIds.push_back(i);
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
    mpImpl->mInactiveSlaveIds.clear();
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
            mpImpl->mInactiveSlaveIds.insert(tls_inactive_slave_dofs.begin(), tls_inactive_slave_dofs.end());
        }
    }

    // Setting the master dofs into the T and C system
    for (auto eq_id : mpImpl->mMasterIds) {
        std::scoped_lock<LockObject> lock(mutexes[eq_id]);
        this->GetConstraintGapVector()[eq_id] = 0.0;
        this->GetRelationMatrix()(eq_id, eq_id) = 1.0;
    }

    // Setting inactive slave dofs in the T and C system
    for (auto eq_id : mpImpl->mInactiveSlaveIds) {
        std::scoped_lock<LockObject> lock(mutexes[eq_id]);
        this->GetConstraintGapVector()[eq_id] = 0.0;
        this->GetRelationMatrix()(eq_id, eq_id) = 1.0;
    }

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Initialize(typename TSparse::MatrixType& rLhs,
                                                                typename TSparse::VectorType& rRhs,
                                                                [[maybe_unused]] typename Base::DofSet::iterator itDofBegin,
                                                                [[maybe_unused]] typename Base::DofSet::iterator itDofEnd)
{
    KRATOS_TRY

    // Output the relation matrix and constraint gaps if the verbosity is high enough.
    if (4 <= mpImpl->mVerbosity) {
        KRATOS_INFO(this->GetName() + ": ")
            << "write relation matrix to "
            << (std::filesystem::current_path() / "relation_matrix.mm")
            << "\n";
        TSparse::WriteMatrixMarketMatrix("relation_matrix.mm", this->GetRelationMatrix(), false);

    KRATOS_INFO(this->GetName() + ": ")
        << "write constraint gaps to "
        << (std::filesystem::current_path() / "constraint_gaps.mm")
        << "\n";
    TSparse::WriteMatrixMarketVector("constraint_gaps.mm", this->GetConstraintGapVector());
    } // if 4 <= verbosity

    // Compute the transposed matrix of the global relation matrix
    {
        // Storage for an intermediate matrix transpose(relation_matrix) * lhs_matrix
        typename TSparse::MatrixType left_multiplied_lhs(this->GetRelationMatrix().size2(), rLhs.size2());
        {
            // Transpose the relation matrix
            typename TSparse::MatrixType transposed_relation_matrix(this->GetRelationMatrix().size2(), this->GetRelationMatrix().size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix(transposed_relation_matrix, this->GetRelationMatrix(), 1.0);

            typename TSparse::VectorType b_modified(rRhs.size());
            BalancedProduct<TSparse,TSparse,TSparse>(transposed_relation_matrix, rRhs, b_modified);
            TSparse::Copy(b_modified, rRhs);

            SparseMatrixMultiplicationUtility::MatrixMultiplication(transposed_relation_matrix, rLhs, left_multiplied_lhs);
        } // deallocate transposed_relation_matrix

        SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, this->GetRelationMatrix(), rLhs);
    } // deallocate left_multiplied_lhs

    // Compute the scale factor for slave DoFs.
    mpImpl->mpDiagonalScaling->template Cache<TSparse>(rLhs);
    typename TSparse::DataType diagonal_scale_factor = mpImpl->mpDiagonalScaling->Evaluate();

    // Apply diagonal values on slaves.
    block_for_each(mpImpl->mSlaveIds, [this, &rLhs, &rRhs, diagonal_scale_factor](const auto iSlave){
        if (mpImpl->mInactiveSlaveIds.find(iSlave) == mpImpl->mInactiveSlaveIds.end()) {
            rLhs(iSlave, iSlave) = diagonal_scale_factor;
            rRhs[iSlave] = 0.0;
        }
    });

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
bool
MasterSlaveConstraintAssembler<TSparse,TDense>::FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                                                     typename TSparse::VectorType& rSolution,
                                                                     typename TSparse::VectorType& rRhs,
                                                                     PMGStatusStream::Report& rReport)
{
    rReport.maybe_constraint_residual = 0;
    rReport.constraints_converged = true;
    return true;
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Apply(typename TSparse::VectorType& rSolution) const
{
    KRATOS_TRY
    typename TSparse::VectorType constrained_solution(rSolution.size());
    TSparse::SetToZero(constrained_solution);
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetRelationMatrix(),
                                             rSolution,
                                             constrained_solution);
    rSolution.swap(constrained_solution);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Clear()
{
    Base::Clear();
    mpImpl->mSlaveIds = decltype(mpImpl->mSlaveIds)();
    mpImpl->mMasterIds = decltype(mpImpl->mMasterIds)();
    mpImpl->mInactiveSlaveIds = decltype(mpImpl->mInactiveSlaveIds)();
}


template <class TSparse, class TDense>
Parameters MasterSlaveConstraintAssembler<TSparse,TDense>::GetDefaultParameters()
{
    return Parameters(R"({
"method" : "master_slave",
"diagonal_scaling" : "norm",
"verbosity" : 1
})");
}


template class MasterSlaveConstraintAssembler<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;

template class MasterSlaveConstraintAssembler<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos
