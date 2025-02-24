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
#include "solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.hpp" // AugmentedLagrangeConstraintAssembler
#include "solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp" // MapRowContribution
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility

// System includes
#include <unordered_set> // std::unordered_set


namespace Kratos {


namespace detail {


/// @brief Checks whether two sorted sets have at least one item in common.
/// @tparam TIndex Set value type.
/// @param itLeftBegin Pointer to the beginning of the first set.
/// @param itLeftEnd Sentinel of the first set.
/// @param itRightBegin Pointer to the beginning of the second set.
/// @param itRightEnd Sentinel of the second set.
/// @return True if the two sets have a non-empty intersection, false otherwise.
template <class TIndex>
bool HasIntersection(const TIndex* itLeftBegin,
                     const TIndex* itLeftEnd,
                     const TIndex* itRightBegin,
                     const TIndex* itRightEnd) noexcept
{
    // Early exit on empty ranges.
    if (itLeftBegin == itLeftEnd or itRightBegin == itRightEnd)
        return false;

    /// @todo Shrink the search ranges to only process overlapping index spaces.

    auto it_left = itLeftBegin;
    auto it_right = itRightBegin;
    while (it_left < itLeftEnd) {
        if (*it_left == *it_right) return true;
        it_right = std::lower_bound(it_right, itRightEnd, *it_left);
        if (it_right == itRightEnd) return false;
        else if (*it_right == *it_left) return true;
        else it_left = std::lower_bound(it_left, itLeftEnd, *it_right);
    } // for it_left in range(itLeftBegin, itLeftEnd)

    return false;
}


} // namespace detail


template <class TSparse, class TDense>
struct AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Impl
{
    using Interface = AugmentedLagrangeConstraintAssembler<TSparse,TDense>;

    /// @brief A map associating slave IDs with constraint indices.
    std::unordered_map<std::size_t,std::size_t> mConstraintIdToIndexMap;

    std::optional<typename TSparse::MatrixType> mMaybeTransposeRelationMatrix;

    /// @details Some models may indirectly set Dirichlet conditions on DoFs through
    ///          multifreedom constraints. For example, if a constraint involves two
    ///          DoFs and one of them has a Dirichlet condition set, the other one
    ///          should theoretically also inherit a scaled version of the condition.
    ///          This is unfortunately not exactly satisfied with penalty or augmented
    ///          lagrange imposition, and trying to force it via penalty factors may
    ///          unnecessarily render the system ill conditioned.
    ///          To handle such cases, the dirichlet conditions are forced to propagate
    ///          through constraints, and this array stores the indices of DoFs that were
    ///          force
    std::vector<typename Dof<typename TDense::DataType>::IndexType> mForcedDirichletDoFs;

    int mVerbosity;

    void PropagateDirichletConditions(typename Interface::DofSet::iterator itDofBegin,
                                      typename Interface::DofSet::iterator itDofEnd,
                                      const typename TSparse::MatrixType& rRelationMatrix)
    {
        KRATOS_TRY

        /// @todo Make this more robust (@matekelemen).

        mForcedDirichletDoFs.clear();
        std::unordered_set<std::size_t> forced_dirichlet_dofs;
        LockObject mutex;

        IndexPartition<std::size_t>(rRelationMatrix.size1()).for_each([&forced_dirichlet_dofs, &mutex, &rRelationMatrix, itDofBegin](const std::size_t i_constraint){
            const auto i_entry_begin = rRelationMatrix.index1_data()[i_constraint];
            const auto i_entry_end = rRelationMatrix.index1_data()[i_constraint +1];

            for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_dof = rRelationMatrix.index2_data()[i_entry];
                const Dof<typename TDense::DataType>& r_dof = *(itDofBegin + i_dof);
                KRATOS_ERROR_IF_NOT(i_dof == r_dof.EquationId());

                if (r_dof.IsFixed()) {
                    KRATOS_ERROR_IF_NOT(i_entry_end - i_entry_begin == 2)
                        << "A Dirichlet condition is set on DoF " << i_dof << " belonging to node " << r_dof.Id() << ", "
                        << "and it also appears in a multifreedom constraint with " << i_entry_end - i_entry_begin - (i_entry_begin == i_entry_end ? 0 : 1) << " other DoFs. "
                        << "In such cases, multifreedom constraints with only 2 DoFs are supported due to performance considerations.";

                    const auto i_forced_dof = i_dof == rRelationMatrix.index2_data()[i_entry_begin]
                        ? rRelationMatrix.index2_data()[i_entry_begin + 1]
                        : rRelationMatrix.index2_data()[i_entry_begin];

                    std::scoped_lock<LockObject> lock(mutex);
                    forced_dirichlet_dofs.insert(i_forced_dof);
                } // if r_dof.IsFixed()
            } // for i_entry in range(i_entry_begin, i_entry_end)
        }); // for i_constraint in range(rRelationMatrix.size1())

        mForcedDirichletDoFs.reserve(forced_dirichlet_dofs.size());
        std::copy(forced_dirichlet_dofs.begin(),
                  forced_dirichlet_dofs.end(),
                  std::back_inserter(mForcedDirichletDoFs));

        KRATOS_CATCH("")
    }

    void ErasePropagatedDirichletConditions(typename Interface::DofSet::iterator itDofBegin,
                                            [[maybe_unused]] typename Interface::DofSet::iterator itDofEnd)
    {
        block_for_each(mForcedDirichletDoFs, [itDofBegin](const auto i_dof){
            const auto it_dof = itDofBegin + i_dof;
            it_dof->FreeDof();
        });
        mForcedDirichletDoFs.clear();
    }
}; // struct AugmentedLagrangeConstraintAssembler::Impl


template <class TSparse, class TDense>

AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler() noexcept
    : AugmentedLagrangeConstraintAssembler(Parameters())
{}


template <class TSparse, class TDense>
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler(Parameters Settings)
    : Base(ConstraintImposition::AugmentedLagrange),
      mpImpl(new Impl{/*mConstraintIdToIndexMap: */std::unordered_map<std::size_t,std::size_t>(),
                      /*mMaybeTransposeRelationMatrix: */std::optional<typename TSparse::MatrixType>(),
                      /*mForcedDirichletDoFs: */std::vector<typename Dof<typename TDense::DataType>::IndexType>(),
                      /*mVerbosity: */1})
{
    Settings.ValidateAndAssignDefaults(AugmentedLagrangeConstraintAssembler::GetDefaultParameters());

    Vector algorithmic_parameters(3);
    algorithmic_parameters[0] = Settings["penalty_factor"].Get<double>();
    algorithmic_parameters[1] = Settings["initial_lagrange_multiplier"].Get<double>();
    algorithmic_parameters[2] = Settings["tolerance"].Get<double>();
    this->SetValue(this->GetAlgorithmicParametersVariable(), algorithmic_parameters);
    this->SetValue(NL_ITERATION_NUMBER, Settings["max_iterations"].Get<int>());
    this->mpImpl->mVerbosity = Settings["verbosity"].Get<int>();
}


template <class TSparse, class TDense>
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::~AugmentedLagrangeConstraintAssembler() = default;


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Allocate(const typename Base::ConstraintArray& rConstraints,
                                                                    const ProcessInfo& rProcessInfo,
                                                                    typename TSparse::MatrixType& rLhs,
                                                                    typename TSparse::VectorType& rSolution,
                                                                    typename TSparse::VectorType& rRhs,
                                                                    typename Base::DofSet& rDofSet)
{
    KRATOS_TRY

    if (rConstraints.empty()) {
        this->GetRelationMatrix() = typename TSparse::MatrixType(0, rLhs.size2(), 0);
        this->GetRelationMatrix().index1_data()[0] = 0;
        this->GetRelationMatrix().set_filled(1, 0);
        return;
    }

    // Build constraint ID => index map.
    mpImpl->mConstraintIdToIndexMap.clear();

    {
        MasterSlaveConstraint::IndexType i_constraint = 0;
        for (const auto& r_constraint : rConstraints) {
            const auto& r_constraint_ids = r_constraint.GetData().GetValue(CONSTRAINT_LABELS);
            for (const auto constraint_id : r_constraint_ids) {
                const auto emplace_result = mpImpl->mConstraintIdToIndexMap.emplace(static_cast<std::size_t>(constraint_id),
                                                                                    i_constraint);
                if (emplace_result.second) ++i_constraint;
            } // for i_slave in slave_ids
        } // for r_constraint in rModelPart.MasterSlaveConstraints
    }

    {
        std::vector<std::unordered_set<IndexType>> indices(mpImpl->mConstraintIdToIndexMap.size());
        std::vector<LockObject> mutexes(indices.size());

        struct TLS {
            MasterSlaveConstraint::EquationIdVectorType slaves, masters;
        };

        block_for_each(rConstraints,
                       TLS(),
                       [&mutexes, &indices, &rProcessInfo, this](const auto& r_constraint, TLS& r_tls) {
            r_constraint.EquationIdVector(r_tls.slaves, r_tls.masters, rProcessInfo);
            KRATOS_ERROR_IF_NOT(r_tls.slaves.empty()) << "AugmentedLagrangeConstraintAssembler supports constraints derived from MultifreedomConstraint only";

            const auto& r_constraint_labels = r_constraint.GetData().GetValue(CONSTRAINT_LABELS);
            for (const auto i_slave : r_constraint_labels) {
                const auto i_constraint = mpImpl->mConstraintIdToIndexMap[static_cast<std::size_t>(i_slave)];
                std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                indices[i_constraint].insert(r_tls.masters.begin(), r_tls.masters.end());
            } // for i_slave in slave_ids
        }); // for r_constraint in rModelPart.MasterSlaveConstraints()

        MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                             rDofSet.size(),
                                                             this->GetRelationMatrix(),
                                                             /*EnsureDiagonal=*/false);

        this->GetConstraintGapVector().resize(mpImpl->mConstraintIdToIndexMap.size(), false);
        TSparse::SetToZero(this->GetConstraintGapVector());
    }

    {
            typename TSparse::MatrixType product;
        {
            typename TSparse::MatrixType transpose;
            SparseMatrixMultiplicationUtility::TransposeMatrix(transpose, this->GetRelationMatrix());
            SparseMatrixMultiplicationUtility::MatrixMultiplication(transpose, this->GetRelationMatrix(), product);
        }
        MergeMatrices<typename TSparse::DataType>(rLhs, product);
    }

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Assemble(const typename Base::ConstraintArray& rConstraints,
                                                                    const ProcessInfo& rProcessInfo,
                                                                    typename Base::DofSet& rDofSet,
                                                                    const bool AssembleLhs,
                                                                    const bool AssembleRhs)
{
    KRATOS_TRY

    /// @todo @p AssembleLhs and @p AssembleRhs are currently ignored,
    ///       and everything gets assembled. Think it through what needs
    ///       to be done and what can be omitted for each case. - matekelemen

    // Function-wide variables.
    std::vector<LockObject> mutexes(mpImpl->mConstraintIdToIndexMap.size());

    // Init.
    TSparse::SetToZero(this->GetRelationMatrix());
    TSparse::SetToZero(this->GetConstraintGapVector());

    struct TLS {
        MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
        std::vector<std::size_t> dof_index_array;
        std::vector<std::size_t> reordered_dof_equation_ids;
        std::vector<MasterSlaveConstraint::MatrixType::value_type> relation_matrix_row;
        MasterSlaveConstraint::MatrixType relation_matrix;
        MasterSlaveConstraint::VectorType constraint_gaps;
    };

    // Constraint assembly.
    block_for_each(rConstraints,
                   TLS(),
                   [&mutexes, &rProcessInfo, this](const MasterSlaveConstraint& r_constraint, TLS& r_tls){
        if (r_constraint.IsActive()) {
            const auto& r_constraint_ids = r_constraint.GetData().GetValue(CONSTRAINT_LABELS);
            const auto& r_dof_equation_ids = r_constraint.GetData().GetValue(CONSTRAINT_DOFS);
            r_tls.relation_matrix = r_constraint.GetData().GetValue(RELATION_MATRIX);
            const auto& r_constraint_gaps = r_constraint.GetData().GetValue(CONSTRAINT_GAPS);

            // Make sure that the constraints assembled are all MultifreedomConstraints
            // and not standard MasterSlaveConstraints, since they're unsupported by
            // constraint assemblers.
            KRATOS_ERROR_IF(r_constraint_ids.empty()) << "AugmentedLagrangeConstraintAssembler only supports MultifreedomConstraints, not MasterSlaveConstraints.";

            r_tls.relation_matrix_row.resize(r_tls.relation_matrix.size2());

            // Assemble local rows into the global relation matrix.
            for (std::size_t i_row=0ul; i_row<r_constraint_ids.size(); ++i_row) {
                const auto it_constraint_index = mpImpl->mConstraintIdToIndexMap.find(static_cast<std::size_t>(r_constraint_ids[i_row]));
                KRATOS_ERROR_IF(it_constraint_index == mpImpl->mConstraintIdToIndexMap.end())
                    << "constraint label " << i_row << " "
                    << "of constraint ID " << r_constraint.Id() << " "
                    << "is not associated with a constraint index.";

                const std::size_t i_constraint = it_constraint_index->second;
                KRATOS_ERROR_IF(this->GetRelationMatrix().size1() <= i_constraint);

                // Indirect sort the local relation matrix' row based on DoF IDs.
                // This step is required because the global relation matrix is in CSR format,
                // and Impl::MapRowContribution expects the column indices to be sorted.
                {
                    r_tls.dof_index_array.resize(r_dof_equation_ids.size());
                    std::iota(r_tls.dof_index_array.begin(), r_tls.dof_index_array.end(), 0ul);
                    std::sort(r_tls.dof_index_array.begin(),
                              r_tls.dof_index_array.end(),
                              [&r_dof_equation_ids](const std::size_t i_left, const std::size_t i_right) -> bool {
                                  return r_dof_equation_ids[i_left] < r_dof_equation_ids[i_right];
                              });

                    // Reorder the current row in the local relation matrix.
                    std::transform(r_tls.dof_index_array.begin(),
                                   r_tls.dof_index_array.end(),
                                   r_tls.relation_matrix_row.begin(),
                                   [&r_tls, i_row](const std::size_t i_column){
                                   return r_tls.relation_matrix(i_row, i_column);
                                   });
                    std::copy(r_tls.relation_matrix_row.begin(),
                              r_tls.relation_matrix_row.end(),
                              r_tls.relation_matrix.data().begin() + i_row * r_tls.relation_matrix.size2());

                    // Reorder DoF indices.
                    r_tls.reordered_dof_equation_ids.resize(r_dof_equation_ids.size());
                    std::transform(r_tls.dof_index_array.begin(),
                                   r_tls.dof_index_array.end(),
                                   r_tls.reordered_dof_equation_ids.begin(),
                                   [&r_dof_equation_ids](const auto i_column){
                                   return r_dof_equation_ids[i_column];
                                   });
                }

                {
                    std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                    MapRowContribution<TSparse,TDense>(this->GetRelationMatrix(),
                                                        r_tls.relation_matrix,
                                                        i_constraint,
                                                        i_row,
                                                        r_tls.reordered_dof_equation_ids);
                }

                AtomicAdd(this->GetConstraintGapVector()[i_constraint],
                          static_cast<typename TSparse::DataType>(r_constraint_gaps[i_row]));
            } // for i_row in range(r_tls.slave_ids.size)
        } // if r_constraint.IsActive
    }); // for r_constraint in rModelPart.MasterSlaveCosntraints

    KRATOS_CATCH("")

    KRATOS_TRY
    mpImpl->PropagateDirichletConditions(rDofSet.begin(), rDofSet.end(), this->GetRelationMatrix());
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Initialize(typename TSparse::MatrixType& rLhs,
                                                                      typename TSparse::VectorType& rRhs)
{
    KRATOS_TRY

    using SparseUtils = SparseMatrixMultiplicationUtility;
    const typename TSparse::DataType penalty_factor = this->GetPenaltyFactor();
    const typename TSparse::DataType initial_lagrange_multiplier = this->GetInitialLagrangeMultiplier();

    {
        const typename TSparse::MatrixType& r_transpose_relation_matrix = this->GetTransposeRelationMatrix();
        KRATOS_ERROR_IF(r_transpose_relation_matrix.size1() != this->GetRelationMatrix().size2()
                        || r_transpose_relation_matrix.size2() != this->GetRelationMatrix().size1()
                        || r_transpose_relation_matrix.nnz() != this->GetRelationMatrix().nnz())
            << "the transpose of the relation matrix is uninitialized or out of date "
            << "(" << this->GetRelationMatrix().size1() << "x" << this->GetRelationMatrix().size2() << "/" << this->GetRelationMatrix().nnz() << ") "
            << "(" << this->GetTransposeRelationMatrix().size1() << "x" << this->GetTransposeRelationMatrix().size2() << "/" << this->GetTransposeRelationMatrix().nnz() << ") ";

        typename TSparse::MatrixType relation_product;
        SparseUtils::MatrixMultiplication(r_transpose_relation_matrix,
                                          this->GetRelationMatrix(),
                                          relation_product);

        // Add terms to the LHS matrix.
        InPlaceMatrixAdd(rLhs,
                         relation_product,
                         penalty_factor);

        // Add terms to the RHS vector.
        typename TSparse::VectorType rhs_term(rRhs.size()),
                                     lagrange_multipliers(this->GetRelationMatrix().size1(), initial_lagrange_multiplier),
                                     constraint_gaps = this->GetConstraintGapVector();
        KRATOS_ERROR_IF_NOT(constraint_gaps.size() == lagrange_multipliers.size());

        TSparse::UnaliasedAdd(lagrange_multipliers,
                              -penalty_factor,
                               constraint_gaps);
        TSparse::Mult(r_transpose_relation_matrix,
                      lagrange_multipliers,
                      rhs_term);
        TSparse::UnaliasedAdd(rRhs,
                              -1.0,
                              rhs_term);
    }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::InitializeSolutionStep(typename TSparse::MatrixType& rLhs,
                                                                                  typename TSparse::VectorType& rSolution,
                                                                                  typename TSparse::VectorType& rRhs,
                                                                                  const std::size_t iIteration)
{
    KRATOS_TRY

    // Compute the constraint residuals.
    typename TSparse::VectorType constraint_residual(this->GetConstraintGapVector().size());
    TSparse::SetToZero(constraint_residual);
    TSparse::Mult(this->GetRelationMatrix(), rSolution, constraint_residual);
    TSparse::UnaliasedAdd(constraint_residual, -1.0, this->GetConstraintGapVector());

    // Update the RHS.
    typename TSparse::VectorType rhs_update(rRhs.size());
    TSparse::SetToZero(rhs_update);
    TSparse::InplaceMult(constraint_residual, -this->GetPenaltyFactor());
    TSparse::Mult(this->GetTransposeRelationMatrix(), constraint_residual, rhs_update);
    TSparse::UnaliasedAdd(rRhs, 1.0, rhs_update);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
typename AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Base::Status
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                                                           typename TSparse::VectorType& rSolution,
                                                                           typename TSparse::VectorType& rRhs,
                                                                           const std::size_t iIteration)
{
    KRATOS_TRY

    // Compute the constraint residuals.
    typename TSparse::VectorType constraint_residual(this->GetConstraintGapVector().size());
    TSparse::SetToZero(constraint_residual);
    TSparse::Mult(this->GetRelationMatrix(), rSolution, constraint_residual);
    TSparse::UnaliasedAdd(constraint_residual, -1.0, this->GetConstraintGapVector());

    // Decide whether to keep looping.
    const auto constraint_norm = TSparse::TwoNorm(constraint_residual);
    if (3 <= this->mpImpl->mVerbosity)
        std::cout << "AugmentedLagrangeConstraintAssembler: iteration " << iIteration
                  << " constraint residual " << (std::stringstream() << std::scientific << std::setprecision(8) << constraint_norm).str()
                  << "\n";

    const bool converged  = constraint_norm <= this->GetTolerance();
    const int max_iterations = this->GetValue(NL_ITERATION_NUMBER);

    if (converged) {
        if (2 <= mpImpl->mVerbosity)
            std::cout << "AugmentedLagrangeConstraintAssembler: constraints converged (residual "
                      << (std::stringstream() << std::scientific << std::setprecision(8) << constraint_norm).str()
                      << ")\n";
        return typename Base::Status {/*finished=*/true, /*converged=*/true};
    } /*if converged*/ else {
        if (static_cast<int>(iIteration) < max_iterations) {
            return typename Base::Status {/*finished=*/false, /*converged=*/false};
        } /*if iIteration < max_iterations*/ else {
            if (1 <= mpImpl->mVerbosity and not converged)
                std::cerr << "AugmentedLagrangeConstraintAssembler: constraints failed to converge (residual "
                          << (std::stringstream() << std::scientific << std::setprecision(8) << constraint_norm).str()
                          << ")\n";
            return typename Base::Status {/*finished=*/true, /*converged=*/false};
        }
    }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Finalize(typename TSparse::MatrixType& rLhs,
                                                                    typename TSparse::VectorType& rSolution,
                                                                    typename TSparse::VectorType& rRhs,
                                                                    typename Base::DofSet& rDofSet)
{
    KRATOS_TRY
    mpImpl->ErasePropagatedDirichletConditions(rDofSet.begin(), rDofSet.end());
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Clear()
{
    Base::Clear();
    mpImpl->mConstraintIdToIndexMap = decltype(mpImpl->mConstraintIdToIndexMap)();
    mpImpl->mMaybeTransposeRelationMatrix = decltype(mpImpl->mMaybeTransposeRelationMatrix)();
}


template <class TSparse, class TDense>
Parameters AugmentedLagrangeConstraintAssembler<TSparse,TDense>::GetDefaultParameters()
{
    return Parameters(R"({
"method" : "augmented_lagrange",
"penalty_factor" : 1e12,
"initial_lagrange_multiplier" : 0.0,
"tolerance" : 1e-8,
"max_iterations" : 1e1,
"verbosity" : 1
})");
}


template <class TSparse, class TDense>
typename TSparse::MatrixType&
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::GetTransposeRelationMatrix()
{
    KRATOS_TRY
    if (not mpImpl->mMaybeTransposeRelationMatrix.has_value()) {
        mpImpl->mMaybeTransposeRelationMatrix.emplace();
        SparseMatrixMultiplicationUtility::TransposeMatrix(mpImpl->mMaybeTransposeRelationMatrix.value(),
                                                           this->GetRelationMatrix());
    }
    return mpImpl->mMaybeTransposeRelationMatrix.value();
    KRATOS_CATCH("")
}


template class AugmentedLagrangeConstraintAssembler<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;

template class AugmentedLagrangeConstraintAssembler<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos
