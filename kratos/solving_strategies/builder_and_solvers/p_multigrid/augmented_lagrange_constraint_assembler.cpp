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


namespace Kratos {


template <class TSparse, class TDense>

AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler() noexcept
    : AugmentedLagrangeConstraintAssembler(Parameters())
{}


template <class TSparse, class TDense>
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler(Parameters Settings)
    : Base(ConstraintImposition::AugmentedLagrange),
      mSlaveToConstraintMap(),
      mMaybeTransposeRelationMatrix(),
      mVerbosity(1)
{
    Settings.ValidateAndAssignDefaults(AugmentedLagrangeConstraintAssembler::GetDefaultParameters());

    Vector algorithmic_parameters(3);
    algorithmic_parameters[0] = Settings["penalty_factor"].Get<double>();
    algorithmic_parameters[1] = Settings["initial_lagrange_multiplier"].Get<double>();
    algorithmic_parameters[2] = Settings["tolerance"].Get<double>();
    this->SetValue(this->GetAlgorithmicParametersVariable(), algorithmic_parameters);
    this->SetValue(NL_ITERATION_NUMBER, Settings["max_iterations"].Get<int>());
    this->mVerbosity = Settings["verbosity"].Get<int>();
}


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

    // Construct a map that associates slave indices with constraint equation indices.
    mSlaveToConstraintMap.clear();

    {
        MasterSlaveConstraint::IndexType i_constraint = 0;
        MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
        for (const auto& r_constraint : rConstraints) {
            r_constraint.EquationIdVector(slave_ids, master_ids, rProcessInfo);
            for (const auto i_slave : slave_ids) {
                const auto emplace_result = mSlaveToConstraintMap.emplace(i_slave, i_constraint);
                if (emplace_result.second) ++i_constraint;
            } // for i_slave in slave_ids
        } // for r_constraint in rModelPart.MasterSlaveConstraints
    }

    {
        struct TLS {
            Element::EquationIdVectorType master_ids, slave_ids;
        }; // struct TLS

        std::vector<std::unordered_set<IndexType>> indices(mSlaveToConstraintMap.size());
        std::vector<LockObject> mutexes(indices.size());

        block_for_each(rConstraints,
                        TLS(),
                        [&mutexes, &indices, &rProcessInfo, this](const auto& r_constraint, TLS& r_tls){
            r_constraint.EquationIdVector(r_tls.slave_ids,
                                          r_tls.master_ids,
                                          rProcessInfo);

            for (const auto i_slave : r_tls.slave_ids) {
                const auto i_constraint = mSlaveToConstraintMap[i_slave];
                std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                indices[i_constraint].insert(r_tls.master_ids.begin(), r_tls.master_ids.end());
                indices[i_constraint].insert(i_slave);
            } // for i_slave in slave_ids
        }); // for r_constraint in rModelPart.MasterSlaveConstraints()

        MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                             rDofSet.size(),
                                                             this->GetRelationMatrix());

        this->GetConstraintGapVector().resize(mSlaveToConstraintMap.size(), false);
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
                                                                    const typename Base::DofSet& rDofSet,
                                                                    const bool AssembleLhs,
                                                                    const bool AssembleRhs)
{
    KRATOS_TRY

    /// @todo @p AssembleLhs and @p AssembleRhs are currently ignored,
    ///       and everything gets assembled. Think it through what needs
    ///       to be done and what can be omitted for each case. - matekelemen

    // Function-wide variables.
    std::vector<LockObject> mutexes(mSlaveToConstraintMap.size());

    // Init.
    TSparse::SetToZero(this->GetRelationMatrix());
    TSparse::SetToZero(this->GetConstraintGapVector());

    // Thread-local storage for constraint assembly.
    struct TLS {
        MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
        typename TDense::MatrixType local_relation_matrix;
        typename TDense::VectorType local_constraint_gap_vector;
    }; // struct TLS

    // Constraint assembly.
    block_for_each(rConstraints,
                    TLS(),
                    [&mutexes, &rProcessInfo, this](const MasterSlaveConstraint& r_constraint, TLS& r_tls){
        if (r_constraint.IsActive()) {
            r_constraint.EquationIdVector(r_tls.slave_ids,
                                          r_tls.master_ids,
                                          rProcessInfo);
            r_constraint.CalculateLocalSystem(r_tls.local_relation_matrix,
                                              r_tls.local_constraint_gap_vector,
                                              rProcessInfo);

            // Slaves also need to be part of the relation matrix in this formulation,
            // so the local relation matrix one must be appended.
            if (r_tls.local_relation_matrix.size2()) {
                r_tls.local_relation_matrix.resize(r_tls.local_relation_matrix.size1(),
                                                   r_tls.local_relation_matrix.size2() + 1,
                                                   true);
                for (std::size_t i_row=0; i_row<r_tls.local_relation_matrix.size1(); ++i_row) {
                    r_tls.local_relation_matrix(i_row, r_tls.local_relation_matrix.size2() - 1) = -1;
                } // for i_row in range(local_relation_matrix.size2)
            } // if local_relation_matrix.size2

            // Assemble local rows into the global relation matrix.
            for (std::size_t i_row=0ul; i_row<r_tls.slave_ids.size(); ++i_row) {
                const auto it_constraint_index = mSlaveToConstraintMap.find(r_tls.slave_ids[i_row]);
                KRATOS_ERROR_IF(it_constraint_index == mSlaveToConstraintMap.end())
                    << "slave " << i_row << " (Dof " << r_tls.slave_ids[i_row] << ") "
                    << "in constraint " << r_constraint.Id() << " "
                    << "is not associated with a constraint index.";

                const std::size_t i_constraint = it_constraint_index->second;
                KRATOS_ERROR_IF(this->GetRelationMatrix().size1() <= i_constraint);

                auto dof_ids = r_tls.master_ids;
                dof_ids.push_back(r_tls.slave_ids[i_row]);
                KRATOS_ERROR_IF_NOT(dof_ids.size() == r_tls.local_relation_matrix.size2());

                // Indirect sort the local relation matrix' row based on DoF IDs.
                // This step is required because the global relation matrix is in CSR format,
                // and Impl::MapRowContribution expects the column indices to be sorted.
                {
                    std::vector<std::size_t> index_array(dof_ids.size());
                    std::iota(index_array.begin(), index_array.end(), 0ul);
                    std::sort(index_array.begin(),
                              index_array.end(),
                              [&dof_ids](const std::size_t i_left, const std::size_t i_right) -> bool {
                                  return dof_ids[i_left] < dof_ids[i_right];
                              });

                    // Reorder the current row in the local relation matrix.
                    std::vector<MasterSlaveConstraint::MatrixType::value_type> row(r_tls.local_relation_matrix.size2());
                    std::transform(index_array.begin(),
                                   index_array.end(),
                                   row.begin(),
                                   [&r_tls, i_row](const std::size_t i_column){
                                   return r_tls.local_relation_matrix(i_row, i_column);
                                   });
                    std::copy(row.begin(),
                              row.end(),
                              r_tls.local_relation_matrix.data().begin() + i_row * r_tls.local_relation_matrix.size2());

                    // Reorder DoF indices.
                    auto swap_dof_ids = dof_ids;
                    std::transform(index_array.begin(),
                                   index_array.end(),
                                   swap_dof_ids.begin(),
                                   [&dof_ids](const auto i_column){
                                   return dof_ids[i_column];
                                   });
                    dof_ids.swap(swap_dof_ids);
                }

                {
                    std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                    MapRowContribution<TSparse,TDense>(this->GetRelationMatrix(),
                                                        r_tls.local_relation_matrix,
                                                        i_constraint,
                                                        i_row,
                                                        dof_ids);
                }
            } // for i_row in range(r_tls.slave_ids.size)
        } // if r_constraint.IsActive
    }); // for r_constraint in rModelPart.MasterSlaveCosntraints
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
            << "the transpose of the relation matrix is uninitialized or out of date";

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
typename AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Base::Status
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                                                           typename TSparse::VectorType& rSolution,
                                                                           typename TSparse::VectorType& rRhs,
                                                                           const std::size_t iIteration)
{
    KRATOS_TRY
    const int max_iterations = this->GetValue(NL_ITERATION_NUMBER);
    // Compute the constraint residuals.
    typename TSparse::VectorType constraint_residual(this->GetConstraintGapVector().size());
    TSparse::SetToZero(constraint_residual);

    KRATOS_TRY
    TSparse::Mult(this->GetRelationMatrix(), rSolution, constraint_residual);
    TSparse::UnaliasedAdd(constraint_residual, -1.0, this->GetConstraintGapVector());
    KRATOS_CATCH("")

    // Decide whether to update the RHS or terminate the loop.
    const auto constraint_norm = TSparse::TwoNorm(constraint_residual);
    if (3 <= this->mVerbosity) {
        std::stringstream residual_stream;
        residual_stream << std::setprecision(8) << std::scientific<< constraint_norm;
        std::cout << "AugmentedLagrangeConstraintAssembler: iteration " << iIteration
                  << " constraint residual " << residual_stream.str() << "\n";
    }

    if (this->GetTolerance() < constraint_norm && static_cast<int>(iIteration) < max_iterations) {
        KRATOS_TRY
        typename TSparse::VectorType rhs_update(rRhs.size());
        TSparse::SetToZero(rhs_update);
        TSparse::InplaceMult(constraint_residual, -this->GetPenaltyFactor());
        TSparse::Mult(this->GetTransposeRelationMatrix(), constraint_residual, rhs_update);
        TSparse::UnaliasedAdd(rRhs, 1.0, rhs_update);
        return typename Base::Status {/*finished=*/false, /*converged=*/false};
        KRATOS_CATCH("")
    } else {
        return typename Base::Status {/*finished=*/true, /*converged=*/constraint_norm <= this->GetTolerance()};
    }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Clear()
{
    Base::Clear();
    mSlaveToConstraintMap = decltype(mSlaveToConstraintMap)();
    mMaybeTransposeRelationMatrix = decltype(mMaybeTransposeRelationMatrix)();
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
    if (not mMaybeTransposeRelationMatrix.has_value()) {
        mMaybeTransposeRelationMatrix.emplace();
        SparseMatrixMultiplicationUtility::TransposeMatrix(mMaybeTransposeRelationMatrix.value(),
                                                           this->GetRelationMatrix());
    }
    return mMaybeTransposeRelationMatrix.value();
    KRATOS_CATCH("")
}


template class AugmentedLagrangeConstraintAssembler<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;

template class AugmentedLagrangeConstraintAssembler<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos
