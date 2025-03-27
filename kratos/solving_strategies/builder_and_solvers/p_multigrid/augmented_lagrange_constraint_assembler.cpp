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


void ProcessMasterSlaveConstraint(std::vector<std::size_t>& rConstraintIndices,
                                  std::vector<std::size_t>& rDofIds,
                                  MasterSlaveConstraint::MatrixType& rRelationMatrix,
                                  [[maybe_unused]] MasterSlaveConstraint::VectorType& rConstraintGaps,
                                  const MasterSlaveConstraint& rConstraint,
                                  const std::vector<std::size_t>& rSlaveIds,
                                  const std::vector<std::size_t>& rMasterIds,
                                  const std::unordered_map<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    // Constraint identifiers are the slave DoFs' IDs.
    rConstraintIndices.resize(rSlaveIds.size());
    std::transform(rSlaveIds.begin(),
                   rSlaveIds.end(),
                   rConstraintIndices.begin(),
                   [&rConstraintIdMap](std::size_t slave_id){
                        return rConstraintIdMap.at(slave_id).first;
                   });

    // Comb DoFs.
    rDofIds.resize(rMasterIds.size() + rSlaveIds.size());
    std::copy(rMasterIds.begin(), rMasterIds.end(), rDofIds.begin());
    std::copy(rSlaveIds.begin(), rSlaveIds.end(), rDofIds.begin() + rMasterIds.size());

    // Reconstruct the constraint equations.
    if (rRelationMatrix.size2()) {
        rRelationMatrix.resize(rRelationMatrix.size1(),
                               rRelationMatrix.size2() + rSlaveIds.size(),
                               /*preserve entries=*/true);
        for (std::size_t i_slave=0ul; i_slave<rSlaveIds.size(); ++i_slave) {
            // Fetch the number of constraint objects defining this constraint equation,
            // because the slave DoF's coefficient must add up to -1 after assembly.
            const auto constraint_object_count = rConstraintIdMap.at(rSlaveIds[i_slave]).second;
            rRelationMatrix(i_slave, rMasterIds.size() + i_slave) = -1.0 / constraint_object_count;
        } // for i_slave in range(rSlaveIds.size())
    } // if rRelationMatrix.size2()

    // No need to move the constraint gaps from RHS to LHS
    // because the slaves were moved to the RHS instead:
    // u^s = T u^m + b
    // =>
    // 0 = A [u^m \\ u^s] + b
    //std::transform(rConstraintGaps.begin(),
    //               rConstraintGaps.end(),
    //               rConstraintGaps.begin(),
    //               [](auto constraint_gap){return -constraint_gap;});
}


void ProcessMultifreedomConstraint(std::vector<std::size_t>& rConstraintIndices,
                                   std::vector<std::size_t>& rDofIds,
                                   const MasterSlaveConstraint& rConstraint,
                                   const std::vector<std::size_t>& rMasterIds,
                                   const std::unordered_map<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    const auto& r_constraint_labels = rConstraint.GetData().GetValue(CONSTRAINT_LABELS);
    rConstraintIndices.resize(r_constraint_labels.size());
    std::transform(r_constraint_labels.begin(),
                   r_constraint_labels.end(),
                   rConstraintIndices.begin(),
                   [&rConstraintIdMap](std::size_t constraint_label){
                        return rConstraintIdMap.at(constraint_label).first;
                   });

    rDofIds = rMasterIds;
}


template <class TSparse, class TDense, class TItDof>
void ApplyDirichletConditions(typename TSparse::MatrixType& rRelationMatrix,
                              typename TSparse::VectorType& rConstraintGaps,
                              TItDof itDofBegin,
                              [[maybe_unused]] TItDof itDofEnd,
                              const std::string& rConstraintAssemblerName,
                              int Verbosity)
{
    KRATOS_TRY

    /// @todo Make applying dirichlet conditions on penalty constraints more robust and efficient.

    std::size_t forced_dirichlet_dof_count,
                total_forced_dirichlet_dof_count = 0ul,
                iteration_count = 0ul;
    std::vector<LockObject> mutexes(rRelationMatrix.size2());

    do {
        forced_dirichlet_dof_count = 0ul;

        IndexPartition<std::size_t>(rRelationMatrix.size1()).for_each([&rRelationMatrix,
                                                                       &rConstraintGaps,
                                                                       &forced_dirichlet_dof_count,
                                                                       &mutexes,
                                                                       itDofBegin](std::size_t i_constraint){
            const auto i_entry_begin = rRelationMatrix.index1_data()[i_constraint];
            const auto i_entry_end = rRelationMatrix.index1_data()[i_constraint + 1];
            auto free_dof_count = i_entry_end - i_entry_begin;

            for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_dof = rRelationMatrix.index2_data()[i_entry];
                const Dof<typename TDense::DataType>& r_dof = *(itDofBegin + i_dof);
                KRATOS_DEBUG_ERROR_IF_NOT(i_dof == r_dof.EquationId());

                std::scoped_lock<LockObject> lock(mutexes[i_dof]);
                if (r_dof.IsFixed()) {
                    --free_dof_count;
                    // A Dirichlet condition is set on the current DoF, which also appears in
                    // a constraint equation => explicitly impose this condition on the constraint
                    // equation.
                    const typename TSparse::DataType constraint_coefficient = rRelationMatrix.value_data()[i_entry];
                    const typename TSparse::DataType dirichlet_condition = r_dof.GetSolutionStepValue();
                    AtomicAdd(rConstraintGaps[i_constraint],
                              static_cast<typename TSparse::DataType>(constraint_coefficient * dirichlet_condition));
                    rRelationMatrix.value_data()[i_entry] = static_cast<typename TSparse::DataType>(0);
                } // if r_dof.IsFixed()
            } // for i_entry in range(i_entry_begin, i_entry_end)

            if (free_dof_count == 1) {
                AtomicAdd(forced_dirichlet_dof_count, static_cast<std::size_t>(1));

                // If there's only 1 DoF left that does not have a Dirichlet condition
                // set on it, the constraint itself becomes equivalent to a Dirichlet
                // condition => find and constrain the last remaining free DoF.
                for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                    const auto i_dof = rRelationMatrix.index2_data()[i_entry];
                    Dof<typename TDense::DataType>& r_dof = *(itDofBegin + i_dof);

                    std::scoped_lock<LockObject> lock(mutexes[i_dof]);
                    if (r_dof.IsFree()) {
                        const typename TSparse::DataType constraint_coefficient = rRelationMatrix.value_data()[i_entry];
                        KRATOS_ERROR_IF_NOT(constraint_coefficient);
                        const typename TSparse::DataType constraint_gap = rConstraintGaps[i_constraint];
                        const auto dirichlet_condition = -constraint_gap / constraint_coefficient;

                        r_dof.GetSolutionStepValue() = dirichlet_condition;
                        r_dof.FixDof();
                        rRelationMatrix.value_data()[i_entry] = static_cast<typename TSparse::DataType>(0);
                        rConstraintGaps[i_constraint] = static_cast<typename TSparse::DataType>(0);
                        break;
                    }
                } // for i_entry in range(i_entry_begin, i_entry_end)
            } // if free_dof_count == 1
        }); // for i_constraint in range(rRelationMatrix.size1())

        ++iteration_count;
        total_forced_dirichlet_dof_count += forced_dirichlet_dof_count;
    } while (forced_dirichlet_dof_count);

    if (2 <= Verbosity && total_forced_dirichlet_dof_count) {
        std::cout << rConstraintAssemblerName << ": "
                  << "propagated Dirichlet conditions to " << total_forced_dirichlet_dof_count << " DoFs "
                  << "in " << iteration_count - 1 << " iterations\n";
    }
    KRATOS_CATCH("")
}


template <class TValue>
std::string FormattedResidual(TValue Residual)
{
    std::stringstream stream;
    stream << std::scientific
           << std::setprecision(8)
           << Residual;
    return stream.str();
}


} // namespace detail


template <class TSparse, class TDense>
struct AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Impl
{
    using Interface = AugmentedLagrangeConstraintAssembler<TSparse,TDense>;

    /// @brief A map associating slave IDs with constraint indices and the number of constraint objects referencing it.
    std::unordered_map<std::size_t,         //< identifier of the constraint equation (slave ID for MasterSlaveConstraint, CONSTRAINT_LABEL for MultifreedomConstraint)
                       std::pair<
                            std::size_t,    //< row index of the constraint equation in the relation matrix
                            std::size_t     //< number of constraint objects defining the equation
                       >
    > mConstraintIdMap;

    std::optional<typename TSparse::MatrixType> mMaybeTransposeRelationMatrix;

    int mVerbosity;
}; // struct AugmentedLagrangeConstraintAssembler::Impl


template <class TSparse, class TDense>

AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler() noexcept
    : AugmentedLagrangeConstraintAssembler(Parameters())
{}


template <class TSparse, class TDense>
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler(Parameters Settings)
    : AugmentedLagrangeConstraintAssembler(Settings, "unnamed")
{
}


template <class TSparse, class TDense>
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::AugmentedLagrangeConstraintAssembler(Parameters Settings,
                                                                                           std::string&& rInstanceName)
    : Base(ConstraintImposition::AugmentedLagrange, std::move(rInstanceName)),
      mpImpl(new Impl{/*mConstraintIdToIndexMap: */std::unordered_map<std::size_t,std::pair<std::size_t,std::size_t>>(),
                      /*mMaybeTransposeRelationMatrix: */std::optional<typename TSparse::MatrixType>(),
                      /*mVerbosity: */1})
{
    KRATOS_INFO_IF(this->Info(), 2 <= mpImpl->mVerbosity && !this->GetName().empty())
        << ": aliased to " << this->GetName();

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
    mpImpl->mConstraintIdMap.clear();

    {
        MasterSlaveConstraint::IndexType i_constraint = 0;
        MasterSlaveConstraint::EquationIdVectorType constraint_labels, master_ids;

        for (const auto& r_constraint : rConstraints) {
            r_constraint.EquationIdVector(constraint_labels, master_ids, rProcessInfo);

            if (constraint_labels.empty()) {
                // Constraint is a MultifreedomConstraint.
                const auto& r_constraint_labels = r_constraint.GetData().GetValue(CONSTRAINT_LABELS);
                constraint_labels.resize(r_constraint_labels.size());
                std::copy(r_constraint_labels.begin(),
                          r_constraint_labels.end(),
                          constraint_labels.begin());
            } /*if constraint_labels.empty()*/

            for (const auto constraint_label : constraint_labels) {
                const auto emplace_result = mpImpl->mConstraintIdMap.emplace(static_cast<std::size_t>(constraint_label),
                                                                             std::make_pair(i_constraint, 0ul));
                ++emplace_result.first->second.second; //< increment the number of constraint objects defining the constraint equation.
                if (emplace_result.second) ++i_constraint;
            } // for constraint_label in constraint_labels
        } // for r_constraint in rModelPart.MasterSlaveConstraints
    }

    {
        std::vector<std::unordered_set<IndexType>> indices(mpImpl->mConstraintIdMap.size());
        std::vector<LockObject> mutexes(mpImpl->mConstraintIdMap.size());

        struct TLS {
            MasterSlaveConstraint::EquationIdVectorType slaves, masters;
            std::vector<std::size_t> constraint_labels;
        };

        block_for_each(rConstraints,
                       TLS(),
                       [&mutexes, &indices, &rProcessInfo, this](const auto& r_constraint, TLS& r_tls) {
            r_constraint.EquationIdVector(r_tls.slaves, r_tls.masters, rProcessInfo);

            if (r_tls.slaves.empty()) {
                // Constraint is a MultifreedomConstraint.
                const auto& r_constraint_labels = r_constraint.GetData().GetValue(CONSTRAINT_LABELS);
                r_tls.constraint_labels.resize(r_constraint_labels.size());
                std::copy(r_constraint_labels.begin(),
                          r_constraint_labels.end(),
                          r_tls.constraint_labels.begin());
            } /*if r_tls.slaves.empty()*/ else {
                // Constraint is a MasterSlaveConstraint.
                r_tls.constraint_labels = r_tls.slaves;
                r_tls.masters.insert(r_tls.masters.end(),
                                     r_tls.slaves.begin(),
                                     r_tls.slaves.end());
            }

            for (const auto i_slave : r_tls.constraint_labels) {
                const auto i_constraint = mpImpl->mConstraintIdMap[static_cast<std::size_t>(i_slave)].first;
                std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                indices[i_constraint].insert(r_tls.masters.begin(), r_tls.masters.end());
            } // for i_slave in slave_ids
        }); // for r_constraint in rModelPart.MasterSlaveConstraints()

        MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                             rDofSet.size(),
                                                             this->GetRelationMatrix(),
                                                             /*EnsureDiagonal=*/false);
        TSparse::SetToZero(this->GetRelationMatrix());

        this->GetConstraintGapVector().resize(this->GetRelationMatrix().size1(), false);
        //TSparse::SetToZero(this->GetConstraintGapVector()); //< unnecessary
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
    std::vector<LockObject> mutexes(mpImpl->mConstraintIdMap.size());

    // Init.
    TSparse::SetToZero(this->GetRelationMatrix());
    TSparse::SetToZero(this->GetConstraintGapVector());

    struct TLS {
        MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
        std::vector<std::size_t> constraint_indices, dof_equation_ids, dof_index_array, reordered_dof_equation_ids;
        std::vector<MasterSlaveConstraint::MatrixType::value_type> relation_matrix_row;
        MasterSlaveConstraint::MatrixType relation_matrix;
        MasterSlaveConstraint::VectorType constraint_gaps;
    };

    // Constraint assembly.
    block_for_each(rConstraints,
                   TLS(),
                   [&mutexes, &rProcessInfo, this](const MasterSlaveConstraint& r_constraint, TLS& r_tls){
        if (r_constraint.IsActive()) {
            r_constraint.EquationIdVector(r_tls.slave_ids,
                                          r_tls.master_ids,
                                          rProcessInfo);
            r_constraint.CalculateLocalSystem(r_tls.relation_matrix,
                                              r_tls.constraint_gaps,
                                              rProcessInfo);

            // MasterSlaveConstraints and MultifreedomConstraints have to be handled separately.
            // MasterSlaveConstraint has both slave and master DoFs, and interprets the relation
            // matrix as a map from master DoFs to slave DoFs. The constraint equation it belongs
            // to is defined by the ID of the slave DoFs (MasterSlaveConstraints sharing slave DoFs
            // belong to the same constraint equation),
            // On the other hand, MultifreedomConstraint only fills the master DoFs, and provides a
            // relation matrix whose rows directly represent coefficients in the constraint equation.
            // Constraint equations are identified by the CONSTRAINT_LABELS variable stored in the
            // DataValueConstainer of the constraint object (constraints with identical CONSTRAINT_LABELS
            // belong to the same constraint equation).
            if (r_tls.slave_ids.empty()) {
                // The constraint is a MultifreedomConstraint.
                detail::ProcessMultifreedomConstraint(r_tls.constraint_indices,
                                                      r_tls.dof_equation_ids,
                                                      r_constraint,
                                                      r_tls.master_ids,
                                                      mpImpl->mConstraintIdMap);
            } /*if r_tls.slave_ids.empty()*/ else {
                // The constraint is a MasterSlaveConstraint.
                detail::ProcessMasterSlaveConstraint(r_tls.constraint_indices,
                                                     r_tls.dof_equation_ids,
                                                     r_tls.relation_matrix,
                                                     r_tls.constraint_gaps,
                                                     r_constraint,
                                                     r_tls.slave_ids,
                                                     r_tls.master_ids,
                                                     mpImpl->mConstraintIdMap);
            } // not r_tls.slave_ids.empty()

            r_tls.relation_matrix_row.resize(r_tls.relation_matrix.size2());

            // Assemble local rows into the global relation matrix.
            for (std::size_t i_row=0ul; i_row<r_tls.constraint_indices.size(); ++i_row) {
                const auto i_constraint = r_tls.constraint_indices[i_row];
                KRATOS_ERROR_IF(this->GetRelationMatrix().size1() <= i_constraint)
                    << "constraint index " << i_constraint
                    << " is out of bounds in relation matrix of size "
                    << "(" << this->GetRelationMatrix().size1() << "x" << this->GetRelationMatrix().size2() << ")";

                // Indirect sort the local relation matrix' row based on DoF IDs.
                // This step is required because the global relation matrix is in CSR format,
                // and Impl::MapRowContribution expects the column indices to be sorted.
                {
                    r_tls.dof_index_array.resize(r_tls.dof_equation_ids.size());
                    std::iota(r_tls.dof_index_array.begin(), r_tls.dof_index_array.end(), 0ul);
                    std::sort(r_tls.dof_index_array.begin(),
                              r_tls.dof_index_array.end(),
                              [&r_tls](const std::size_t i_left, const std::size_t i_right) -> bool {
                                  return r_tls.dof_equation_ids[i_left] < r_tls.dof_equation_ids[i_right];
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
                    r_tls.reordered_dof_equation_ids.resize(r_tls.dof_equation_ids.size());
                    std::transform(r_tls.dof_index_array.begin(),
                                   r_tls.dof_index_array.end(),
                                   r_tls.reordered_dof_equation_ids.begin(),
                                   [&r_tls](const auto i_column){
                                   return r_tls.dof_equation_ids[i_column];
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
                          static_cast<typename TSparse::DataType>(r_tls.constraint_gaps[i_row]));
            } // for i_row in range(r_tls.slave_ids.size)
        } // if r_constraint.IsActive
    }); // for r_constraint in rModelPart.MasterSlaveCosntraints

    // Propagate obvious dirichlet conditions.
    detail::ApplyDirichletConditions<TSparse,TDense>(this->GetRelationMatrix(),
                                                     this->GetConstraintGapVector(),
                                                     rDofSet.begin(),
                                                     rDofSet.end(),
                                                     this->GetName(),
                                                     mpImpl->mVerbosity);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Initialize(typename TSparse::MatrixType& rLhs,
                                                                      typename TSparse::VectorType& rRhs,
                                                                      typename Base::DofSet::iterator itDofBegin,
                                                                      typename Base::DofSet::iterator itDofEnd)
{
    KRATOS_TRY

    // Reset the relation matrix' transpose to make sure that propagated dirichlet
    // conditions are reflected in the transpose as well.
    mpImpl->mMaybeTransposeRelationMatrix.reset();

    //// Propagate obvious dirichlet conditions.
    //detail::ApplyDirichletConditions<TSparse,TDense>(this->GetRelationMatrix(),
    //                                                 this->GetConstraintGapVector(),
    //                                                 itDofBegin,
    //                                                 itDofEnd,
    //                                                 this->GetName(),
    //                                                 mpImpl->mVerbosity);

    // Fetch function-wide constants.
    using SparseUtils = SparseMatrixMultiplicationUtility;
    const typename TSparse::DataType penalty_factor = this->GetPenaltyFactor();
    const typename TSparse::DataType initial_lagrange_multiplier = this->GetInitialLagrangeMultiplier();

    // Apply initial penalty- and lagrange terms to the left- and right hand sides.
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
                                     lagrange_multipliers(this->GetRelationMatrix().size1(),
                                                          initial_lagrange_multiplier);

        TSparse::UnaliasedAdd(lagrange_multipliers,
                              -penalty_factor,
                              this->GetConstraintGapVector());

        TSparse::SetToZero(rhs_term);
        BalancedProduct<TSparse,TSparse,TSparse>(r_transpose_relation_matrix,
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
    typename TSparse::VectorType constraint_residual = this->GetConstraintGapVector();
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetRelationMatrix(), rSolution, constraint_residual);

    // Update the RHS.
    TSparse::InplaceMult(constraint_residual, -this->GetPenaltyFactor());
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetTransposeRelationMatrix(), constraint_residual, rRhs);

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
    typename TSparse::VectorType constraint_residual = this->GetConstraintGapVector();
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetRelationMatrix(), rSolution, constraint_residual);

    // Decide whether to keep looping.
    const auto constraint_norm = TSparse::TwoNorm(constraint_residual);
    if (3 <= this->mpImpl->mVerbosity)
        std::cout << this->GetName() << ": iteration " << iIteration
                  << " constraint residual " << detail::FormattedResidual(constraint_norm)
                  << "\n";

    const bool converged  = constraint_norm <= this->GetTolerance();
    const int max_iterations = this->GetValue(NL_ITERATION_NUMBER);

    if (converged) {
        if (2 <= mpImpl->mVerbosity)
            std::cout << this->GetName() << ": constraints converged (residual "
                      << detail::FormattedResidual(constraint_norm)
                      << ")\n";
        return typename Base::Status {/*finished=*/true, /*converged=*/true};
    } /*if converged*/ else {
        if (static_cast<int>(iIteration) < max_iterations) {
            return typename Base::Status {/*finished=*/false, /*converged=*/false};
        } /*if iIteration < max_iterations*/ else {
            if (1 <= mpImpl->mVerbosity && !converged)
                std::cerr << this->GetName() << ": constraints failed to converge (residual "
                          << detail::FormattedResidual(constraint_norm)
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
}


template <class TSparse, class TDense>
void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::Clear()
{
    Base::Clear();
    mpImpl->mConstraintIdMap = decltype(mpImpl->mConstraintIdMap)();
    mpImpl->mMaybeTransposeRelationMatrix = decltype(mpImpl->mMaybeTransposeRelationMatrix)();
}


template <class TSparse, class TDense>
Parameters AugmentedLagrangeConstraintAssembler<TSparse,TDense>::GetDefaultParameters()
{
    return Parameters(R"({
"method" : "augmented_lagrange",
"penalty_factor" : 1e12,
"initial_lagrange_multiplier" : 0.0,
"tolerance" : 1e-6,
"max_iterations" : 1e1,
"verbosity" : 1
})");
}


template <class TSparse, class TDense>
typename TSparse::MatrixType&
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::GetTransposeRelationMatrix()
{
    KRATOS_TRY
    if (!mpImpl->mMaybeTransposeRelationMatrix.has_value()) {
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
