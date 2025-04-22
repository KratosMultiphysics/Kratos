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
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp" // ProcessMasterSlaveConstraint, ProcessMultifreedomConstraint
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // Scaling
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility

// STL includes
#include <filesystem> // std::filesystem::current_path


namespace Kratos {


namespace detail {


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

    std::unique_ptr<Scaling> mpPenaltyFunctor;

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
                      /*mpScaling*/ nullptr,
                      /*mVerbosity: */1})
{
    KRATOS_INFO_IF(this->Info(), 2 <= mpImpl->mVerbosity && !this->GetName().empty())
        << ": aliased to " << this->GetName();

    // Parse the penalty factor and validate other settings.
    Parameters default_parameters = this->GetDefaultParameters();
    Parameters default_penalty_factor = default_parameters["penalty_factor"].Clone();
    default_parameters.RemoveValue("penalty_factor");
    std::optional<Parameters> maybe_penalty_factor = Settings.Has("penalty_factor") ? Settings["penalty_factor"].Clone() : std::optional<Parameters>();
    if (maybe_penalty_factor.has_value()) Settings.RemoveValue("penalty_factor");
    Settings.ValidateAndAssignDefaults(default_parameters);
    Settings.AddValue("penalty_factor", maybe_penalty_factor.has_value() ? maybe_penalty_factor.value() : default_penalty_factor);
    mpImpl->mpPenaltyFunctor = std::make_unique<Scaling>(Settings["penalty_factor"]);

    // Parse other algorithmic settings.
    Vector algorithmic_parameters(2);
    algorithmic_parameters[0] = Settings["initial_lagrange_multiplier"].Get<double>();
    algorithmic_parameters[1] = Settings["tolerance"].Get<double>();
    this->SetValue(this->GetAlgorithmicParametersVariable(), algorithmic_parameters);
    this->SetValue(NL_ITERATION_NUMBER, Settings["max_iterations"].Get<int>());

    // Parse miscellaneous settings.
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

    this->Clear();
    if (rConstraints.empty()) {
        this->GetRelationMatrix() = typename TSparse::MatrixType(0, rLhs.size2(), 0);
        this->GetRelationMatrix().index1_data()[0] = 0;
        this->GetRelationMatrix().set_filled(1, 0);
        return;
    }

    detail::MakeRelationTopology<TSparse,TDense>(rDofSet.size(),
                                                 rConstraints,
                                                 rProcessInfo,
                                                 this->GetRelationMatrix(),
                                                 this->GetConstraintGapVector(),
                                                 mpImpl->mConstraintIdMap);

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
    detail::AssembleRelationMatrix<TSparse,TDense>(rConstraints,
                                                   rProcessInfo,
                                                   this->GetRelationMatrix(),
                                                   this->GetConstraintGapVector(),
                                                   mpImpl->mConstraintIdMap);

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

    // Reset the relation matrix' transpose to make sure that propagated dirichlet
    // conditions are reflected in the transpose as well.
    mpImpl->mMaybeTransposeRelationMatrix.reset();

    // Fetch function-wide constants.
    using SparseUtils = SparseMatrixMultiplicationUtility;
    mpImpl->mpPenaltyFunctor->template Cache<TSparse>(rLhs);
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
                                                                                  typename TSparse::VectorType& rRhs)
{
    KRATOS_TRY

    // Compute the constraint residuals.
    typename TSparse::VectorType constraint_residual = this->GetConstraintGapVector();
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetRelationMatrix(), rSolution, constraint_residual);

    // Update the RHS.
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetTransposeRelationMatrix(),
                                             constraint_residual,
                                             rRhs,
                                             static_cast<typename TSparse::DataType>(-this->GetPenaltyFactor()));

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
bool
AugmentedLagrangeConstraintAssembler<TSparse,TDense>::FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                                                           typename TSparse::VectorType& rSolution,
                                                                           typename TSparse::VectorType& rRhs,
                                                                           PMGStatusStream::Report& rReport,
                                                                           PMGStatusStream& rStream)
{
    KRATOS_TRY

    // Compute the constraint residuals.
    typename TSparse::VectorType constraint_residual = this->GetConstraintGapVector();
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetRelationMatrix(), rSolution, constraint_residual);

    // Update status.
    rReport.maybe_constraint_residual = TSparse::TwoNorm(constraint_residual);
    rReport.constraints_converged = rReport.maybe_constraint_residual.value() <= this->GetTolerance();

    // Decide on whether to keep looping and log state.
    if (rReport.constraints_converged) {
        rStream.Submit(rReport.Tag(2), mpImpl->mVerbosity);
        return true;
    } /*if converged*/ else {
        const int max_iterations = this->GetValue(NL_ITERATION_NUMBER);
        if (static_cast<int>(rReport.constraint_iteration) + 1 < max_iterations) {
            rStream.Submit(rReport.Tag(3), mpImpl->mVerbosity);
            return false;
        } /*if iIteration < max_iterations*/ else {
            rStream.Submit(rReport.Tag(2), mpImpl->mVerbosity);
            return true;
        }
    } /*if not converged*/
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
"penalty_factor" : "norm",
"initial_lagrange_multiplier" : 0.0,
"tolerance" : 1e-6,
"max_iterations" : 1e1,
"verbosity" : 1
})");
}


template <class TSparse, class TDense>
typename TSparse::DataType AugmentedLagrangeConstraintAssembler<TSparse,TDense>::GetPenaltyFactor() const
{
    return mpImpl->mpPenaltyFunctor->Evaluate();
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
