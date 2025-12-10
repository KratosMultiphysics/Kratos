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
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp" // ProcessMasterSlaveConstraint, ProcessMultifreedomConstraint, detail::MakeRelationTopology
#include "factories/linear_solver_factory.h" // LinearSolver, LinearSolverFactory
#include "includes/define.h" // KRATOS_TRY, KRATOS_CATCH
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility


namespace Kratos {


template <class TSparse, class TDense>
struct MasterSlaveConstraintAssembler<TSparse,TDense>::Impl
{
    std::vector<const Dof<typename TDense::DataType>*> mSlaves;

    std::unique_ptr<Scaling> mpDiagonalScaling;

    /// @brief A map associating slave IDs with constraint indices and the number of constraint objects referencing it.
    CSRHashMap<std::size_t,         //< identifier of the constraint equation (slave ID for MasterSlaveConstraint, CONSTRAINT_LABEL for MultifreedomConstraint)
               std::pair<
                    std::size_t,    //< row index of the constraint equation in the relation matrix
                    std::size_t     //< number of constraint objects defining the equation
              >
    > mConstraintIdMap;

    /// @brief Index map from the dependent space to the independent space.
    std::vector<std::optional<typename TSparse::IndexType>> mIndependentMap;

    typename TSparse::MatrixType mMasterSlaveRelations;

    struct LinearSystem {
        typename TSparse::MatrixType mLhs;
        typename TSparse::VectorType mSolution, mRhs;
        typename Base::DofSet mDofSet;
    }; // struct LinearSystem

    std::optional<LinearSystem> mDependentSystem;

    int mVerbosity;

    std::pair<typename TSparse::MatrixType&, LinearSolver<TSparse,TDense>&> GetIndependentTransform()
    {
        if (!mIndependentTransform.has_value()) {
            KRATOS_TRY
            typename TSparse::MatrixType transpose_master_slave_relations;
            SparseMatrixMultiplicationUtility::TransposeMatrix(transpose_master_slave_relations, mMasterSlaveRelations);
            mIndependentTransform.emplace();

            typename TSparse::MatrixType& r_lhs = mIndependentTransform.value().first;
            SparseMatrixMultiplicationUtility::MatrixMultiplication(transpose_master_slave_relations,
                                                                    mMasterSlaveRelations,
                                                                    r_lhs);

            mIndependentTransform.value().second = LinearSolverFactory<TSparse,TDense>().Create(mIndependentSolverSettings);
            typename TSparse::VectorType dummy;
            mIndependentTransform.value().second->Initialize(r_lhs, dummy, dummy);
            mIndependentTransform.value().second->InitializeSolutionStep(r_lhs, dummy, dummy);
            KRATOS_CATCH("")
        }

        return std::pair<typename TSparse::MatrixType&,LinearSolver<TSparse,TDense>&>(mIndependentTransform.value().first, *mIndependentTransform.value().second);
    }

    std::pair<typename TSparse::MatrixType&, LinearSolver<TSparse,TDense>&> GetDependentTransform()
    {
        if (!mDependentTransform.has_value()) {
            KRATOS_TRY
            typename TSparse::MatrixType transpose_master_slave_relations;
            SparseMatrixMultiplicationUtility::TransposeMatrix(transpose_master_slave_relations, mMasterSlaveRelations);
            mDependentTransform.emplace();

            typename TSparse::MatrixType& r_lhs = mDependentTransform.value().first;
            SparseMatrixMultiplicationUtility::MatrixMultiplication(mMasterSlaveRelations,
                                                                    transpose_master_slave_relations,
                                                                    r_lhs);

            mDependentTransform.value().second = LinearSolverFactory<TSparse,TDense>().Create(mDependentSolverSettings);
            typename TSparse::VectorType dummy;
            mDependentTransform.value().second->Initialize(r_lhs, dummy, dummy);
            mDependentTransform.value().second->InitializeSolutionStep(r_lhs, dummy, dummy);
            KRATOS_CATCH("")
        }

        return std::pair<typename TSparse::MatrixType&,LinearSolver<TSparse,TDense>&>(mDependentTransform.value().first, *mDependentTransform.value().second);
    }

    Parameters mDependentSolverSettings, mIndependentSolverSettings;

    std::optional<std::pair<
        typename TSparse::MatrixType,
        typename LinearSolver<TSparse,TDense>::Pointer
    >> mDependentTransform, mIndependentTransform;
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
    mpImpl->mDependentSolverSettings = Settings["dependent_solver_settings"];
    mpImpl->mIndependentSolverSettings = Settings["independent_solver_settings"];
    mpImpl->mVerbosity = Settings["verbosity"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
MasterSlaveConstraintAssembler<TSparse,TDense>::~MasterSlaveConstraintAssembler() = default;


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::AllocateConstraints(PointerVectorSet<MasterSlaveConstraint,IndexedObject>::const_iterator itConstraintBegin,
                                                                         PointerVectorSet<MasterSlaveConstraint,IndexedObject>::const_iterator itConstraintEnd,
                                                                         const ProcessInfo& rProcessInfo,
                                                                         typename Base::DofSet::const_iterator itDofBegin,
                                                                         typename Base::DofSet::const_iterator itDofEnd)
{
    KRATOS_TRY
    this->Clear();
    mpImpl->mDependentSystem.emplace();
    KRATOS_CATCH("")

    KRATOS_TRY
    const typename TSparse::IndexType dependent_system_size = std::distance(itDofBegin, itDofEnd);

    if (itConstraintBegin == itConstraintEnd) {
        this->GetRelationMatrix() = typename TSparse::MatrixType(0, dependent_system_size, 0);
        this->GetRelationMatrix().index1_data()[0] = 0;
        std::fill(this->GetRelationMatrix().index1_data().begin(),
                  this->GetRelationMatrix().index1_data().end(),
                  static_cast<typename TSparse::IndexType>(0));
        this->GetRelationMatrix().set_filled(1, 0);
    } else {
        detail::MakeRelationTopology<TSparse,TDense>(dependent_system_size,
                                                     itConstraintBegin,
                                                     itConstraintEnd,
                                                     rProcessInfo,
                                                     this->GetRelationMatrix(),
                                                     this->GetConstraintGapVector(),
                                                     mpImpl->mConstraintIdMap);
    }

    // Compute the modified relation matrix (T) for master-slave imposition.
    // This involves a partitioning of the DoFs (u) into:
    // 0) DoFs that don't participate in any constraints    (w)
    // 1) master DoFs                                       (m)
    // 2) slave DoFs                                        (s)
    //
    // Assuming the constraint equations take the following form:
    // A u + g = 0
    //
    // After partitioning u:
    // + ----------- +  + - +
    // | 0  A_m  A_s |  | w | + g = 0
    // + ----------- +  | m |
    //                  | s |
    //                  + - +
    //
    // The modified relation matrix then expresses the dependent set of DoFs (w, m, s)
    // in terms of the independent set (w, m):
    // + - +   + --------------- + + - +   + ----------- +
    // | w | = | I       0       | | w | + |      0      |
    // | m |   | 0       I       | | m |   |      0      |
    // | s |   | 0  -A_s^{-1}A_m | + - +   | -A_s^{-1} g |
    // + - +   + --------------- +         + ----------- +
    //
    // Note that in general, the partitioning is not contiguous so you won't
    // get the nice block structure depicted above.

    KRATOS_ERROR_IF_NOT(this->GetRelationMatrix().size1() <= dependent_system_size)
        << "Number of constraint equations (" << this->GetRelationMatrix().size1() << ") "
        << "is greater than the number of DoFs (" << dependent_system_size << ")";

    KRATOS_ERROR_IF(this->GetRelationMatrix().nnz() < this->GetRelationMatrix().size1())
        << "Assembled relation matrix has less entries (" << this->GetRelationMatrix().nnz() << ") "
        << "than constraint equations (" << this->GetRelationMatrix().size1() << ")";

    const auto master_slave_system_size = dependent_system_size - this->GetRelationMatrix().size1();
    const auto master_slave_entry_count = master_slave_system_size              //< 1s on the diagonal of independent DoFs
                                        + this->GetRelationMatrix().nnz()       //< entries in the relation matrix
                                        - this->GetRelationMatrix().size1();    //< one unique slave per constraint equation
    mpImpl->mMasterSlaveRelations = typename TSparse::MatrixType(dependent_system_size,
                                                                 master_slave_system_size,
                                                                 master_slave_entry_count);

    // Compute row extents.
    // This happens in 2 stages: write the number of entries in each row into
    // an intermediate vector, then compute the cumulative sum onto the row
    // extents.
    std::vector<std::optional<typename TSparse::IndexType>> row_sizes(dependent_system_size);

    mpImpl->mIndependentMap.clear();
    mpImpl->mIndependentMap.resize(dependent_system_size);

    //block_for_each(mpImpl->mConstraintIdMap, [&row_sizes, this] (const auto& r_constraint_equation_info) {
    for (const auto& r_constraint_equation_info : mpImpl->mConstraintIdMap) {
        const auto i_slave = r_constraint_equation_info.first;
        mpImpl->mSlaves.push_back(&*(itDofBegin + i_slave));
        const auto i_constraint = r_constraint_equation_info.second.first;
        const auto i_entry_begin = this->GetRelationMatrix().index1_data()[i_constraint];
        const auto i_entry_end = this->GetRelationMatrix().index1_data()[i_constraint + 1];
        KRATOS_ERROR_IF_NOT(i_entry_end - i_entry_begin) << "empty constraint equation at index " << i_constraint;
        const auto master_slave_row_size = i_entry_end - i_entry_begin - 1;
        row_sizes[i_slave] = master_slave_row_size;
    }
    //});

    mpImpl->mMasterSlaveRelations.index1_data()[0] = static_cast<typename TSparse::IndexType>(0);
    typename TSparse::IndexType i_independent_dof = 0;
    for (typename TSparse::IndexType i_row=0; i_row<dependent_system_size; ++i_row) {
        if (row_sizes[i_row].has_value()) {
            mpImpl->mMasterSlaveRelations.index1_data()[i_row + 1] = mpImpl->mMasterSlaveRelations.index1_data()[i_row]
                                                                   + row_sizes[i_row].value();
        } else {
            mpImpl->mIndependentMap[i_row] = i_independent_dof++;
            mpImpl->mMasterSlaveRelations.index1_data()[i_row + 1] = mpImpl->mMasterSlaveRelations.index1_data()[i_row]
                                                                   + 1;
        }
    } // for i_row in range(dependent_system_size)
    KRATOS_ERROR_IF_NOT(mpImpl->mMasterSlaveRelations.index1_data()[dependent_system_size] == master_slave_entry_count)
        << "internal size mismatch: "
        << mpImpl-> mMasterSlaveRelations.index1_data()[dependent_system_size]
        << " != "
        << master_slave_entry_count;

    // Copy column indices.
    IndexPartition<typename TSparse::IndexType>(dependent_system_size).for_each([this] (auto i_dof) {
        const auto it_constraint_equation_info = mpImpl->mConstraintIdMap.find(i_dof);
        if (it_constraint_equation_info == mpImpl->mConstraintIdMap.end()) {
            // DoF is not a slave.
            const auto i_entry = mpImpl->mMasterSlaveRelations.index1_data()[i_dof];
            mpImpl->mMasterSlaveRelations.index2_data()[i_entry] = mpImpl->mIndependentMap[i_dof].value();
            mpImpl->mMasterSlaveRelations.value_data()[i_entry] = static_cast<typename TSparse::DataType>(1);
        } else {
            // Dof is a slave.
            const auto i_slave = it_constraint_equation_info->first;
            const auto i_constraint = it_constraint_equation_info->second.first;
            const auto i_relation_entry_begin = this->GetRelationMatrix().index1_data()[i_constraint];
            const auto i_relation_entry_end = this->GetRelationMatrix().index1_data()[i_constraint + 1];
            const auto i_master_slave_entry_begin = mpImpl->mMasterSlaveRelations.index1_data()[i_slave];

            typename TSparse::IndexType i_master_slave_entry = i_master_slave_entry_begin;
            for (typename TSparse::IndexType i_relation_entry=i_relation_entry_begin; i_relation_entry<i_relation_entry_end; ++i_relation_entry) {
                const auto i_relation_column = this->GetRelationMatrix().index2_data()[i_relation_entry];
                if (i_relation_column != i_slave) {
                    mpImpl->mMasterSlaveRelations.index2_data()[i_master_slave_entry] = mpImpl->mIndependentMap[i_relation_column].value();
                    ++i_master_slave_entry;
                }
            } // for i_relation_entry in range(i_relation_entry_begin, i_relation_entry_end)
        }
    });

    mpImpl->mMasterSlaveRelations.set_filled(dependent_system_size + 1, master_slave_entry_count);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Assemble(const typename Base::ConstraintArray& rConstraints,
                                                              const ProcessInfo& rProcessInfo,
                                                              [[maybe_unused]] typename Base::DofSet::const_iterator itDofBegin,
                                                              [[maybe_unused]] typename Base::DofSet::const_iterator itDofEnd,
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
        this->GetHessian(),
        this->GetConstraintGapVector(),
        mpImpl->mConstraintIdMap);

    // Map relation matrix to master-slave relation matrix.
    for (const auto& r_constraint_equation_info : mpImpl->mConstraintIdMap) {
        const auto i_slave = r_constraint_equation_info.first;
        const auto i_constraint = r_constraint_equation_info.second.first;
        const auto i_relation_entry_begin = this->GetRelationMatrix().index1_data()[i_constraint];
        const auto i_relation_entry_end = this->GetRelationMatrix().index1_data()[i_constraint + 1];
        const auto i_master_slave_entry_begin = mpImpl->mMasterSlaveRelations.index1_data()[i_slave];

        // Find the slave coefficient.
        auto it_slave = std::lower_bound(this->GetRelationMatrix().index2_data().begin() + i_relation_entry_begin,
                                         this->GetRelationMatrix().index2_data().begin() + i_relation_entry_end,
                                         i_slave);
        KRATOS_ERROR_IF(it_slave == this->GetRelationMatrix().index2_data().begin() + i_relation_entry_end
                     || *it_slave != i_slave)
            << "no entry in relation matrix for slave DoF " << i_slave << " "
            << "of constraint equation " << i_constraint;
        const typename TSparse::IndexType i_slave_in_master_slave = std::distance(this->GetRelationMatrix().index2_data().begin(),
                                                                                  it_slave);
        const typename TSparse::DataType slave_coefficient = this->GetRelationMatrix().value_data()[i_slave_in_master_slave];
        KRATOS_ERROR_IF_NOT(slave_coefficient) << "coefficient of slave DoF " << i_slave << " in constraint equation "
                                               << i_constraint << " is null.";
        const typename TSparse::DataType coefficient_scale = static_cast<typename TSparse::DataType>(-1) / slave_coefficient;

        // Assemble entries from the relation matrix to the master-slave relation matrix.
        typename TSparse::IndexType i_master_slave_entry = i_master_slave_entry_begin;
        for (typename TSparse::IndexType i_relation_entry=i_relation_entry_begin; i_relation_entry<i_relation_entry_end; ++i_relation_entry) {
            const auto i_relation_column = this->GetRelationMatrix().index2_data()[i_relation_entry];
            if (i_relation_column != i_slave) {
                mpImpl->mMasterSlaveRelations.value_data()[i_master_slave_entry] = coefficient_scale * this->GetRelationMatrix().value_data()[i_relation_entry];
                ++i_master_slave_entry;
            }
        } // for i_entry in range(i_relation_entry_begin, i_relation_entry_end)
    }

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Initialize(typename TSparse::MatrixType& rLhs,
                                                                typename TSparse::VectorType& rSolution,
                                                                typename TSparse::VectorType& rRhs,
                                                                typename Base::DofSet& rDofs)
{
    KRATOS_TRY

    // Output the relation matrix and constraint gaps if the verbosity is high enough.
    if (4 <= mpImpl->mVerbosity) {
        KRATOS_INFO(this->GetName())
            << "write relation matrix to "
            << (std::filesystem::current_path() / "relation_matrix.mm")
            << "\n";
        TSparse::WriteMatrixMarketMatrix("relation_matrix.mm", mpImpl->mMasterSlaveRelations, false);

        KRATOS_INFO(this->GetName())
            << "write constraint gaps to "
            << (std::filesystem::current_path() / "constraint_gaps.mm")
            << "\n";
        TSparse::WriteMatrixMarketVector("constraint_gaps.mm", this->GetConstraintGapVector());
    } // if 4 <= verbosity

    // Compute the transposed matrix of the global relation matrix
    {
        // Initialize independent RHS.
        mpImpl->mDependentSystem.value().mRhs.resize(mpImpl->mMasterSlaveRelations.size2());
        TSparse::SetToZero(mpImpl->mDependentSystem.value().mRhs);

        // Initialize independent solution.
        mpImpl->mDependentSystem.value().mSolution.resize(mpImpl->mMasterSlaveRelations.size2());
        TSparse::SetToZero(mpImpl->mDependentSystem.value().mSolution);

        // Storage for an intermediate matrix transpose(relation_matrix) * lhs_matrix
        typename TSparse::MatrixType left_multiplied_lhs(mpImpl->mMasterSlaveRelations.size2(), rLhs.size2());
        {
            // Transpose the relation matrix
            typename TSparse::MatrixType transposed_relation_matrix;
            SparseMatrixMultiplicationUtility::TransposeMatrix(transposed_relation_matrix, mpImpl->mMasterSlaveRelations, 1.0);

            // Compute initial independent solution.
            BalancedProduct<TSparse,TSparse,TSparse>(transposed_relation_matrix,
                                                     rSolution,
                                                     mpImpl->mDependentSystem.value().mSolution);
            rSolution.swap(mpImpl->mDependentSystem.value().mSolution);

            // Compute independent RHS.
            BalancedProduct<TSparse,TSparse,TSparse>(transposed_relation_matrix,
                                                     rRhs,
                                                     mpImpl->mDependentSystem.value().mRhs);
            rRhs.swap(mpImpl->mDependentSystem.value().mRhs);

            SparseMatrixMultiplicationUtility::MatrixMultiplication(transposed_relation_matrix, rLhs, left_multiplied_lhs);
        } // deallocate transposed_relation_matrix

        SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs,
                                                                mpImpl->mMasterSlaveRelations,
                                                                mpImpl->mDependentSystem.value().mLhs);
        rLhs.swap(mpImpl->mDependentSystem.value().mLhs);
    } // deallocate left_multiplied_lhs

    // Swap DoFs.
    auto& r_dependent_dofs = mpImpl->mDependentSystem.value().mDofSet;
    r_dependent_dofs.clear();
    typename TSparse::IndexType i_independent_dof = 0;
    for (Dof<typename TDense::DataType>& r_dof : rDofs)
        if (mpImpl->mIndependentMap[r_dof.EquationId()].has_value()) {
            r_dof.SetEquationId(i_independent_dof++);
            r_dependent_dofs.insert(r_dependent_dofs.end(), &r_dof);
        }
    rDofs.swap(r_dependent_dofs);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
bool
MasterSlaveConstraintAssembler<TSparse,TDense>::FinalizeConstraintIteration(typename TSparse::MatrixType& rLhs,
                                                                            typename TSparse::VectorType& rSolution,
                                                                            typename TSparse::VectorType& rRhs,
                                                                            typename Base::DofSet::iterator itDofBegin,
                                                                            typename Base::DofSet::iterator itDofEnd,
                                                                            PMGStatusStream::Report& rReport,
                                                                            PMGStatusStream& rStream)
{
    BalancedProduct<TSparse,TSparse,TSparse>(mpImpl->mMasterSlaveRelations,
                                             rSolution,
                                             mpImpl->mDependentSystem.value().mSolution);
    typename TSparse::VectorType constraint_residual = this->GetConstraintGapVector();
    BalancedProduct<TSparse,TSparse,TSparse>(this->GetRelationMatrix(),
                                             mpImpl->mDependentSystem.value().mSolution,
                                             constraint_residual,
                                             static_cast<typename TSparse::DataType>(-1));

    rReport.maybe_constraint_residual = TSparse::TwoNorm(constraint_residual);
    rReport.constraints_converged = true;
    rStream.Submit(rReport.Tag(2), mpImpl->mVerbosity);
    return true;
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Finalize([[maybe_unused]] typename TSparse::MatrixType& rLhs,
                                                              [[maybe_unused]] typename TSparse::VectorType& rSolution,
                                                              [[maybe_unused]] typename TSparse::VectorType& rRhs,
                                                              [[maybe_unused]] typename Base::DofSet& rDofSet)
{
    // Restore dependent system.
    rLhs.swap(mpImpl->mDependentSystem.value().mLhs);
    rSolution.swap(mpImpl->mDependentSystem.value().mSolution);
    rRhs.swap(mpImpl->mDependentSystem.value().mRhs);

    // Restore dependent DoFs.
    rDofSet.swap(mpImpl->mDependentSystem.value().mDofSet);
    typename TSparse::IndexType i_dependent_dof = 0;
    for (Dof<typename TDense::DataType>& r_dof : rDofSet)
        r_dof.SetEquationId(i_dependent_dof++);

    // Release the independent system.
    mpImpl->mDependentSystem.reset();
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::ComputeIndependentResidual(typename TSparse::VectorType& rResidual) const
{
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::ComputeDependentResidual(typename TSparse::VectorType& rResidual) const
{
    KRATOS_TRY

    typename TSparse::VectorType dependent_residual(mpImpl->mMasterSlaveRelations.size1());
    TSparse::SetToZero(dependent_residual);

    // Compute right hand side.
    typename TSparse::VectorType rhs(dependent_residual.size());
    TSparse::SetToZero(rhs);
    BalancedProduct<TSparse,TSparse,TSparse>(mpImpl->mMasterSlaveRelations,
                                             rResidual,
                                             rhs);

    // Solve linear system.
    auto [r_lhs, r_linear_solver] = mpImpl->GetDependentTransform();
    r_linear_solver.PerformSolutionStep(r_lhs, dependent_residual, rhs);
    rResidual.swap(dependent_residual);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::ComputeIndependentSolution(typename TSparse::VectorType& rSolution) const
{
    KRATOS_TRY
    auto [r_lhs, r_linear_solver] = mpImpl->GetIndependentTransform();

    // Compute reverse RHS.
    // T^T (u - b)
    typename TSparse::VectorType rhs;
    rhs.resize(r_lhs.size1(), false);
    TSparse::SetToZero(rhs);

    typename TSparse::VectorType tmp;
    tmp.resize(rSolution.size());
    TSparse::Copy(rSolution, tmp);

    /// @todo TSparse::UnaliasedAdd(tmp, static_cast<TSparse::DataType>(-1), where are the master-slave constraint gaps?)
    TSparse::TransposeMult(mpImpl->mMasterSlaveRelations, tmp, rhs);

    // Solve the system.
    typename TSparse::VectorType independent_solution(r_lhs.size1());
    TSparse::SetToZero(independent_solution);

    r_linear_solver.PerformSolutionStep(r_lhs, independent_solution, rhs);
    rSolution.swap(independent_solution);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::ComputeDependentSolution(typename TSparse::VectorType& rSolution) const
{
    KRATOS_TRY
    typename TSparse::VectorType dependent_solution(mpImpl->mMasterSlaveRelations.size1());
    TSparse::SetToZero(dependent_solution);
    BalancedProduct<TSparse,TSparse,TSparse>(mpImpl->mMasterSlaveRelations,
                                             rSolution,
                                             dependent_solution);
    /// @todo TSparse::UnaliasedAdd(rOutput, where are the master-slave constraint gaps?);

    rSolution.swap(dependent_solution);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
const typename MasterSlaveConstraintAssembler<TSparse,TDense>::Base::DofSet&
MasterSlaveConstraintAssembler<TSparse,TDense>::GetDependentDofs(const typename Base::DofSet& rIndependentDofSet) const noexcept
{
    return mpImpl->mDependentSystem.value().mDofSet;
}


template <class TSparse, class TDense>
void MasterSlaveConstraintAssembler<TSparse,TDense>::Clear()
{
    Base::Clear();
    mpImpl->mSlaves = decltype(mpImpl->mSlaves)();
    mpImpl->mConstraintIdMap = decltype(mpImpl->mConstraintIdMap)();
    mpImpl->mMasterSlaveRelations = decltype(mpImpl->mMasterSlaveRelations)();
    mpImpl->mDependentSystem.reset();
}


template <class TSparse, class TDense>
Parameters MasterSlaveConstraintAssembler<TSparse,TDense>::GetDefaultParameters()
{
    return Parameters(R"({
        "method" : "master_slave",
        "diagonal_scaling" : "norm",
        "dependent_solver_settings" : {
            "solver_type"   : "amgcl",
            "smoother_type" : "gauss_seidel",
            "krylov_type"   : "gmres"
        },
        "independent_solver_settings" : {
            "solver_type" : "amgcl",
            "smoother_type" : "ilu0",
            "krylov_type" : "cg"
        },
        "verbosity" : 1
    })");
}


template class MasterSlaveConstraintAssembler<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;

template class MasterSlaveConstraintAssembler<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos
