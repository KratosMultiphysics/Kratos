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
#include "solving_strategies/builder_and_solvers/p_multigrid_builder_and_solver.hpp" // PMultigridBuilderAndSolver
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver
#include "includes/model_part.h" // ModelPart
#include "utilities/dof_utilities/block_build_dof_array_utility.h" // BlockBuildDofArrayUtility
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility
#include "utilities/atomic_utilities.h" // AtomicAdd
#include "utilities/reduction_utilities.h"
#include "utilities/proxies.h" // MakeProxy
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE, KRATOS_PROFILE_SCOPE_MILLI
#include "solving_strategies/builder_and_solvers/p_multigrid_utilities.hpp" // MakePRestrictionOperator, MakeSparseTopology

// System includes
#include <algorithm> // std::lower_bound, std::sort, std::transform
#include <optional> // std::optional
#include <limits> // std::numeric_limits
#include <unordered_set> // std::unordered_set
#include <unordered_map> // std::unordered_map
#include <sstream> // std::stringstream
#include <numeric> // std::iota


namespace Kratos {


enum class DiagonalScaling
{
    None        = 0,
    AbsMax      = 1,
    Norm        = 2,
    Constant    = 3
}; // enum class DiagonalScaling


enum class ConstraintImposition
{
    None                = 0,
    MasterSlave         = 1,
    Lagrange            = 2,
    AugmentedLagrange   = 3,
    Penalty             = 4
}; // enum class ConstraintImposition


// --------------------------------------------------------- //
// PMG Hierarchy
// --------------------------------------------------------- //


template <class TSparse,
          class TDense,
          class TSolver>
class PGrid
{
public:
    PGrid(const unsigned Level)
        : mpRestrictionOperator(new typename TSparse::MatrixType),
          mpLhs(new typename TSparse::MatrixType),
          mpRhs(new typename TSparse::VectorType)
    {}


    PGrid()
        : PGrid(0u)
    {}


    template <class TParentSparse>
    void MakeLhsTopology(const ModelPart& rModelPart,
                         const typename TParentSparse::MatrixType& rParentLhs)
    {
        KRATOS_TRY
        // The restriction operator immediately constructs the linear equivalent,
        // because one-level coarsening strategies are not supported yet. The problem
        // is that when going from a generic polynomial level q to some other lower polynomial
        // level p!=1, new nodes must be introduced that do not exist in the original fine mesh.
        // This wouldn't be a huge problem by itself, but deciding on DoF indices on the coarse
        // level would be quite painful and require to keep track of the coarse grid's topology
        // in some manner.
        MakePRestrictionOperator<
            std::numeric_limits<unsigned>::max(),
            typename TSparse::DataType>(rModelPart,
                                        *mpRestrictionOperator,
                                        rParentLhs.size1());
        KRATOS_CATCH("")
    }


    template <bool AssembleLHS,
              bool AssembleRHS,
              class TParentSparse>
    void Assemble(const ModelPart& rModelPart,
                  const typename TParentSparse::MatrixType* pParentLhs,
                  const typename TParentSparse::VectorType* pParentRhs)
    {
        KRATOS_TRY

        // Assemble the coarse LHS matrix if requested.
        if constexpr (AssembleLHS) {
            KRATOS_ERROR_IF_NOT(pParentLhs);

            const auto check_matrix_multiplication = [](const typename TSparse::MatrixType& r_left,
                                                        const typename TSparse::MatrixType& r_right,
                                                        const typename TSparse::MatrixType& r_out) {
                KRATOS_ERROR_IF((r_left.size2() != r_right.size1()) ||
                                (r_out.size1() != r_left.size1()) ||
                                (r_out.size2() != r_right.size2()))
                    << "invalid sizes for matrix product: "
                    << "(" << r_left.size1() << "x" << r_left.size2() << ")"
                    << " x "
                    << "(" << r_right.size1() << "x" << r_right.size2() << ")"
                    << " => "
                    << "(" << r_out.size1() << "x" << r_out.size2() << ")";
            };

            typename TSparse::MatrixType prolongation_operator(mpRestrictionOperator->size2(), mpRestrictionOperator->size1()),
                                         left_multiplied_lhs(mpRestrictionOperator->size1(), pParentLhs->size2());

            check_matrix_multiplication(*mpRestrictionOperator, *pParentLhs, left_multiplied_lhs);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(*mpRestrictionOperator, *pParentLhs, left_multiplied_lhs);
            SparseMatrixMultiplicationUtility::TransposeMatrix(prolongation_operator, *mpRestrictionOperator, 1.0);

            mpLhs->resize(mpRestrictionOperator->size1(), mpRestrictionOperator->size1(), false);
            check_matrix_multiplication(left_multiplied_lhs, prolongation_operator, *mpLhs);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, prolongation_operator, *mpLhs);
        } // if AssembleLHS

        // Assemble the coarse RHS vector if requested.
        if constexpr (AssembleRHS) {
            // Sanity checks
            KRATOS_ERROR_IF_NOT(pParentRhs);
            KRATOS_ERROR_IF_NOT(mpRestrictionOperator->size2() == pParentRhs->size())
                << "expecting an RHS vector of size " << mpRestrictionOperator->size2()
                << " but got " << pParentRhs->size();

            mpRhs->resize(mpRestrictionOperator->size1(), false);
            TSparse::Mult(*mpRestrictionOperator, *pParentRhs, *mpRhs);
        } // if AssembleRHS

        KRATOS_CATCH("")
    }

private:
    std::shared_ptr<typename TSparse::MatrixType> mpRestrictionOperator;

    std::shared_ptr<typename TSparse::MatrixType> mpLhs;

    std::shared_ptr<typename TSparse::VectorType> mpRhs;

    std::optional<std::unique_ptr<PGrid>> mMaybeChild;
}; // class PGrid


// --------------------------------------------------------- //
// PIMPL
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
struct PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Impl
{
    using Interface = PMultigridBuilderAndSolver<TSparse,TDense,TSolver>;

    // --------------------------------------------------------- //
    // Member Variables
    // --------------------------------------------------------- //

    Interface* mpInterface;

    typename Interface::TSystemMatrixType mRelationMatrix;

    typename Interface::TSystemVectorType mConstraintGapVector;

    std::vector<std::size_t> mSlaveIds;

    std::vector<std::size_t> mMasterIds;

    std::unordered_set<std::size_t> mInactiveSlaveIds;

    /// @brief A map associating slave IDs with constraint indices.
    std::unordered_map<std::size_t,std::size_t> mSlaveToConstraintMap;

    /// @brief Penalty factor for augmented lagrangian and penalty constraint imposition.
    typename TSparse::DataType mPenaltyFactor;

    typename TSparse::DataType mDiagonalScaleFactor;

    PGrid<TSparse,TDense,TSolver> mHierarchy;

    DiagonalScaling mDiagonalScaling;

    ConstraintImposition mConstraintImposition;

    int mMaxDepth;

    int mVerbosity;

    // --------------------------------------------------------- //
    // Special Member Functions
    // --------------------------------------------------------- //

    Impl(Interface* pInterface)
        : mpInterface(pInterface),
          mRelationMatrix(),
          mConstraintGapVector(),
          mSlaveIds(),
          mMasterIds(),
          mInactiveSlaveIds(),
          mPenaltyFactor(1e3),
          mDiagonalScaleFactor(1),
          mHierarchy(),
          mDiagonalScaling(DiagonalScaling::None),
          mConstraintImposition(ConstraintImposition::MasterSlave),
          mMaxDepth(-1),
          mVerbosity(0)
    {}

    // --------------------------------------------------------- //
    // Solution
    // --------------------------------------------------------- //

    /// @brief Initialize the linear solver and solve the provided system.
    void Solve(typename Interface::TSystemMatrixType& rLhs,
               typename Interface::TSystemVectorType& rSolution,
               typename Interface::TSystemVectorType& rRhs,
               ModelPart& rModelPart,
               Interface& rInterface)
    {
        KRATOS_TRY
        if (rInterface.GetLinearSolver().AdditionalPhysicalDataIsNeeded()) {
            rInterface.GetLinearSolver().ProvideAdditionalData(rLhs,
                                                               rSolution,
                                                               rRhs,
                                                               rInterface.mDofSet,
                                                               rModelPart);
        }
        KRATOS_CATCH("")

        if(!rModelPart.MasterSlaveConstraints().empty()) {
            switch (this->mConstraintImposition) {
                case ConstraintImposition::MasterSlave: {
                    KRATOS_TRY
                    typename Interface::TSystemVectorType modified_solution(rRhs.size());
                    TSparse::SetToZero(modified_solution);
                    rInterface.GetLinearSolver().Solve(rLhs, modified_solution, rRhs);

                    // Recover the solution of the original problem.
                    TSparse::Mult(mRelationMatrix, modified_solution, rSolution);
                    KRATOS_CATCH("")
                    break;
                } // case ConstraintImposition::MasterSlave

                case ConstraintImposition::AugmentedLagrange: {
                    KRATOS_TRY
                    typename Interface::TSystemMatrixType transpose_relation_matrix;
                    SparseMatrixMultiplicationUtility::TransposeMatrix(transpose_relation_matrix, mRelationMatrix);

                    for (std::size_t i_iteration=0ul; i_iteration<10; ++i_iteration) {
                        // Solve the current system.
                        rInterface.GetLinearSolver().Solve(rLhs, rSolution, rRhs);

                        // Compute the constraint residuals.
                        typename Interface::TSystemVectorType constraint_residual(mConstraintGapVector.size());
                        TSparse::SetToZero(constraint_residual);

                        KRATOS_TRY
                        TSparse::Mult(mRelationMatrix, rSolution, constraint_residual);
                        TSparse::UnaliasedAdd(constraint_residual, -1.0, mConstraintGapVector);
                        KRATOS_CATCH("")

                        // Decide whether to update the RHS or terminate the loop.
                        const auto constraint_norm = TSparse::TwoNorm(constraint_residual);
                        if (3 <= this->mVerbosity) {
                            std::stringstream residual_stream;
                            residual_stream << std::setprecision(8) << std::scientific<< constraint_norm;
                            std::cout << mpInterface->Info() << ": iteration " << i_iteration << " constraint residual " << residual_stream.str() << "\n";
                        }

                        if (1e-6 < constraint_norm) {
                            KRATOS_TRY
                            typename Interface::TSystemVectorType rhs_update(rRhs.size());
                            TSparse::SetToZero(rhs_update);
                            TSparse::InplaceMult(constraint_residual, -mPenaltyFactor);
                            TSparse::Mult(transpose_relation_matrix, constraint_residual, rhs_update);
                            TSparse::UnaliasedAdd(rRhs, 1.0, rhs_update);
                            KRATOS_CATCH("")
                        } else {
                            break;
                        }
                    } // for i_iteration
                    KRATOS_CATCH("")
                    break;
                } // case ConstraintImposition::AugmentedLagrange

                default: {
                    KRATOS_ERROR << "Unsupported constraint imposition method " << (int)this->mConstraintImposition;
                }
            } // switch constraint imposition
        } /*if rModelPart.MasterSlaveConstraints()*/ else {
            KRATOS_TRY
            rInterface.GetLinearSolver().Solve(rLhs, rSolution, rRhs);
            KRATOS_CATCH("")
        } // if !rModelPart.MasterSlaveConstraints()
    }


    // --------------------------------------------------------- //
    // Topology
    // --------------------------------------------------------- //


    template <class TProxy>
    static void CollectDoFs(const TProxy& rEntities,
                            const ProcessInfo& rProcessInfo,
                            typename Interface::TSchemeType& rScheme,
                            LockObject* pLockBegin,
                            [[maybe_unused]] LockObject* pLockEnd,
                            std::unordered_set<std::size_t>* pRowSetBegin,
                            [[maybe_unused]] std::unordered_set<std::size_t>* pRowSetEnd)
    {
        KRATOS_TRY
        KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
        using TLS = Element::EquationIdVectorType;
        block_for_each(rEntities,
                       TLS(),
                       [&rScheme, &rProcessInfo, pLockBegin, pRowSetBegin](const auto& rEntity, TLS& rTls){
            if (rEntity.GetEntity().IsActive()) {
                rScheme.EquationId(rEntity.GetEntity(), rTls, rProcessInfo);
                for (const auto equation_id : rTls) {
                    [[maybe_unused]] std::scoped_lock<LockObject> lock(pLockBegin[equation_id]);
                    auto& r_row_indices = pRowSetBegin[equation_id];
                    r_row_indices.insert(rTls.begin(), rTls.end());
                } // for equation_id in rTls
            } // if rEntity.IsActive
        });
        KRATOS_CATCH("")
    }


    void MakeLhsTopology(const typename Interface::TSchemeType::Pointer& rpScheme,
                         typename Interface::TSystemMatrixType& rLhs,
                         ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        std::vector<LockObject> mutexes(mpInterface->GetEquationSystemSize());
        std::vector<std::unordered_set<std::size_t> > indices(mpInterface->GetEquationSystemSize());

        // Collect DoFs from elements.
        Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Element>(rModelPart),
                          rModelPart.GetProcessInfo(),
                          *rpScheme,
                          mutexes.data(),
                          mutexes.data() + mutexes.size(),
                          indices.data(),
                          indices.data() + indices.size());

        // Collect DoFs from conditions.
        Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Condition>(rModelPart),
                          rModelPart.GetProcessInfo(),
                          *rpScheme,
                          mutexes.data(),
                          mutexes.data() + mutexes.size(),
                          indices.data(),
                          indices.data() + indices.size());

        // Collect DoFs from multifreedom constraints.
        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            switch (this->mConstraintImposition) {
                case ConstraintImposition::MasterSlave: {
                    KRATOS_TRY
                    struct TLS {
                        Element::EquationIdVectorType master_ids, slave_ids;
                    };

                    block_for_each(rModelPart.MasterSlaveConstraints(), TLS(), [&](MasterSlaveConstraint& r_constraint, TLS& r_tls){
                        r_constraint.EquationIdVector(r_tls.slave_ids, r_tls.master_ids, r_current_process_info);

                        for (std::size_t i = 0; i < r_tls.slave_ids.size(); i++) {
                            auto& row_indices = indices[r_tls.slave_ids[i]];
                            std::scoped_lock<LockObject> lock(mutexes[r_tls.slave_ids[i]]);
                            row_indices.insert(r_tls.slave_ids[i]);
                        }

                        for (std::size_t i = 0; i < r_tls.master_ids.size(); i++) {
                            auto& row_indices = indices[r_tls.master_ids[i]];
                            std::scoped_lock<LockObject> lock(mutexes[r_tls.master_ids[i]]);
                            row_indices.insert(r_tls.master_ids[i]);
                        }
                    });
                    KRATOS_CATCH("")
                    break;
                } // case ConstraintImposition::MasterSlave

                case ConstraintImposition::AugmentedLagrange: {
                    KRATOS_TRY
                    struct TLS {
                        Element::EquationIdVectorType master_ids, slave_ids;
                    };

                    block_for_each(rModelPart.MasterSlaveConstraints(), TLS(), [&](MasterSlaveConstraint& r_constraint, TLS& r_tls){
                        r_constraint.EquationIdVector(r_tls.slave_ids, r_tls.master_ids, r_current_process_info);
                        for (const auto i_row : r_tls.slave_ids) {
                            for (const auto i_column : r_tls.master_ids) {
                                std::scoped_lock<LockObject> lock(mutexes[i_column]);
                                indices[i_column].insert(r_tls.master_ids.begin(), r_tls.master_ids.end());
                                indices[i_column].insert(i_row);
                            } // for i_column in r_tls.master_ids
                            std::scoped_lock<LockObject> lock(mutexes[i_row]);
                            indices[i_row].insert(r_tls.master_ids.begin(), r_tls.master_ids.end());
                            indices[i_row].insert(i_row);
                        } // for i_row in r_tls.master_ids
                    });
                    KRATOS_CATCH("")
                    break;
                } // ConstraintImposition::AugmentedLagrange

                default: {
                    KRATOS_ERROR << "Unsupported constraint imposition method" << (int)this->mConstraintImposition;
                }
            } // swith constraint imposition
        } // if rModelPart.MasterSlaveConstraints

        // Destroy mutexes.
        mutexes = std::vector<LockObject>();

        // Compute and allocate LHS topology.
        MakeSparseTopology<typename TSparse::DataType>(indices,
                                                       indices.size(),
                                                       rLhs);

        // Construct the coarse hierarhy's topology.
        mHierarchy.template MakeLhsTopology<TSparse>(rModelPart, rLhs);
        KRATOS_CATCH("")
    }


    /// @brief Compute the sparsity pattern of the multifreedom constraints' relation matrix assuming imposition via master-slave elimination.
    void MakeMasterSlaveConstraintTopology(const ModelPart& rModelPart)
    {
        KRATOS_TRY
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        std::vector<std::unordered_set<IndexType>> indices(mpInterface->mDofSet.size());
        std::vector<LockObject> mutexes(indices.size());
        const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();

        #pragma omp parallel
        {
            Element::EquationIdVectorType slave_ids;
            Element::EquationIdVectorType master_ids;
            std::unordered_map<IndexType, std::unordered_set<IndexType>> temp_indices;

            #pragma omp for nowait
            for (int i_const = 0; i_const < static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i_const) {
                auto it_const = it_const_begin + i_const;
                it_const->EquationIdVector(slave_ids, master_ids, r_process_info);

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

        MakeSparseTopology<typename TSparse::DataType>(indices,
                                                       indices.size(),
                                                       mRelationMatrix);
        mConstraintGapVector.resize(indices.size(), false);

        KRATOS_CATCH("")
    }


    /// @brief Compute the sparsity pattern of the multifreedom constraints' relation matrix assuming imposition via augmented lagrange.
    void MakeAugmentedLagrangeConstraintTopology(const ModelPart& rModelPart)
    {
        KRATOS_TRY

        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Construct a map that associates slave indices with constraint equation indices.
        mSlaveToConstraintMap.clear();

        {
            MasterSlaveConstraint::IndexType i_constraint = 0;
            MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
            for (const auto& r_constraint : rModelPart.MasterSlaveConstraints()) {
                r_constraint.EquationIdVector(slave_ids, master_ids, r_process_info);
                for (const auto i_slave : slave_ids) {
                    const auto emplace_result = mSlaveToConstraintMap.emplace(i_slave, i_constraint);
                    if (emplace_result.second) ++i_constraint;
                } // for i_slave in slave_ids
            } // for r_constraint in rModelPart.MasterSlaveConstraints
        }

        std::vector<std::unordered_set<IndexType>> indices(mSlaveToConstraintMap.size());
        std::vector<LockObject> mutexes(indices.size());

        struct TLS {
            Element::EquationIdVectorType master_ids, slave_ids;
        }; // struct TLS

        block_for_each(rModelPart.MasterSlaveConstraints(),
                       TLS(),
                       [&mutexes, &indices, &r_process_info, this](const auto& r_constraint, TLS& r_tls){
            r_constraint.EquationIdVector(r_tls.slave_ids,
                                          r_tls.master_ids,
                                          r_process_info);

            for (const auto i_slave : r_tls.slave_ids) {
                const auto i_constraint = mSlaveToConstraintMap[i_slave];
                std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                indices[i_constraint].insert(r_tls.master_ids.begin(), r_tls.master_ids.end());
                indices[i_constraint].insert(i_slave);
            } // for i_slave in slave_ids
        }); // for r_constraint in rModelPart.MasterSlaveConstraints()

        MakeSparseTopology<typename TSparse::DataType>(indices,
                                                       mpInterface->mDofSet.size(),
                                                       mRelationMatrix);
        mConstraintGapVector.resize(indices.size(), false);

        KRATOS_CATCH("")
    }


    /// @brief Compute the sparsity pattern of the multifreedom constraints' relation matrix for the requested imposition method.
    void MakeConstraintTopology(const ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY
        if (!rModelPart.MasterSlaveConstraints().empty()) {
            switch (this->mConstraintImposition) {
                case ConstraintImposition::MasterSlave: {
                    this->MakeMasterSlaveConstraintTopology(rModelPart);
                    break;
                } // case ConstraintImposition::MasterSlave
                case ConstraintImposition::AugmentedLagrange: {
                    this->MakeAugmentedLagrangeConstraintTopology(rModelPart);
                    break;
                } // case ConstraintImposition::AugmentedLagrange
                default: {
                    KRATOS_ERROR << "Unsupported constraint imposition method: " << (int)this->mConstraintImposition;
                }
            } // switch constraint imposition
        }// if rModelPart.MasterSlaveConstraints()
        KRATOS_CATCH("")
    }


    // --------------------------------------------------------- //
    // Mapping
    // --------------------------------------------------------- //


    /// @brief Map LHS contributions from a local dense matrix to a sparse global matrix.
    static void MapRowContribution(typename Interface::TSystemMatrixType& rLhs,
                                   const typename Interface::LocalSystemMatrixType& rLocalLhs,
                                   const unsigned iRow,
                                   const unsigned iLocalRow,
                                   const Element::EquationIdVectorType& rEquationIds) noexcept
    {
        auto& r_entries = rLhs.value_data();
        const auto& r_row_extents = rLhs.index1_data();
        const auto& r_column_indices = rLhs.index2_data();

        const auto i_row_begin = r_row_extents[iRow];
        const auto i_row_end = r_row_extents[iRow + 1];

        for (unsigned int i_local_column=0; i_local_column<rEquationIds.size(); i_local_column++) {
            const unsigned int iColumn = rEquationIds[i_local_column];
            const auto it = std::lower_bound(r_column_indices.begin() + i_row_begin,
                                             r_column_indices.begin() + i_row_end,
                                             iColumn);
            r_entries[std::distance(r_column_indices.begin(), it)] += rLocalLhs(iLocalRow, i_local_column);
        }
    }


    /// @brief Map RHS contributions from local to global space.
    static void MapContribution(typename Interface::TSystemVectorType& rRhs,
                                const typename Interface::LocalSystemVectorType& rContribution,
                                const Element::EquationIdVectorType& rEquationIds) noexcept
    {
        const unsigned local_size = rContribution.size();

        for (unsigned i_local = 0; i_local < local_size; i_local++) {
            const unsigned i_global = rEquationIds[i_local];
            AtomicAdd(rRhs[i_global], rContribution[i_local]);
        }
    }


    /// @brief Map LHS contributions from local to global space.
    static void MapContribution(typename Interface::TSystemMatrixType& rLhs,
                                const typename Interface::LocalSystemMatrixType& rContribution,
                                const Element::EquationIdVectorType& rEquationIds,
                                LockObject* pLockBegin) noexcept
    {
        const std::size_t local_size = rContribution.size1();
        for (IndexType i_local = 0; i_local < local_size; i_local++) {
            const IndexType i_global = rEquationIds[i_local];
            std::scoped_lock<LockObject> lock(pLockBegin[i_global]);
            Impl::MapRowContribution(rLhs, rContribution, i_global, i_local, rEquationIds);
        }
    }


    /// @brief Map LHS and RHS contributions from local to global space.
    static void MapContribution(typename Interface::TSystemMatrixType& rLhs,
                                typename Interface::TSystemVectorType& rRhs,
                                const typename Interface::LocalSystemMatrixType& rLhsContribution,
                                const typename Interface::LocalSystemVectorType& rRhsContribution,
                                const Element::EquationIdVectorType& rEquationIds,
                                LockObject* pLockBegin) noexcept
    {
        const unsigned local_size = rLhsContribution.size1();
        for (unsigned i_local = 0; i_local < local_size; i_local++) {
            const unsigned i_global = rEquationIds[i_local];
            std::scoped_lock<LockObject> lock(pLockBegin[i_global]);
            rRhs[i_global] += rRhsContribution[i_local];
            Impl::MapRowContribution(rLhs,
                                     rLhsContribution,
                                     i_global,
                                     i_local,
                                     rEquationIds);
        }
    }


    /// @brief Compute the local contributions of an @ref Element or @ref Condition and assemble them into the global system.
    /// @details Whether from elements or conditions, whether assembling into the LHS, RHS, or both,
    ///          all cases are handled in this one function to reduce code duplication and have logically
    ///          coherent parts of the code in one place.
    template <bool AssembleLHS,
              bool AssembleRHS,
              class TEntity>
    static void MapEntityContribution(TEntity& rEntity,
                                      typename Interface::TSchemeType& rScheme,
                                      const ProcessInfo& rProcessInfo,
                                      Element::EquationIdVectorType& rEquationIndices,
                                      typename Interface::LocalSystemMatrixType* pLhsContribution,
                                      typename Interface::LocalSystemVectorType* pRhsContribution,
                                      typename Interface::TSystemMatrixType* pLhs,
                                      typename Interface::TSystemVectorType* pRhs,
                                      LockObject* pLockBegin)
    {
        if constexpr (AssembleLHS || AssembleRHS) {
            if (rEntity.IsActive()) {
                if constexpr (AssembleLHS) {
                    if constexpr (AssembleRHS) {
                        // Assemble LHS and RHS.
                        rScheme.CalculateSystemContributions(rEntity,
                                                             *pLhsContribution,
                                                             *pRhsContribution,
                                                             rEquationIndices,
                                                             rProcessInfo);
                        Impl::MapContribution(*pLhs,
                                              *pRhs,
                                              *pLhsContribution,
                                              *pRhsContribution,
                                              rEquationIndices,
                                              pLockBegin);
                    } /*if AssembleRHS*/ else {
                        // Assemble LHS only.
                        rScheme.CalculateLHSContribution(rEntity,
                                                        *pLhsContribution,
                                                        rEquationIndices,
                                                        rProcessInfo);
                        Impl::MapContribution(*pLhs,
                                              *pLhsContribution,
                                              rEquationIndices,
                                              pLockBegin);
                    } // if !AssembleRHS
                } /*if AssembleLHS*/ else {
                    // Assemble RHS only.
                    rScheme.CalculateRHSContribution(rEntity,
                                                     *pRhsContribution,
                                                     rEquationIndices,
                                                     rProcessInfo);
                    Impl::MapContribution(*pRhs,
                                          *pRhsContribution,
                                          rEquationIndices);
                } // if !AssembleLHS
            } // if rEntity.IsActive
        } // if AssembleLHS or AssembleRHS
    }


    // --------------------------------------------------------- //
    // Assembly
    // --------------------------------------------------------- //


    /// @brief Compute and assemble local contributions from elements and conditions into the global system.
    /// @details This function body mainly handles looping over elements/conditions and its parallelization.
    ///          The actual local system calculation, as well as assembly, is deferred to @ref MapEntityContribution.
    template <bool AssembleLHS,
              bool AssembleRHS>
    void Assemble(ModelPart& rModelPart,
                  typename Interface::TSchemeType& rScheme,
                  std::optional<typename Interface::TSystemMatrixType*> pMaybeLhs,
                  std::optional<typename Interface::TSystemVectorType*> pMaybeRhs)
    {
        KRATOS_TRY

        // Sanity checks.
        if constexpr (AssembleLHS) {
            KRATOS_ERROR_IF_NOT(pMaybeLhs.has_value() && pMaybeLhs.value() != nullptr)
                << "Requested an assembly of the left hand side, but no matrix is provided to assemble into.";
        }

        if constexpr (AssembleRHS) {
            KRATOS_ERROR_IF_NOT(pMaybeRhs.has_value() && pMaybeRhs.value() != nullptr)
                << "Requested an assembly of the right hand side, but no vector is provided to assemble into.";
        }

        // Global variables.
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        typename Interface::TSystemMatrixType* pLhs = pMaybeLhs.has_value() ? pMaybeLhs.value() : nullptr;
        typename Interface::TSystemVectorType* pRhs = pMaybeRhs.has_value() ? pMaybeRhs.value() : nullptr;
        auto p_locks = std::make_unique<std::vector<LockObject>>(AssembleLHS ? pLhs->size1() : 0ul);

        [[maybe_unused]] const int element_count = rModelPart.Elements().size();
        [[maybe_unused]] const int condition_count = rModelPart.Conditions().size();
        [[maybe_unused]] const auto it_element_begin = rModelPart.Elements().begin();
        [[maybe_unused]] const auto it_condition_begin = rModelPart.Conditions().begin();

        // Collect contributions from constitutive entities.
        #pragma omp parallel
        {
            // Thread-local variables.
            Element::EquationIdVectorType equation_indices;
            typename Interface::LocalSystemMatrixType lhs_contribution;
            typename Interface::LocalSystemVectorType rhs_contribution;

            // Collect contributions from elements.
            #pragma omp for schedule(guided, 512) nowait
            for (int i_entity=0; i_entity<element_count; ++i_entity) {
                Impl::MapEntityContribution<AssembleLHS,AssembleRHS>(*(it_element_begin + i_entity),
                                                                     rScheme,
                                                                     r_process_info,
                                                                     equation_indices,
                                                                     &lhs_contribution,
                                                                     &rhs_contribution,
                                                                     pLhs,
                                                                     pRhs,
                                                                     p_locks->data());
            } // pragma omp for

            // Collect contributions from conditions.
            #pragma omp for schedule(guided, 512)
            for (int i_entity=0; i_entity<condition_count; ++i_entity) {
                Impl::MapEntityContribution<AssembleLHS,AssembleRHS>(*(it_condition_begin + i_entity),
                                                                     rScheme,
                                                                     r_process_info,
                                                                     equation_indices,
                                                                     &lhs_contribution,
                                                                     &rhs_contribution,
                                                                     pLhs,
                                                                     pRhs,
                                                                     p_locks->data());
            } // pragma omp for
        } // pragma omp parallel

        // Assemble coarse hierarchy.
        mHierarchy.template Assemble<AssembleLHS, AssembleRHS,TSparse>(rModelPart,
                                                                       pLhs,
                                                                       pRhs);

        KRATOS_CATCH("")
    }


    /// @brief Compute the relation matrix and constraint gaps for multifreedom constraints assuming imposition via master-slave elimination.
    void MakeMasterSlaveConstraints(const ModelPart& rModelPart)
    {
        KRATOS_TRY

        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        std::vector<LockObject> mutexes(mpInterface->mDofSet.size());

        // Declare local containers.
        typename Interface::LocalSystemMatrixType local_relation_matrix;
        typename Interface::LocalSystemVectorType local_constraint_gap_vector;
        Element::EquationIdVectorType slave_equation_ids, master_equation_ids;
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        #pragma omp parallel firstprivate(local_relation_matrix, local_constraint_gap_vector, slave_equation_ids, master_equation_ids)
        {
            std::unordered_set<IndexType> tls_inactive_slave_dofs;

            #pragma omp for
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

                // If the constraint is active
                if (it_const->IsActive()) {
                    it_const->EquationIdVector(slave_equation_ids,
                                               master_equation_ids,
                                               r_process_info);
                    it_const->CalculateLocalSystem(local_relation_matrix,
                                                   local_constraint_gap_vector,
                                                   r_process_info);

                    for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                        const IndexType i_global = slave_equation_ids[i];

                        // Assemble matrix row.
                        {
                            std::scoped_lock<LockObject> lock(mutexes[i_global]);
                            Impl::MapRowContribution(mRelationMatrix,
                                                     local_relation_matrix,
                                                     i_global,
                                                     i,
                                                     master_equation_ids);
                        }

                        // Assemble constant vector
                        KRATOS_ERROR_IF_NOT(i_global < mConstraintGapVector.size())
                            << "constraint gap vector size " << mConstraintGapVector.size()
                            << " <= insertion index " << i_global << "\n";
                        const double constant_value = local_constraint_gap_vector[i];
                        AtomicAdd(mConstraintGapVector[i_global], constant_value);
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
            mConstraintGapVector[eq_id] = 0.0;
            mRelationMatrix(eq_id, eq_id) = 1.0;
        }

        // Setting inactive slave dofs in the T and C system
        for (auto eq_id : mInactiveSlaveIds) {
            std::scoped_lock<LockObject> lock(mutexes[eq_id]);
            mConstraintGapVector[eq_id] = 0.0;
            mRelationMatrix(eq_id, eq_id) = 1.0;
        }

        KRATOS_CATCH("")
    }


    /// @brief Compute the relation matrix and constraint gaps for multifreedom constraints assuming imposition via augmented lagrange.
    void MakeAugmentedLagrangeConstraints(const ModelPart& rModelPart)
    {
        KRATOS_TRY

        // Function-wide variables.
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        std::vector<LockObject> mutexes(mSlaveToConstraintMap.size());

        // Thread-local storage for constraint assembly.
        struct TLS {
            MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
            typename Interface::LocalSystemMatrixType local_relation_matrix;
            typename Interface::LocalSystemVectorType local_constraint_gap_vector;
        }; // struct TLS

        // Constraint assembly.
        block_for_each(rModelPart.MasterSlaveConstraints(),
                        TLS(),
                        [&mutexes, &r_process_info, this](const MasterSlaveConstraint& r_constraint, TLS& r_tls){
            if (r_constraint.IsActive()) {
                r_constraint.EquationIdVector(r_tls.slave_ids,
                                              r_tls.master_ids,
                                              r_process_info);
                r_constraint.CalculateLocalSystem(r_tls.local_relation_matrix,
                                                  r_tls.local_constraint_gap_vector,
                                                  r_process_info);

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
                    KRATOS_ERROR_IF(mRelationMatrix.size1() <= i_constraint);

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
                        Impl::MapRowContribution(mRelationMatrix,
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


    /// @brief Compute the relation matrix and constraint gap vector for multifreedom constraints using the requested imposition method.
    void MakeConstraints(const ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);

        // Sanity checks.
        KRATOS_ERROR_IF_NOT(mRelationMatrix.size1() == mConstraintGapVector.size());

        KRATOS_TRY

        // Init.
        mInactiveSlaveIds.clear();
        TSparse::SetToZero(mRelationMatrix);
        TSparse::SetToZero(mConstraintGapVector);

        switch (this->mConstraintImposition) {
            case ConstraintImposition::MasterSlave: {
                this->MakeMasterSlaveConstraints(rModelPart);
                break;
            } // case ConstraintImposition::MasterSlave
            case ConstraintImposition::AugmentedLagrange: {
                this->MakeAugmentedLagrangeConstraints(rModelPart);
                break;
            } // case ConstraintImposition::AugmentedLagrange
            default: {
                KRATOS_ERROR << "Unsupported constraint imposition method: " << (int)this->mConstraintImposition;
            }
        } // switch constraint imposition

        KRATOS_CATCH("")
    }


    static void GetDiagonalScaleFactor(typename TSparse::DataType& rDiagonalScaleFactor,
                                       const typename TSparse::MatrixType& rMatrix,
                                       const DiagonalScaling ScalingStrategy,
                                       [[maybe_unused]] const ProcessInfo& rProcessInfo)
    {
        KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
        KRATOS_TRY
        switch (ScalingStrategy) {
            case DiagonalScaling::None:
                rDiagonalScaleFactor = 1;
                break;

            case DiagonalScaling::AbsMax:
                rDiagonalScaleFactor = IndexPartition(rMatrix.size1()).template for_each<MaxReduction<typename TSparse::DataType>>(
                    [&rMatrix](std::size_t iRow) -> typename TSparse::DataType {
                        const auto itBegin = rMatrix.index2_data().begin() + rMatrix.index1_data()[iRow];
                        const auto itEnd = rMatrix.index2_data().begin() + rMatrix.index1_data()[iRow + 1];

                        // Look for the diagonal entry in the current row.
                        const auto itColumnIndex = std::lower_bound(itBegin, itEnd, iRow);
                        KRATOS_ERROR_IF(itBegin == itEnd || *itColumnIndex != iRow)
                            << "row " << iRow << " has no diagonal entry";

                        const auto diagonal_entry = rMatrix.value_data()[std::distance(itBegin, itColumnIndex)];
                        return std::abs(diagonal_entry);
                });
                break;

            case DiagonalScaling::Norm:
                rDiagonalScaleFactor = IndexPartition(rMatrix.size1()).template for_each<AbsMaxReduction<typename TSparse::DataType>>(
                    [&rMatrix](std::size_t iRow) -> typename TSparse::DataType {
                        const auto itBegin = rMatrix.index2_data().begin() + rMatrix.index1_data()[iRow];
                        const auto itEnd = rMatrix.index2_data().begin() + rMatrix.index1_data()[iRow + 1];

                        // Look for the diagonal entry in the current row.
                        const auto itColumnIndex = std::lower_bound(itBegin, itEnd, iRow);
                        KRATOS_ERROR_IF(itBegin == itEnd || *itColumnIndex != iRow)
                            << "row " << iRow << " has no diagonal entry";

                        const auto diagonal_entry = rMatrix.value_data()[std::distance(itBegin, itColumnIndex)];
                        return diagonal_entry * diagonal_entry;
                });
                rDiagonalScaleFactor = std::sqrt(rDiagonalScaleFactor);
                break;

            case DiagonalScaling::Constant:
                rDiagonalScaleFactor = rProcessInfo.GetValue(BUILD_SCALE_FACTOR);
                break;

            default: {
                KRATOS_ERROR << "unsupported diagonal scaling (" << (int)ScalingStrategy << ')';
            }
        } // switch ScalingStrategy
        KRATOS_CATCH("")
    }


    // --------------------------------------------------------- //
    // Constraint Imposition
    // --------------------------------------------------------- //


    void ImposeMasterSlaveConstraints([[maybe_unused]] const ModelPart& rModelPart,
                                      typename Interface::TSystemMatrixType& rLhs,
                                      typename Interface::TSystemVectorType& rRhs)
    {
        KRATOS_TRY
        // Compute the transposed matrix of the global relation matrix
        {
            // Storage for an intermediate matrix transpose(relation_matrix) * lhs_matrix
            typename Interface::TSystemMatrixType left_multiplied_lhs(this->mRelationMatrix.size2(), rLhs.size2());
            {
                // Transpose the relation matrix
                typename Interface::TSystemMatrixType transposed_relation_matrix(this->mRelationMatrix.size2(), this->mRelationMatrix.size1());
                SparseMatrixMultiplicationUtility::TransposeMatrix(transposed_relation_matrix, this->mRelationMatrix, 1.0);

                typename Interface::TSystemVectorType b_modified(rRhs.size());
                TSparse::Mult(transposed_relation_matrix, rRhs, b_modified);
                TSparse::Copy(b_modified, rRhs);

                SparseMatrixMultiplicationUtility::MatrixMultiplication(transposed_relation_matrix, rLhs, left_multiplied_lhs);
            } // deallocate transposed_relation_matrix

            SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, this->mRelationMatrix, rLhs);
        } // deallocate left_multiplied_lhs

        // Compute the scale factor for slave DoFs.
        Impl::GetDiagonalScaleFactor(this->mDiagonalScaleFactor,
                                    rLhs,
                                    this->mDiagonalScaling,
                                    rModelPart.GetProcessInfo());

        // Apply diagonal values on slaves
        block_for_each(this->mSlaveIds, [this, &rLhs, &rRhs](const auto iSlave){
            if (this->mInactiveSlaveIds.find(iSlave) == this->mInactiveSlaveIds.end()) {
                rLhs(iSlave, iSlave) = this->mDiagonalScaleFactor;
                rRhs[iSlave] = 0.0;
            }
        });
        KRATOS_CATCH("")
    }


    void ImposeAugmentedLagrangeConstraints([[maybe_unused]] const ModelPart& rModelPart,
                                            typename Interface::TSystemMatrixType& rLhs,
                                            typename Interface::TSystemVectorType& rRhs)
    {
        KRATOS_TRY

        using SparseUtils = SparseMatrixMultiplicationUtility;

        {
            typename Interface::TSystemMatrixType transposed_relation_matrix;
            SparseUtils::TransposeMatrix(transposed_relation_matrix, mRelationMatrix);

            typename Interface::TSystemMatrixType relation_product;
            SparseUtils::MatrixMultiplication(transposed_relation_matrix,
                                              mRelationMatrix,
                                              relation_product);

            // Add terms to the LHS matrix.
            SparseUtils::InPlaceMatrixAdd(rLhs,
                                          relation_product,
                                          mPenaltyFactor);

            // Add terms to the RHS vector.
            typename Interface::TSystemVectorType rhs_term(rRhs.size()),
                                                  lagrange_multipliers(mRelationMatrix.size1(), 0),
                                                  constraint_gaps = mConstraintGapVector;
            KRATOS_ERROR_IF_NOT(constraint_gaps.size() == lagrange_multipliers.size());

            TSparse::UnaliasedAdd(lagrange_multipliers,
                                  -mPenaltyFactor,
                                  constraint_gaps);
            TSparse::Mult(transposed_relation_matrix,
                          lagrange_multipliers,
                          rhs_term);
            TSparse::UnaliasedAdd(rRhs,
                                  -1.0,
                                  rhs_term);
        }

        KRATOS_CATCH("")
    }


    void ImposeConstraints([[maybe_unused]] typename Interface::TSchemeType& rScheme,
                           const ModelPart& rModelPart,
                           typename Interface::TSystemMatrixType& rLhs,
                           typename Interface::TSystemVectorType& rRhs)
    {
        KRATOS_TRY
        if (!rModelPart.MasterSlaveConstraints().empty()) {
            this->MakeConstraints(rModelPart);

            switch (this->mConstraintImposition) {
                case ConstraintImposition::MasterSlave: {
                this->ImposeMasterSlaveConstraints(rModelPart,
                                                   rLhs,
                                                   rRhs);
                break;
            } // case ConstraintImposition::MasterSlave
            case ConstraintImposition::AugmentedLagrange: {
                this->ImposeAugmentedLagrangeConstraints(rModelPart,
                                                         rLhs,
                                                         rRhs);
                break;
            } // case ConstraintImposition::AugmentedLagrange
            default: {
                KRATOS_ERROR << "Unsupported constraint imposition method: " << (int)this->mConstraintImposition;
            }
        } // switch constraint imposition
    } // if !rModelPart.MasterSlaveConstraints().empty()
        KRATOS_CATCH("")
    }

}; // class PMultigridBuilderAndSolver::Impl


// --------------------------------------------------------- //
// Constructors
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::~PMultigridBuilderAndSolver() = default;


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::PMultigridBuilderAndSolver()
    : Interface(),
      mpImpl(new Impl(this))
{
}


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::PMultigridBuilderAndSolver(const typename TSolver::Pointer& pSolver,
                                                                               Parameters Settings)
    : Interface(pSolver),
      mpImpl(new Impl(this))
{
    Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());
    this->AssignSettings(Settings);
}


template <class TSparse, class TDense, class TSolver>
typename PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Interface::Pointer
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Create(typename TSolver::Pointer pSolver,
                                                           Parameters Settings) const
{
    KRATOS_TRY
    return typename Interface::Pointer(new PMultigridBuilderAndSolver(pSolver, Settings));
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Allocation and Initialization
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::SetUpDofSet(typename Interface::TSchemeType::Pointer pScheme,
                                                                     ModelPart& rModelPart)
{
    KRATOS_TRY;

    BlockBuildDofArrayUtility::SetUpDofArray(rModelPart,
                                             this->mDofSet,
                                             this->GetEchoLevel(),
                                             this->GetCalculateReactionsFlag());

    KRATOS_CATCH("");
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::SetUpSystem(ModelPart& rModelPart)
{
    this->mEquationSystemSize = this->mDofSet.size();

    KRATOS_TRY
    // Set equation indices of DoFs.
    IndexPartition<std::size_t>(this->mDofSet.size()).for_each([&, this](std::size_t Index){
        (this->mDofSet.begin() + Index)->SetEquationId(Index);
    });
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ResizeAndInitializeVectors(typename Interface::TSchemeType::Pointer pScheme,
                                                                                    typename Interface::TSystemMatrixPointerType& rpLhs,
                                                                                    typename Interface::TSystemVectorPointerType& rpSolution,
                                                                                    typename Interface::TSystemVectorPointerType& rpRhs,
                                                                                    ModelPart& rModelPart)
{
    KRATOS_TRY

    // Construct empty containers if necessary.
    if (!rpLhs)
        rpLhs.reset(new typename Interface::TSystemMatrixType);

    if (!rpSolution)
        rpSolution.reset(new typename Interface::TSystemVectorType);

    if (!rpRhs)
        rpRhs.reset(new typename Interface::TSystemVectorType);

    // Construct LHS topology if necessary or requested.
    if (rpLhs->size1() == 0 || this->GetReshapeMatrixFlag() == true) {
        rpLhs->resize(this->mEquationSystemSize, this->mEquationSystemSize, false);
        mpImpl->MakeLhsTopology(pScheme, *rpLhs, rModelPart);
    } else {
        if (rpLhs->size1() != this->mEquationSystemSize || rpLhs->size2() != this->mEquationSystemSize) {
            KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permitted."<<std::endl;
            rpLhs->resize(this->mEquationSystemSize, this->mEquationSystemSize, false);
            mpImpl->MakeLhsTopology(pScheme, *rpLhs, rModelPart);
        }
    }

    // Zero out system vectors.
    if (rpSolution->size() != this->mEquationSystemSize)
        rpSolution->resize(this->mEquationSystemSize, false);
    TSparse::SetToZero(*rpSolution);

    if (rpRhs->size() != this->mEquationSystemSize)
        rpRhs->resize(this->mEquationSystemSize, false);
    TSparse::SetToZero(*rpRhs);

    // Construct constraint topology.
    mpImpl->MakeConstraintTopology(rModelPart);

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Hooks
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse, TDense, TSolver>::InitializeSolutionStep(ModelPart& rModelPart,
                                                                                  typename Interface::TSystemMatrixType& rLhs,
                                                                                  typename Interface::TSystemVectorType& rSolution,
                                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    //
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Assembly
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Build(typename Interface::TSchemeType::Pointer pScheme,
                                                               ModelPart& rModelPart,
                                                               typename Interface::TSystemMatrixType& rLhs,
                                                               typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF(!pScheme) << "missing scheme" << std::endl;
    KRATOS_TRY
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/true>(rModelPart,
                                                                         *pScheme,
                                                                         &rLhs,
                                                                         &rRhs);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildLHS(typename Interface::TSchemeType::Pointer pScheme,
                                                                  ModelPart& rModelPart,
                                                                  typename Interface::TSystemMatrixType& rLhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    KRATOS_ERROR_IF(!pScheme) << "missing scheme" << std::endl;
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/false>(rModelPart,
                                                                          *pScheme,
                                                                          &rLhs,
                                                                          nullptr);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildRHS(typename Interface::TSchemeType::Pointer pScheme,
                                                                  ModelPart& rModelPart,
                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    mpImpl->template Assemble</*AssembleLHS=*/false,/*AssembleRHS=*/true>(rModelPart,
                                                                          *pScheme,
                                                                          nullptr,
                                                                          &rRhs);
    block_for_each(this->mDofSet, [&rRhs](Dof<double>& rDof){
        if (rDof.IsFixed())
            rRhs[rDof.EquationId()] = 0.0;
    });

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Constraint Imposition
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyDirichletConditions(typename Interface::TSchemeType::Pointer pScheme,
                                                                                  ModelPart& rModelPart,
                                                                                  typename Interface::TSystemMatrixType& rLhs,
                                                                                  typename Interface::TSystemVectorType& rSolution,
                                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    const std::size_t system_size = rLhs.size1();
    Vector scaling_factors (system_size);

    const auto it_dof_iterator_begin = this->mDofSet.begin();

    // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
    IndexPartition<std::size_t>(this->mDofSet.size()).for_each([&](std::size_t Index){
        auto it_dof_iterator = it_dof_iterator_begin + Index;
        if (it_dof_iterator->IsFixed()) {
            scaling_factors[Index] = 0.0;
        } else {
            scaling_factors[Index] = 1.0;
        }
    });

    //GetDiagonalScaleFactor<TSparse>(mpImpl->mDiagonalScaleFactor,
    //                                rLhs,
    //                                mpImpl->mDiagonalScaling,
    //                                rModelPart.GetProcessInfo());

    double* Avalues = rLhs.value_data().begin();
    std::size_t* Arow_indices = rLhs.index1_data().begin();
    std::size_t* Acol_indices = rLhs.index2_data().begin();

    IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end = Arow_indices[Index+1];
        const double k_factor = scaling_factors[Index];
        if (k_factor == 0.0) {
            // Zero out the whole row, except the diagonal
            for (std::size_t j = col_begin; j < col_end; ++j)
                if (Acol_indices[j] != Index )
                    Avalues[j] = 0.0;
            // Zero out the RHS
            rRhs[Index] = 0.0;
        } else {
            // Zero out the column which is associated with the zero'ed row
            for (std::size_t j = col_begin; j < col_end; ++j)
                if(scaling_factors[ Acol_indices[j] ] == 0 )
                    Avalues[j] = 0.0;
        }
    });
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyRHSConstraints(typename Interface::TSchemeType::Pointer pScheme,
                                                                             ModelPart& rModelPart,
                                                                             typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY

    if (rModelPart.MasterSlaveConstraints().size() != 0) {
        mpImpl->MakeConstraints(rModelPart);

        switch (mpImpl->mConstraintImposition) {
            case ConstraintImposition::MasterSlave: {
                // We compute the transposed matrix of the global relation matrix
                typename Interface::TSystemMatrixType T_transpose_matrix(mpImpl->mRelationMatrix.size2(),
                                                                        mpImpl->mRelationMatrix.size1());
                SparseMatrixMultiplicationUtility::TransposeMatrix<typename Interface::TSystemMatrixType, typename Interface::TSystemMatrixType>(T_transpose_matrix, mpImpl->mRelationMatrix, 1.0);

                typename Interface::TSystemVectorType modified_rhs(rRhs.size());
                TSparse::Mult(T_transpose_matrix, rRhs, modified_rhs);
                TSparse::Copy(modified_rhs, rRhs);

                // Apply diagonal values on slaves
                IndexPartition<std::size_t>(mpImpl->mSlaveIds.size()).for_each([&](std::size_t Index){
                    const IndexType slave_equation_id = mpImpl->mSlaveIds[Index];
                    if (mpImpl->mInactiveSlaveIds.find(slave_equation_id) == mpImpl->mInactiveSlaveIds.end()) {
                        rRhs[slave_equation_id] = 0.0;
                    }
                });

                break;
            } // case ConstraintImposition::MasterSlave

            case ConstraintImposition::AugmentedLagrange: {
                KRATOS_ERROR << "Not supported yet";
                break;
            } // ConstraintImposition::AugmentedLagrange

            default: {
                KRATOS_ERROR << "Unsupported constraint imposition method: " << (int)mpImpl->mConstraintImposition;
            }
        } // switch (mpImpl->mConstraintImposition)
    } // if rModelPart.MasterSlaveConstraints

    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyConstraints(typename Interface::TSchemeType::Pointer pScheme,
                                                                          ModelPart& rModelPart,
                                                                          typename Interface::TSystemMatrixType& rLhs,
                                                                          typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    mpImpl->ImposeConstraints(*pScheme,
                              rModelPart,
                              rLhs,
                              rRhs);
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Solution
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::SystemSolve(typename Interface::TSystemMatrixType& rLhs,
                                                                     typename Interface::TSystemVectorType& rSolution,
                                                                     typename Interface::TSystemVectorType& rRhs)
{
    throw std::runtime_error("");
    KRATOS_TRY
    this->mpLinearSystemSolver->Solve(rLhs, rSolution, rRhs);
    KRATOS_CATCH("")

    KRATOS_TRY
    switch (mpImpl->mConstraintImposition) {
        case ConstraintImposition::MasterSlave: {
            if(mpImpl->mRelationMatrix.size1() != 0) { // If there are master-slave constraints
                // Recover solution of the original problem
                typename Interface::TSystemVectorType modified_solution = rSolution;

                // Recover solution of the original problem
                TSparse::Mult(mpImpl->mRelationMatrix, modified_solution, rSolution);
            }
            break;
        } // case ConstraintImposition::MasterSlave

        case ConstraintImposition::AugmentedLagrange: {
            KRATOS_ERROR << "Not implemented yet";
            break;
        } // case ConstraintImposition::AugmentedLagrange

        default: {
            KRATOS_ERROR << "Unsupported constraint imposition method: " << (int)mpImpl->mConstraintImposition;
        }
    } // switch constraint imposition
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Compound Assembly and Solution
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildAndSolve(typename Interface::TSchemeType::Pointer pScheme,
                                                                       ModelPart& rModelPart,
                                                                       typename Interface::TSystemMatrixType& rLhs,
                                                                       typename Interface::TSystemVectorType& rSolution,
                                                                       typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY

    // Assemble unconstrained system.
    Build(pScheme, rModelPart, rLhs, rRhs);

    // Apply multifreedom constraints.
    ApplyConstraints(pScheme, rModelPart, rLhs, rRhs);

    // Apply Dirichlet conditions.
    ApplyDirichletConditions(pScheme, rModelPart, rLhs, rSolution, rRhs);

    // Solve constrained assembled system.
    mpImpl->Solve(rLhs, rSolution, rRhs, rModelPart, *this);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildRHSAndSolve(typename Interface::TSchemeType::Pointer pScheme,
                                                                          ModelPart& rModelPart,
                                                                          typename Interface::TSystemMatrixType& rLhs,
                                                                          typename Interface::TSystemVectorType& rSolution,
                                                                          typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    BuildRHS(pScheme, rModelPart, rRhs);

    if(rModelPart.MasterSlaveConstraints().size() != 0) {
        ApplyRHSConstraints(pScheme, rModelPart, rRhs);
    }

    ApplyDirichletConditions(pScheme, rModelPart, rLhs, rSolution, rRhs);
    mpImpl->Solve(rLhs, rSolution, rRhs, rModelPart, *this);
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Postprocessing
// --------------------------------------------------------- //



template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::CalculateReactions(typename Interface::TSchemeType::Pointer pScheme,
                                                                            ModelPart& rModelPart,
                                                                            typename Interface::TSystemMatrixType& rLhs,
                                                                            typename Interface::TSystemVectorType& rSolution,
                                                                            typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    TSparse::SetToZero(rRhs);
    mpImpl->template Assemble</*AssembleLHS=*/false,/*AssembleRHS=*/true>(rModelPart,
                                                                          *pScheme,
                                                                          nullptr,
                                                                          &rRhs);
    block_for_each(this->mDofSet, [&rRhs](Dof<double>& rDof){
        rDof.GetSolutionStepReactionValue() = -rRhs[rDof.EquationId()];
    });
}


// --------------------------------------------------------- //
// Misc
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::AssignSettings(const Parameters Settings)
{
    KRATOS_TRY
    Settings.ValidateDefaults(this->GetDefaultParameters());
    KRATOS_CATCH("")

    KRATOS_TRY
    Interface::AssignSettings(Settings);
    KRATOS_CATCH("")

    // Set the scaling strategy for the diagonal entries of constrained DoFs.
    KRATOS_TRY
    if (Settings["diagonal_scaling"].Is<std::string>()) {
        const auto diagonal_scaling_strategy = Settings["diagonal_scaling"].Get<std::string>();
        if (diagonal_scaling_strategy == "none") {
            mpImpl->mDiagonalScaling = DiagonalScaling::None;
        } else if (diagonal_scaling_strategy == "abs_max") {
            mpImpl->mDiagonalScaling = DiagonalScaling::AbsMax;
        } else if (diagonal_scaling_strategy == "norm") {
            mpImpl->mDiagonalScaling = DiagonalScaling::Norm;
        } else {
            KRATOS_ERROR << "unsupported setting for \"diagonal_scaling\": "
                         << diagonal_scaling_strategy << ". Options are:\n"
                         << "- \"none\"\n"
                         << "- \"abs_max\"\n"
                         << "- \"norm\""
                         << "- a floating point constant\n";
        }
    } /*if Settings["diagonal_scaling"].Is<std::string>()*/ else if (Settings["diagonal_scaling"].IsNumber()) {
        mpImpl->mDiagonalScaling = DiagonalScaling::Constant;
        mpImpl->mDiagonalScaleFactor = Settings["diagonal_scaling"].Get<double>();
    }
    KRATOS_CATCH("")

    // Set multifreedom constraint imposition strategy.
    KRATOS_TRY
    Parameters constraint_imposition_settings = Settings["constraint_imposition"];
    const std::string constraint_imposition_name = constraint_imposition_settings["method"].Get<std::string>();
    if (constraint_imposition_name == "master_slave_elimination") {
        mpImpl->mConstraintImposition = ConstraintImposition::MasterSlave;
    } else if (constraint_imposition_name == "augmented_lagrange") {
        Parameters default_parameters(R"({"method" : "augmented_lagrange", "penalty_factor" : 1e3})");
        constraint_imposition_settings.ValidateAndAssignDefaults(default_parameters);
        mpImpl->mConstraintImposition = ConstraintImposition::AugmentedLagrange;
        mpImpl->mPenaltyFactor = constraint_imposition_settings["penalty_factor"].Get<double>();
    } else {
        std::stringstream message;
        message << "Unsupported constraint imposition \"" << constraint_imposition_name << "\". Options are:\n";
        message << "\t\"master_slave_elimination\n\"";
        message << "\t\"augmented_lagrange\"";
        KRATOS_ERROR << message.str();
    }
    KRATOS_CATCH("")

    // Other settings.
    KRATOS_TRY
    mpImpl->mMaxDepth = Settings["max_depth"].Get<int>();
    mpImpl->mVerbosity = Settings["verbosity"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Clear()
{
    Interface::Clear();
    mpImpl->mSlaveIds = decltype(mpImpl->mSlaveIds)();
    mpImpl->mMasterIds = decltype(mpImpl->mMasterIds)();
    mpImpl->mInactiveSlaveIds = decltype(mpImpl->mInactiveSlaveIds)();
    mpImpl->mRelationMatrix = decltype(mpImpl->mRelationMatrix)();
    mpImpl->mConstraintGapVector = decltype(mpImpl->mConstraintGapVector)();
}


template <class TSparse, class TDense, class TSolver>
Parameters PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetDefaultParameters() const
{
    Parameters parameters = Parameters(R"({
"name"              : "p_multigrid",
"diagonal_scaling"  : "abs_max",
"max_depth"         : -1,
"smoother_settings" : {
    "solver_type" : "gauss_seidel"
},
"solver_settings"   : {
    "solver_type" : "amgcl"
},
"constraint_imposition" : {
    "method" : "master_slave_elimination"
},
"verbosity"         : 0
})");
    parameters.RecursivelyAddMissingParameters(Interface::GetDefaultParameters());
    return parameters;
}


template <class TSparse, class TDense, class TSolver>
std::string PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Info() const
{
    return "PMultigridBuilderAndSolver";
}


template <class TSparse, class TDense, class TSolver>
std::size_t PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetEquationSystemSize() const noexcept
{
    return Interface::mEquationSystemSize;
}

template <class TSparse, class TDense, class TSolver>
TSolver& PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetLinearSolver() noexcept
{
    return *Interface::mpLinearSystemSolver;
}


// --------------------------------------------------------- //
// Template Instantiations
// --------------------------------------------------------- //


template class PMultigridBuilderAndSolver<TUblasSparseSpace<double>,
                                          TUblasDenseSpace<double>,
                                          LinearSolver<TUblasSparseSpace<double>,
                                                       TUblasDenseSpace<double>>>;


} // namespace Kratos
