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
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility
#include "utilities/reduction_utilities.h"
#include "utilities/proxies.h" // MakeProxy
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE, KRATOS_PROFILE_SCOPE_MILLI
#include "solving_strategies/builder_and_solvers/p_multigrid_utilities.hpp" // MakePRestrictionOperator, MakeSparseTopology

// System includes
#include <algorithm> // std::lower_bound
#include <optional> // std::optional
#include <limits> // std::numeric_limits


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
    None        = 0,
    MasterSlave = 1,
    Lagrange    = 2,
    Penalty     = 3
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

            mpLhs->resize(mpRestrictionOperator->size1(), mpRestrictionOperator->size1());
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

            mpRhs->resize(mpRestrictionOperator->size1());
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

    Impl(Interface* pInterface)
        : mpInterface(pInterface),
          mRelationMatrix(),
          mConstraintGapVector(),
          mSlaveIds(),
          mMasterIds(),
          mInactiveSlaveIds(),
          mDiagonalScaleFactor(1),
          mHierarchy(),
          mDiagonalScaling(DiagonalScaling::None),
          mMaxDepth(-1),
          mVerbosity(0)
    {}

    /// @brief Initialize the linear solver and solve the provided system.
    void Solve(typename Interface::TSystemMatrixType& rLhs,
               typename Interface::TSystemVectorType& rSolution,
               typename Interface::TSystemVectorType& rRhs,
               ModelPart& rModelPart,
               Interface& rInterface)
    {
        if(!rModelPart.MasterSlaveConstraints().empty()) {
            typename Interface::TSystemVectorType modified_solution(rRhs.size());
            TSparse::SetToZero(modified_solution);

            if (rInterface.GetLinearSolver().AdditionalPhysicalDataIsNeeded()) {
                rInterface.GetLinearSolver().ProvideAdditionalData(rLhs,
                                                                   modified_solution,
                                                                   rRhs,
                                                                   rInterface.mDofSet,
                                                                   rModelPart);
            }
            rInterface.GetLinearSolver().Solve(rLhs, modified_solution, rRhs);

            //recover solution of the original problem
            TSparse::Mult(mRelationMatrix, modified_solution, rSolution);
        } else {
            if (rInterface.GetLinearSolver().AdditionalPhysicalDataIsNeeded()) {
                rInterface.GetLinearSolver().ProvideAdditionalData(rLhs,
                                                                   rSolution,
                                                                   rRhs,
                                                                   rInterface.mDofSet,
                                                                   rModelPart);
            }
            rInterface.GetLinearSolver().Solve(rLhs, rSolution, rRhs);
        }
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
    }


    void MakeLhsTopology(const typename Interface::TSchemeType::Pointer& rpScheme,
                         typename Interface::TSystemMatrixType& rLhs,
                         ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        std::vector<LockObject> locks(mpInterface->GetEquationSystemSize());
        std::vector<std::unordered_set<std::size_t> > indices(mpInterface->GetEquationSystemSize());

        // Collect DoFs from elements.
        Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Element>(rModelPart),
                          rModelPart.GetProcessInfo(),
                          *rpScheme,
                          locks.data(),
                          locks.data() + locks.size(),
                          indices.data(),
                          indices.data() + indices.size());

        // Collect DoFs from conditions.
        Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Condition>(rModelPart),
                          rModelPart.GetProcessInfo(),
                          *rpScheme,
                          locks.data(),
                          locks.data() + locks.size(),
                          indices.data(),
                          indices.data() + indices.size());

        // Collect DoFs from multifreedom constraints.
        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            struct TLS
            {
                Element::EquationIdVectorType master_ids, slave_ids;
            };

            block_for_each(rModelPart.MasterSlaveConstraints(), TLS(), [&](MasterSlaveConstraint& rConst, TLS& rTls){
                rConst.EquationIdVector(rTls.slave_ids, rTls.master_ids, r_current_process_info);

                for (std::size_t i = 0; i < rTls.slave_ids.size(); i++) {
                    auto& row_indices = indices[rTls.slave_ids[i]];
                    std::scoped_lock<LockObject> lock(locks[rTls.slave_ids[i]]);
                    row_indices.insert(rTls.slave_ids[i]);
                }

                for (std::size_t i = 0; i < rTls.master_ids.size(); i++) {
                    auto& row_indices = indices[rTls.master_ids[i]];
                    std::scoped_lock<LockObject> lock(locks[rTls.master_ids[i]]);
                    row_indices.insert(rTls.master_ids[i]);
                }
            });
        } // if rModelPart.MasterSlaveConstraints

        // Destroy locks.
        locks = std::vector<LockObject>();

        // Compute and allocate LHS topology.
        MakeSparseTopology<typename TSparse::DataType>(indices,
                                                       indices.size(),
                                                       rLhs);

        // Construct the coarse hierarhy's topology.
        mHierarchy.template MakeLhsTopology<TSparse>(rModelPart, rLhs);
        KRATOS_CATCH("")
    }


    void MakeConstraintTopology(ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY
        if (!rModelPart.MasterSlaveConstraints().empty()) {
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // Global variables.
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            std::vector<std::unordered_set<IndexType>> indices(mpInterface->mDofSet.size());
            std::vector<LockObject> locks(indices.size());

            #pragma omp parallel
            {
                Element::EquationIdVectorType slave_ids;
                Element::EquationIdVectorType master_ids;
                std::unordered_map<IndexType, std::unordered_set<IndexType>> temp_indices;

                #pragma omp for schedule(guided, 512) nowait
                for (int i_const = 0; i_const < static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i_const) {
                    auto it_const = it_const_begin + i_const;
                    it_const->EquationIdVector(slave_ids, master_ids, r_current_process_info);

                    // Slave DoFs
                    for (auto &id_i : slave_ids) {
                        temp_indices[id_i].insert(master_ids.begin(), master_ids.end());
                    }
                }

                // Merging temporary indices into global rows.
                for (auto& pair_temp_indices : temp_indices) {
                    std::scoped_lock<LockObject> lock(locks[pair_temp_indices.first]);
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
            mConstraintGapVector.resize(indices.size());
        } // if rModelPart.MasterSlaveConstraints()
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


    void MakeConstraints(ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mRelationMatrix.size1() == mConstraintGapVector.size());
        TSparse::SetToZero(mRelationMatrix);
        TSparse::SetToZero(mConstraintGapVector);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Declare local containers.
        typename Interface::LocalSystemMatrixType local_relation_matrix;
        typename Interface::LocalSystemVectorType local_constraint_gap_vector;
        Element::EquationIdVectorType slave_equation_ids, master_equation_ids;
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        // We clear the set
        mInactiveSlaveIds.clear();

        #pragma omp parallel firstprivate(local_relation_matrix, local_constraint_gap_vector, slave_equation_ids, master_equation_ids)
        {
            std::unordered_set<IndexType> tls_inactive_slave_dofs;

            #pragma omp for schedule(guided, 512)
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;
                it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);

                // If the constraint is active
                if (it_const->IsActive()) {
                    it_const->CalculateLocalSystem(local_relation_matrix,
                                                   local_constraint_gap_vector,
                                                   r_current_process_info);

                    for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                        const IndexType i_global = slave_equation_ids[i];

                        // Assemble matrix row
                        Impl::MapRowContribution(mRelationMatrix,
                                                 local_relation_matrix,
                                                 i_global,
                                                 i,
                                                 master_equation_ids);

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
            mConstraintGapVector[eq_id] = 0.0;
            mRelationMatrix(eq_id, eq_id) = 1.0;
        }

        // Setting inactive slave dofs in the T and C system
        for (auto eq_id : mInactiveSlaveIds) {
            mConstraintGapVector[eq_id] = 0.0;
            mRelationMatrix(eq_id, eq_id) = 1.0;
        }

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
    // Member Variables
    // --------------------------------------------------------- //

    Interface* mpInterface;

    typename Interface::TSystemMatrixType mRelationMatrix;

    typename Interface::TSystemVectorType mConstraintGapVector;

    std::vector<std::size_t> mSlaveIds;

    std::vector<std::size_t> mMasterIds;

    std::unordered_set<std::size_t> mInactiveSlaveIds;

    typename TSparse::DataType mDiagonalScaleFactor;

    PGrid<TSparse,TDense,TSolver> mHierarchy;

    DiagonalScaling mDiagonalScaling;

    int mMaxDepth;

    int mVerbosity;
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
            rpLhs->resize(this->mEquationSystemSize, this->mEquationSystemSize, true);
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
    }

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

    if (!rModelPart.MasterSlaveConstraints().empty()) {
        mpImpl->MakeConstraints(rModelPart);

        // We compute the transposed matrix of the global relation matrix
        {
            // Storage for an intermediate matrix transpose(relation_matrix) * lhs_matrix
            typename Interface::TSystemMatrixType left_multiplied_lhs(mpImpl->mRelationMatrix.size2(), rLhs.size2());
            {
                // Transpose the relation matrix
                typename Interface::TSystemMatrixType transposed_relation_matrix(mpImpl->mRelationMatrix.size2(), mpImpl->mRelationMatrix.size1());
                SparseMatrixMultiplicationUtility::TransposeMatrix(transposed_relation_matrix, mpImpl->mRelationMatrix, 1.0);

                typename Interface::TSystemVectorType b_modified(rRhs.size());
                TSparse::Mult(transposed_relation_matrix, rRhs, b_modified);
                TSparse::Copy(b_modified, rRhs);

                SparseMatrixMultiplicationUtility::MatrixMultiplication(transposed_relation_matrix, rLhs, left_multiplied_lhs);
            } // deallocate transposed_relation_matrix

            SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, mpImpl->mRelationMatrix, rLhs);
        } // deallocate left_multiplied_lhs

        // Compute the scale factor for slave DoFs.
        Impl::GetDiagonalScaleFactor(mpImpl->mDiagonalScaleFactor,
                                     rLhs,
                                     mpImpl->mDiagonalScaling,
                                     rModelPart.GetProcessInfo());

        // Apply diagonal values on slaves
        block_for_each(mpImpl->mSlaveIds, [this, &rLhs, &rRhs](const auto iSlave){
            if (mpImpl->mInactiveSlaveIds.find(iSlave) == mpImpl->mInactiveSlaveIds.end()) {
                rLhs(iSlave, iSlave) = mpImpl->mDiagonalScaleFactor;
                rRhs[iSlave] = 0.0;
            }
        });
    } // if !rModelPart.MasterSlaveConstraints().empty()

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
    if(mpImpl->mRelationMatrix.size1() != 0) { // If there are master-slave constraints
        // Recover solution of the original problem
        typename Interface::TSystemVectorType modified_solution = rSolution;

        // Recover solution of the original problem
        TSparse::Mult(mpImpl->mRelationMatrix, modified_solution, rSolution);
    }
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
