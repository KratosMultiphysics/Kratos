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

#pragma once

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler.hpp" // ConstraintAssembler
#include "solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp" // MakeSparseTopology
#include "includes/master_slave_constraint.h" // MasterSlaveConstraint
#include "includes/variables.h" // CONSTRAINT_LABELS

// STL includes
#include <vector> // std::vector


namespace Kratos::detail {


/// @brief Compute the constraint gradient and gaps of a @ref LinearMasterSlaveConstraint.
/// @details This utility function computes the constraint gradient of constraint gaps of a
///          constraint that does not inherit from @ref MultifreedomConstraint. It also collects
///          the constraint indices of the returned entries, as well as the IDs of the @ref Dof "DoFs"
///          involved.
/// @param rConstraintIndices Output vector containing the indices of each returned constraint equation.
/// @param rDofIds Output vector of @ref Dof IDs that participate in the returned constraint equations.
/// @param rRelationMatrix Output matrix containing the gradients of each constraint equation defined by
///                        the input constraint. Each row represents a gradient of the constraint equation
///                        whose constraint index is defined by the corresponding entry in @p rConstraintIndices.
///                        Columns represent coefficients of DoFs that are identified by their IDs at the
///                        corresponding position in @p rDofIds.
/// @param rConstraintGaps Output vector of constraint gaps related to constraint equations defined by
///                        @p rConstraintIndices.
/// @param rConstraint Constraint to extract values from. It must not be an instance inheriting from
///                    @ref MultifreedomConstraint.
/// @param rSlaveDofIds List of slave @ref Dof "DoFs"' IDs.
/// @param rMasterDofIds List of master @ref Dof "DoFs"' IDs.
/// @param rConstraintIdMap Map relating slave DoF IDs to constraint equation indices, as well as the number
///                         of constraint objects defining the constraint equation.
/// @see ProcessMultifreedomConstraint
inline void ProcessMasterSlaveConstraint(std::vector<std::size_t>& rConstraintIndices,
                                         std::vector<std::size_t>& rDofIds,
                                         MasterSlaveConstraint::MatrixType& rRelationMatrix,
                                         [[maybe_unused]] MasterSlaveConstraint::VectorType& rConstraintGaps,
                                         const MasterSlaveConstraint& rConstraint,
                                         const std::vector<std::size_t>& rSlaveDofIds,
                                         const std::vector<std::size_t>& rMasterDofIds,
                                         const CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    // Constraint identifiers are the slave DoFs' IDs.
    rConstraintIndices.resize(rSlaveDofIds.size());
    std::transform(rSlaveDofIds.begin(),
                   rSlaveDofIds.end(),
                   rConstraintIndices.begin(),
                   [&rConstraintIdMap](std::size_t slave_id){
                        return rConstraintIdMap.at(slave_id).first;
                   });

    // Comb DoFs.
    rDofIds.resize(rMasterDofIds.size() + rSlaveDofIds.size());
    std::copy(rMasterDofIds.begin(), rMasterDofIds.end(), rDofIds.begin());
    std::copy(rSlaveDofIds.begin(), rSlaveDofIds.end(), rDofIds.begin() + rMasterDofIds.size());

    // Reconstruct the constraint equations.
    if (rRelationMatrix.size2()) {
        rRelationMatrix.resize(rRelationMatrix.size1(),
                               rRelationMatrix.size2() + rSlaveDofIds.size(),
                               /*preserve entries=*/true);
        for (std::size_t i_slave=0ul; i_slave<rSlaveDofIds.size(); ++i_slave) {
            // Fetch the number of constraint objects defining this constraint equation,
            // because the slave DoF's coefficient must add up to -1 after assembly.
            const auto constraint_object_count = rConstraintIdMap.at(rSlaveDofIds[i_slave]).second;
            rRelationMatrix(i_slave, rMasterDofIds.size() + i_slave) = -1.0 / constraint_object_count;
        } // for i_slave in range(rSlaveDofIds.size())
    } // if rRelationMatrix.size2()

    // No need to move the constraint gaps from RHS to LHS
    // because the slaves were moved to the RHS instead:
    // u^s = T u^m + b
    // =>
    // 0 = A [u^m \\ u^s] + b
}


/// @brief Collect @ref Dof IDs and constraint equation indices of a @ref MultifreedomConstraint.
/// @param rConstraintIndices Output vector containing the indices of each returned constraint equation.
/// @param rDofIds Output vector of @ref Dof IDs that participate in the returned constraint equations.
/// @param rConstraint Constraint to extract values from. It must not be an instance inheriting from
///                    @ref MultifreedomConstraint.
/// @param rMasterDofIds List of master @ref Dof "DoFs"' IDs.
/// @param rConstraintIdMap Map relating slave DoF IDs to constraint equation indices, as well as the number
///                         of constraint objects defining the constraint equation.
/// @see ProcessMasterSlaveConstraint
inline void ProcessMultifreedomConstraint(std::vector<std::size_t>& rConstraintIndices,
                                          std::vector<std::size_t>& rDofIds,
                                          const MasterSlaveConstraint& rConstraint,
                                          const std::vector<std::size_t>& rMasterDofIds,
                                          const CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    const auto& r_constraint_labels = rConstraint.GetData().GetValue(CONSTRAINT_LABELS);
    rConstraintIndices.resize(r_constraint_labels.size());
    std::transform(r_constraint_labels.begin(),
                   r_constraint_labels.end(),
                   rConstraintIndices.begin(),
                   [&rConstraintIdMap](std::size_t constraint_label){
                        return rConstraintIdMap.at(constraint_label).first;
                   });

    rDofIds = rMasterDofIds;
}


template <class TSparse, class TDense>
void MakeRelationTopology(std::size_t SystemSize,
                          typename ConstraintAssembler<TSparse,TDense>::ConstraintArray::const_iterator itConstraintBegin,
                          typename ConstraintAssembler<TSparse,TDense>::ConstraintArray::const_iterator itConstraintEnd,
                          const ProcessInfo& rProcessInfo,
                          typename TSparse::MatrixType& rRelationMatrix,
                          typename TSparse::VectorType& rConstraintGaps,
                          CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    KRATOS_TRY

    // Build constraint ID => index map.
    rConstraintIdMap.clear();

    {
        MasterSlaveConstraint::IndexType i_constraint = 0;
        MasterSlaveConstraint::EquationIdVectorType constraint_labels, master_ids;

        for (auto it_constraint=itConstraintBegin; it_constraint!=itConstraintEnd; ++it_constraint) {
            const auto& r_constraint = *it_constraint;
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
                const auto emplace_result = rConstraintIdMap.emplace(static_cast<std::size_t>(constraint_label),
                                                                     std::make_pair(i_constraint, 0ul));
                ++emplace_result.first->second.second; //< increment the number of constraint objects defining the constraint equation.
                if (emplace_result.second) ++i_constraint;
            } // for constraint_label in constraint_labels
        } // for r_constraint in rModelPart.MasterSlaveConstraints
    }

    {
        std::vector<CSRHashSet<IndexType>> indices(rConstraintIdMap.size());
        std::vector<LockObject> mutexes(rConstraintIdMap.size());

        struct TLS {
            MasterSlaveConstraint::EquationIdVectorType slaves, masters;
            std::vector<std::size_t> constraint_labels;
        };

        block_for_each(itConstraintBegin,
                       itConstraintEnd,
                       TLS(),
                       [&mutexes, &indices, &rProcessInfo, &rConstraintIdMap](const auto& r_constraint, TLS& r_tls) {
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
                const auto i_constraint = rConstraintIdMap[static_cast<std::size_t>(i_slave)].first;
                std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                indices[i_constraint].insert(r_tls.masters.begin(), r_tls.masters.end());
            } // for i_slave in slave_ids
        }); // for r_constraint in rModelPart.MasterSlaveConstraints()

        MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                             SystemSize,
                                                             rRelationMatrix,
                                                             /*EnsureDiagonal=*/false);
        TSparse::SetToZero(rRelationMatrix);
        rConstraintGaps.resize(rRelationMatrix.size1(), false);
        //TSparse::SetToZero(rConstraintGaps); //< unnecessary
    }

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void AssembleRelationMatrix(const typename ConstraintAssembler<TSparse,TDense>::ConstraintArray& rConstraints,
                            const ProcessInfo& rProcessInfo,
                            typename TSparse::MatrixType& rRelationMatrix,
                            typename TSparse::MatrixType& rHessian,
                            typename TSparse::VectorType& rConstraintGaps,
                            CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    KRATOS_TRY

    // Function-wide variables.
    std::vector<LockObject> constraint_mutexes(rConstraintIdMap.size());
    std::vector<LockObject> dof_mutexes(rHessian.size1());

    // Init.
    TSparse::SetToZero(rRelationMatrix);
    TSparse::SetToZero(rHessian);
    TSparse::SetToZero(rConstraintGaps);

    struct TLS {
        MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
        std::vector<std::size_t> constraint_indices, dof_equation_ids, dof_index_array, reordered_dof_equation_ids;
        std::vector<MasterSlaveConstraint::MatrixType::value_type> matrix_row;
        MasterSlaveConstraint::MatrixType relation_matrix, hessian;
        MasterSlaveConstraint::VectorType constraint_gaps;
    };

    // Constraint assembly.
    block_for_each(rConstraints,
                   TLS(),
                   [&constraint_mutexes,
                    &dof_mutexes,
                    &rProcessInfo,
                    &rConstraintIdMap,
                    &rRelationMatrix,
                    &rHessian,
                    &rConstraintGaps](const MasterSlaveConstraint& r_constraint,
                                                                                                               TLS& r_tls){
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
                                                      rConstraintIdMap);
                const auto& r_hessian = r_constraint.GetData().GetValue(GEOMETRIC_STIFFNESS_MATRIX);
                r_tls.hessian.resize(r_hessian.size1(), r_hessian.size2());
                std::copy(r_hessian.data().begin(),
                          r_hessian.data().end(),
                          r_tls.hessian.data().begin());
            } /*if r_tls.slave_ids.empty()*/ else {
                // The constraint is a MasterSlaveConstraint.
                detail::ProcessMasterSlaveConstraint(r_tls.constraint_indices,
                                                     r_tls.dof_equation_ids,
                                                     r_tls.relation_matrix,
                                                     r_tls.constraint_gaps,
                                                     r_constraint,
                                                     r_tls.slave_ids,
                                                     r_tls.master_ids,
                                                     rConstraintIdMap);

                // The standard MasterSlaveConstraint can only represent linear constraints,
                // whose Hessian vanish, so there's no need assemble them into the global
                // constraint Hessian.
                r_tls.hessian.resize(0, 0);
            } // not r_tls.slave_ids.empty()

            KRATOS_ERROR_IF_NOT(r_tls.constraint_indices.size() == r_tls.relation_matrix.size1())
                << "constraint " << r_constraint.Id() << " is ill-formed: "
                << "has " << r_tls.constraint_indices.size() << " constraint equations, but its relation matrix is "
                << r_tls.relation_matrix.size1() << "x" << r_tls.relation_matrix.size2();

            KRATOS_ERROR_IF_NOT(r_tls.dof_equation_ids.size() == r_tls.relation_matrix.size2())
                << "constraint " << r_constraint.Id() << " is ill-formed: "
                << "defined on " << r_tls.dof_equation_ids.size() << " DoFs, but its relation matrix is "
                << r_tls.relation_matrix.size1() << "x" << r_tls.relation_matrix.size2();

            KRATOS_ERROR_IF_NOT(r_tls.constraint_gaps.size() == r_tls.relation_matrix.size1())
                << "constraint " << r_constraint.Id() << " is ill formed: "
                << "relation matrix is " << r_tls.relation_matrix.size1() << "x" << r_tls.relation_matrix.size2()
                << " but the constraint gap vector is of size " << r_tls.constraint_gaps.size();

            r_tls.matrix_row.resize(r_tls.relation_matrix.size2());

            // Sort DoFs based on their equation IDs because the CSR format expects rows to be sorted.
            r_tls.dof_index_array.resize(r_tls.dof_equation_ids.size());
            std::iota(r_tls.dof_index_array.begin(), r_tls.dof_index_array.end(), 0ul);
            std::sort(r_tls.dof_index_array.begin(),
                      r_tls.dof_index_array.end(),
                      [&r_tls](const std::size_t i_left, const std::size_t i_right) -> bool {
                          return r_tls.dof_equation_ids[i_left] < r_tls.dof_equation_ids[i_right];
                      });

            r_tls.reordered_dof_equation_ids.resize(r_tls.dof_equation_ids.size());
            std::transform(r_tls.dof_index_array.begin(),
                           r_tls.dof_index_array.end(),
                           r_tls.reordered_dof_equation_ids.begin(),
                           [&r_tls](const auto i_column){
                               return r_tls.dof_equation_ids[i_column];
                           });

            // Assemble local rows into the global relation matrix.
            for (std::size_t i_row=0ul; i_row<r_tls.constraint_indices.size(); ++i_row) {
                const auto i_constraint = r_tls.constraint_indices[i_row];
                KRATOS_ERROR_IF(rRelationMatrix.size1() <= i_constraint)
                    << "constraint index " << i_constraint
                    << " is out of bounds in relation matrix of size "
                    << "(" << rRelationMatrix.size1() << "x" << rRelationMatrix.size2() << ")";

                // Indirect sort the local relation matrix' row based on DoF IDs.
                // This step is required because the global relation matrix is in CSR format,
                // and Impl::MapRowContribution expects the column indices to be sorted.
                {
                    // Reorder the current row in the local relation matrix.
                    std::transform(r_tls.dof_index_array.begin(),
                                   r_tls.dof_index_array.end(),
                                   r_tls.matrix_row.begin(),
                                   [&r_tls, i_row](const std::size_t i_column){
                                        return r_tls.relation_matrix(i_row, i_column);
                                   });
                    std::copy(r_tls.matrix_row.begin(),
                              r_tls.matrix_row.end(),
                              (r_tls.relation_matrix.begin1() + i_row).begin());
                }

                {
                    std::scoped_lock<LockObject> lock(constraint_mutexes[i_constraint]);
                    MapRowContribution<TSparse,TDense>(rRelationMatrix,
                                                       r_tls.relation_matrix,
                                                       i_constraint,
                                                       i_row,
                                                       r_tls.reordered_dof_equation_ids);
                }

                AtomicAdd(rConstraintGaps[i_constraint],
                          static_cast<typename TSparse::DataType>(r_tls.constraint_gaps[i_row]));
            } // for i_row in range(r_tls.slave_ids.size)

            // Map contributions to the Hessian.
            if (r_tls.hessian.size1() && r_tls.hessian.size2()) {
                for (std::size_t i_row=0ul; i_row<r_tls.hessian.size1(); ++i_row) {
                    const auto i_reordered_row = r_tls.dof_index_array[i_row];
                    const auto i_row_dof = r_tls.reordered_dof_equation_ids[i_reordered_row];

                    // Reorder the current row in the local hessian.
                    std::transform(r_tls.dof_index_array.begin(),
                                   r_tls.dof_index_array.end(),
                                   r_tls.matrix_row.begin(),
                                   [&r_tls, i_reordered_row](const std::size_t i_column){
                                       return r_tls.hessian(i_reordered_row, i_column);
                                   });

                    std::copy(r_tls.matrix_row.begin(),
                              r_tls.matrix_row.end(),
                              (r_tls.hessian.begin1() + i_reordered_row).begin());

                    std::scoped_lock<LockObject> lock(dof_mutexes[i_row_dof]);
                    MapRowContribution<TSparse,TDense>(rHessian,
                                                       r_tls.hessian,
                                                       i_row_dof,
                                                       i_reordered_row,
                                                       r_tls.reordered_dof_equation_ids);
                } // for i_row in range(r_tls.hessian.size1())
            } // if the hessian does not vanish
        } // if r_constraint.IsActive
    }); // for r_constraint in rModelPart.MasterSlaveConstraints

    KRATOS_CATCH("")
}


} // namespace Kratos::detail
