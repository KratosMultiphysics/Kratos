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
#include <unordered_map> // std::unordered_map


namespace Kratos::detail {


inline void ProcessMasterSlaveConstraint(std::vector<std::size_t>& rConstraintIndices,
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


inline void ProcessMultifreedomConstraint(std::vector<std::size_t>& rConstraintIndices,
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


template <class TSparse, class TDense>
void MakeRelationTopology(std::size_t SystemSize,
                          const typename ConstraintAssembler<TSparse,TDense>::ConstraintArray& rConstraints,
                          const ProcessInfo& rProcessInfo,
                          typename TSparse::MatrixType& rRelationMatrix,
                          typename TSparse::VectorType& rConstraintGaps,
                          std::unordered_map<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    KRATOS_TRY

    // Build constraint ID => index map.
    rConstraintIdMap.clear();

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
                const auto emplace_result = rConstraintIdMap.emplace(static_cast<std::size_t>(constraint_label),
                                                                     std::make_pair(i_constraint, 0ul));
                ++emplace_result.first->second.second; //< increment the number of constraint objects defining the constraint equation.
                if (emplace_result.second) ++i_constraint;
            } // for constraint_label in constraint_labels
        } // for r_constraint in rModelPart.MasterSlaveConstraints
    }

    {
        std::vector<std::unordered_set<IndexType>> indices(rConstraintIdMap.size());
        std::vector<LockObject> mutexes(rConstraintIdMap.size());

        struct TLS {
            MasterSlaveConstraint::EquationIdVectorType slaves, masters;
            std::vector<std::size_t> constraint_labels;
        };

        block_for_each(rConstraints,
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
                            typename TSparse::VectorType& rConstraintGaps,
                            std::unordered_map<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)
{
    KRATOS_TRY

    // Function-wide variables.
    std::vector<LockObject> mutexes(rConstraintIdMap.size());

    // Init.
    TSparse::SetToZero(rRelationMatrix);
    TSparse::SetToZero(rConstraintGaps);

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
                   [&mutexes, &rProcessInfo, &rConstraintIdMap, &rRelationMatrix, &rConstraintGaps](const MasterSlaveConstraint& r_constraint, TLS& r_tls){
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
            } // not r_tls.slave_ids.empty()

            r_tls.relation_matrix_row.resize(r_tls.relation_matrix.size2());

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
                    MapRowContribution<TSparse,TDense>(rRelationMatrix,
                                                       r_tls.relation_matrix,
                                                       i_constraint,
                                                       i_row,
                                                       r_tls.reordered_dof_equation_ids);
                }

                AtomicAdd(rConstraintGaps[i_constraint],
                          static_cast<typename TSparse::DataType>(r_tls.constraint_gaps[i_row]));
            } // for i_row in range(r_tls.slave_ids.size)
        } // if r_constraint.IsActive
    }); // for r_constraint in rModelPart.MasterSlaveCosntraints

    KRATOS_CATCH("")
}


} // namespace Kratos::detail
