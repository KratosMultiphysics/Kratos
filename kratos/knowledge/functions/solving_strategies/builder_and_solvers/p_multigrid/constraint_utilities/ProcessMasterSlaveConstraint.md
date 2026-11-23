---
type: C++ Function
title: ProcessMasterSlaveConstraint
resource: solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp#L49-L90
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/AssembleRelationMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void ProcessMasterSlaveConstraint(std::vector<std::size_t>& rConstraintIndices, std::vector<std::size_t>& rDofIds, MasterSlaveConstraint::MatrixType& rRelationMatrix, [[maybe_unused]] MasterSlaveConstraint::VectorType& rConstraintGaps, const MasterSlaveConstraint& rConstraint, const std::vector<std::size_t>& rSlaveDofIds, const std::vector<std::size_t>& rMasterDofIds, const CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)`

# Called by

- [AssembleRelationMatrix](../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/AssembleRelationMatrix.md)