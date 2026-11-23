---
type: C++ Function
title: AssembleRelationMatrix
resource: solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp#L212-L397
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/ProcessMultifreedomConstraint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/ProcessMasterSlaveConstraint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssembleRelationMatrix(const typename ConstraintAssembler<TSparse,TDense>::ConstraintArray& rConstraints, const ProcessInfo& rProcessInfo, typename TSparse::MatrixType& rRelationMatrix, typename TSparse::MatrixType& rHessian, typename TSparse::VectorType& rConstraintGaps, CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)`

# Calls

- [ProcessMultifreedomConstraint](../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/ProcessMultifreedomConstraint.md)
- [GetData](../../../../../functions/sources/gid_io/GidIOBase/GetData.md)
- [ProcessMasterSlaveConstraint](../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/ProcessMasterSlaveConstraint.md)