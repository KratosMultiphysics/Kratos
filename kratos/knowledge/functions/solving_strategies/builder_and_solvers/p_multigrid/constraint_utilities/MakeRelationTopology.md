---
type: C++ Function
title: MakeRelationTopology
resource: solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp#L122-L208
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MakeRelationTopology(std::size_t SystemSize, typename ConstraintAssembler<TSparse,TDense>::ConstraintArray::const_iterator itConstraintBegin, typename ConstraintAssembler<TSparse,TDense>::ConstraintArray::const_iterator itConstraintEnd, const ProcessInfo& rProcessInfo, typename TSparse::MatrixType& rRelationMatrix, typename TSparse::VectorType& rConstraintGaps, CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)`

# Calls

- [GetData](../../../../../functions/sources/gid_io/GidIOBase/GetData.md)