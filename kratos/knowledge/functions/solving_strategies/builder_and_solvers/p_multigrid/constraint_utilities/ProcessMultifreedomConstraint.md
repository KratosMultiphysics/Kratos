---
type: C++ Function
title: ProcessMultifreedomConstraint
resource: solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities.hpp#L102-L118
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/AssembleRelationMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void ProcessMultifreedomConstraint(std::vector<std::size_t>& rConstraintIndices, std::vector<std::size_t>& rDofIds, const MasterSlaveConstraint& rConstraint, const std::vector<std::size_t>& rMasterDofIds, const CSRHashMap<std::size_t,std::pair<std::size_t,std::size_t>>& rConstraintIdMap)`

# Calls

- [GetData](../../../../../functions/sources/gid_io/GidIOBase/GetData.md)

# Called by

- [AssembleRelationMatrix](../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/AssembleRelationMatrix.md)