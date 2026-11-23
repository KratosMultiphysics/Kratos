---
type: C++ Method
title: pSetVector
resource: future/containers/linear_system.h#L266-L272
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/containers/linear_system/LinearSystem/LinearSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/python/add_containers_to_python/AddContainersToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void pSetVector( typename VectorType::Pointer pVector, DenseVectorTag Tag)`

# Called by

- [LinearSystem](../../../../../functions/future/containers/linear_system/LinearSystem/LinearSystem.md)
- [AddContainersToPython](../../../../../functions/future/python/add_containers_to_python/AddContainersToPython.md)
- [ApplyBlockBuildMasterSlaveConstraints](../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints.md)