---
type: C++ Method
title: pSetLinearOperator
resource: future/containers/linear_system.h#L294-L300
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
---

# Signature

`void pSetLinearOperator( typename LinearOperatorType::UniquePointer&& pLinearOperator, SparseMatrixTag Tag)`

# Called by

- [LinearSystem](../../../../../functions/future/containers/linear_system/LinearSystem/LinearSystem.md)
- [AddContainersToPython](../../../../../functions/future/python/add_containers_to_python/AddContainersToPython.md)