---
type: C++ Function
title: AddContainersToPython
resource: future/python/add_containers_to_python.cpp#L34-L85
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/containers/linear_system/LinearSystem/pSetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system/LinearSystem/pSetLinearOperator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddContainersToPython(py::module& m)`

# Calls

- [pSetVector](../../../../functions/future/containers/linear_system/LinearSystem/pSetVector.md)
- [pSetLinearOperator](../../../../functions/future/containers/linear_system/LinearSystem/pSetLinearOperator.md)
- [pGetVector](../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)