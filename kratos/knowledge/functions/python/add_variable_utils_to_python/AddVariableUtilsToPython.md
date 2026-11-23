---
type: C++ Function
title: AddVariableUtilsToPython
resource: python/add_variable_utils_to_python.cpp#L215-L560
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/VariableUtils/SelectNodeList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddVariableUtilsToPython(pybind11::module &m)`

# Calls

- [SelectNodeList](../../../functions/utilities/variable_utils/VariableUtils/SelectNodeList.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)