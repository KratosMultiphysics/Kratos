---
type: C++ Function
title: AddKratosApplicationToPython
resource: python/add_kratos_application_to_python.cpp#L41-L53
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/add_kratos_application_to_python/RegisterToPythonApplicationVariables
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddKratosApplicationToPython(pybind11::module& m)`

# Calls

- [RegisterToPythonApplicationVariables](../../../functions/python/add_kratos_application_to_python/RegisterToPythonApplicationVariables.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)