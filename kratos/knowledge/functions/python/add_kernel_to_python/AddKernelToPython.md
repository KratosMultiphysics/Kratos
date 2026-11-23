---
type: C++ Function
title: AddKernelToPython
resource: python/add_kernel_to_python.cpp#L96-L172
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kernel/Kernel/SetPythonVersion
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddKernelToPython(pybind11::module& m)`

# Calls

- [SetPythonVersion](../../../functions/sources/kernel/Kernel/SetPythonVersion.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)