---
type: C++ Function
title: AddMatrixToPython
resource: python/add_matrix_to_python.cpp#L68-L226
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddMatrixToPython(pybind11::module& m)`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)