---
type: C++ Function
title: CreateMatrixInterface
resource: python/add_matrix_to_python.cpp#L29-L66
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`py::class_< TMatrixType > CreateMatrixInterface(pybind11::module& m, std::string Name )`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)