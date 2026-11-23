---
type: C++ Function
title: AddSearchStrategiesToPython
resource: python/add_search_strategies_to_python.cpp#L149-L609
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/add_search_strategies_to_python/GenerateSpatialSearchSolutionTuple
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddSearchStrategiesToPython(pybind11::module& m)`

# Calls

- [GenerateSpatialSearchSolutionTuple](../../../functions/python/add_search_strategies_to_python/GenerateSpatialSearchSolutionTuple.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)