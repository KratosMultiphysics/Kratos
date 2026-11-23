---
type: C++ Function
title: AddModelPartToPython
resource: python/add_model_part_to_python.cpp#L742-L998
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_model_part_to_python/GetNonHistoricalVariablesNames
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddModelPartToPython(pybind11::module& m)`

# Calls

- [RecursivelyHasProperties](../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [CreateNewProperties](../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [GetModel](../../../functions/python_scripts/project/Project/GetModel.md)
- [GetNonHistoricalVariablesNames](../../../functions/python/add_model_part_to_python/GetNonHistoricalVariablesNames.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)