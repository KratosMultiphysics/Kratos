---
type: Python Method
title: __init__
resource: python_scripts/base_convergence_criteria_factory.py#L6-L45
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, convergence_criterion_parameters): # Note that all the convergence settings are introduced via a Kratos parameters object.`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)