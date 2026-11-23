---
type: Python Method
title: __init__
resource: python_scripts/unv_output_process.py#L14-L51
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/unv_output/UnvOutput/UnvOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/InitializeMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [UnvOutput](../../../../functions/input_output/unv_output/UnvOutput/UnvOutput.md)
- [InitializeMesh](../../../../functions/includes/gid_io/InitializeMesh.md)