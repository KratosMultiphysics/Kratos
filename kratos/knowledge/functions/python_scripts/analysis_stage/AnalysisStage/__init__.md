---
type: Python Method
title: __init__
resource: python_scripts/analysis_stage.py#L12-L43
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_AddVariables
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, project_parameters):`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [_AddVariables](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_AddVariables.md)