---
type: Python Method
title: CreateStage
resource: python_scripts/orchestrators/orchestrator.py#L49-L120
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CreateStage(self, stage_name : str) -> AnalysisStage:`

# Calls

- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetModel](../../../../../functions/python_scripts/project/Project/GetModel.md)

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)