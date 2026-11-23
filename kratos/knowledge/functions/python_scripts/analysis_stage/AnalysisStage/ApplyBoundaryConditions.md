---
type: Python Method
title: ApplyBoundaryConditions
resource: python_scripts/analysis_stage.py#L232-L236
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ApplyBoundaryConditions(self):`

# Calls

- [_GetListOfProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfProcesses.md)

# Called by

- [InitializeSolutionStep](../../../../functions/python_scripts/analysis_stage/AnalysisStage/InitializeSolutionStep.md)