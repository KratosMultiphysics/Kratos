---
type: Python Function
title: IssueDeprecationWarning
resource: python_scripts/kratos_utilities.py#L5-L6
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/ModifyAfterSolverInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_CheckDeprecatedOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modeler_factory/KratosModelerFactory/ConstructListOfModelers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def IssueDeprecationWarning(label, *message):`

# Called by

- [ModifyAfterSolverInitialize](../../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyAfterSolverInitialize.md)
- [_CheckDeprecatedOutputProcesses](../../../functions/python_scripts/analysis_stage/AnalysisStage/_CheckDeprecatedOutputProcesses.md)
- [ConstructListOfModelers](../../../functions/python_scripts/modeler_factory/KratosModelerFactory/ConstructListOfModelers.md)