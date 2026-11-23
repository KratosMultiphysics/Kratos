---
type: Python Method
title: PrintAnalysisStageProgressInformation
resource: python_scripts/analysis_stage.py#L154-L156
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetSimulationName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def PrintAnalysisStageProgressInformation(self):`

# Calls

- [_GetSimulationName](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSimulationName.md)

# Called by

- [InitializeSolutionStep](../../../../functions/python_scripts/analysis_stage/AnalysisStage/InitializeSolutionStep.md)