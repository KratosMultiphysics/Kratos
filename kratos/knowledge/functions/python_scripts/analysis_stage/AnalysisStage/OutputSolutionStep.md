---
type: Python Method
title: OutputSolutionStep
resource: python_scripts/analysis_stage.py#L173-L188
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/RunSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def OutputSolutionStep(self):`

# Calls

- [_GetListOfOutputProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfOutputProcesses.md)
- [_GetListOfProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfProcesses.md)

# Called by

- [RunSolutionLoop](../../../../functions/python_scripts/analysis_stage/AnalysisStage/RunSolutionLoop.md)