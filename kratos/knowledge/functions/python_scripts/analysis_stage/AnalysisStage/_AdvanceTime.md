---
type: Python Method
title: _AdvanceTime
resource: python_scripts/analysis_stage.py#L281-L285
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_solver/PythonSolver/AdvanceInTime
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/RunSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _AdvanceTime(self):`

# Calls

- [_GetSolver](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSolver.md)
- [AdvanceInTime](../../../../functions/python_scripts/python_solver/PythonSolver/AdvanceInTime.md)

# Called by

- [RunSolutionLoop](../../../../functions/python_scripts/analysis_stage/AnalysisStage/RunSolutionLoop.md)