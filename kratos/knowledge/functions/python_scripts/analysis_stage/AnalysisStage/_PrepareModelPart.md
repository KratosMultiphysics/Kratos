---
type: Python Method
title: _PrepareModelPart
resource: python_scripts/analysis_stage.py#L293-L299
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_solver/PythonSolver/PrepareModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _PrepareModelPart(self):`

# Calls

- [_GetSolver](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSolver.md)
- [PrepareModelPart](../../../../functions/python_scripts/python_solver/PythonSolver/PrepareModelPart.md)

# Called by

- [Initialize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Initialize.md)