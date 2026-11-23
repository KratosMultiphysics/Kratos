---
type: Python Method
title: _AddDofs
resource: python_scripts/analysis_stage.py#L301-L305
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetSolver
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _AddDofs(self):`

# Calls

- [_GetSolver](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSolver.md)

# Called by

- [Initialize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Initialize.md)