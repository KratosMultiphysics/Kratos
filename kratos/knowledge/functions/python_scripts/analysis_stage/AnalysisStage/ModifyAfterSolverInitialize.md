---
type: Python Method
title: ModifyAfterSolverInitialize
resource: python_scripts/analysis_stage.py#L228-L230
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/IssueDeprecationWarning
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ModifyAfterSolverInitialize(self):`

# Calls

- [IssueDeprecationWarning](../../../../functions/python_scripts/kratos_utilities/IssueDeprecationWarning.md)

# Called by

- [Initialize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Initialize.md)