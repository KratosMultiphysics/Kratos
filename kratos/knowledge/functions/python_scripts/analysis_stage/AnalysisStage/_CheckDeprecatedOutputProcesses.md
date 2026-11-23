---
type: Python Method
title: _CheckDeprecatedOutputProcesses
resource: python_scripts/analysis_stage.py#L441-L449
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/IssueDeprecationWarning
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/__CreateListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _CheckDeprecatedOutputProcesses(self, list_of_processes):`

# Calls

- [IssueDeprecationWarning](../../../../functions/python_scripts/kratos_utilities/IssueDeprecationWarning.md)

# Called by

- [__CreateListOfProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/__CreateListOfProcesses.md)