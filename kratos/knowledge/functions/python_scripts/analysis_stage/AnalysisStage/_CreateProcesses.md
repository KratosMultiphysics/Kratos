---
type: Python Method
title: _CreateProcesses
resource: python_scripts/analysis_stage.py#L392-L427
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/process_factory/KratosProcessFactory/ConstructListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/__CreateListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _CreateProcesses(self, parameter_name, initialization_order):`

# Calls

- [ConstructListOfProcesses](../../../../functions/python_scripts/process_factory/KratosProcessFactory/ConstructListOfProcesses.md)

# Called by

- [__CreateListOfProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/__CreateListOfProcesses.md)