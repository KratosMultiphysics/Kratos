---
type: Python Method
title: __CreateListOfProcesses
resource: python_scripts/analysis_stage.py#L456-L465
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetOrderOfProcessesInitialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_CreateProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_CheckDeprecatedOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetOrderOfOutputProcessesInitialization
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CreateListOfProcesses(self):`

# Calls

- [_GetOrderOfProcessesInitialization](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetOrderOfProcessesInitialization.md)
- [_CreateProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateProcesses.md)
- [_CheckDeprecatedOutputProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_CheckDeprecatedOutputProcesses.md)
- [_GetOrderOfOutputProcessesInitialization](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetOrderOfOutputProcessesInitialization.md)

# Called by

- [Initialize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Initialize.md)