---
type: Python Method
title: _ModelersSetupGeometryModel
resource: python_scripts/analysis_stage.py#L315-L322
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetSimulationName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _ModelersSetupGeometryModel(self): # Import or generate geometry models from external input.`

# Calls

- [_GetListOfModelers](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfModelers.md)
- [_GetSimulationName](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSimulationName.md)

# Called by

- [Initialize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Initialize.md)