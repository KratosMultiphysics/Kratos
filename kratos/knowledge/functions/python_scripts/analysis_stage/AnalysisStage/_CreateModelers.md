---
type: Python Method
title: _CreateModelers
resource: python_scripts/analysis_stage.py#L350-L371
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/__BackwardCompatibleModelersCreation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modeler_factory/KratosModelerFactory/ConstructListOfModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructListOfItems
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _CreateModelers(self):`

# Calls

- [__BackwardCompatibleModelersCreation](../../../../functions/python_scripts/analysis_stage/AnalysisStage/__BackwardCompatibleModelersCreation.md)
- [ConstructListOfModelers](../../../../functions/python_scripts/modeler_factory/KratosModelerFactory/ConstructListOfModelers.md)
- [ConstructListOfItems](../../../../functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructListOfItems.md)

# Called by

- [Initialize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Initialize.md)