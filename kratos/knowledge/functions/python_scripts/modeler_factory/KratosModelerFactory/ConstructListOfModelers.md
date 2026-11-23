---
type: Python Method
title: ConstructListOfModelers
resource: python_scripts/modeler_factory.py#L9-L35
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/IssueDeprecationWarning
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_CreateModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_cad_tessellation_modeler/run_modelers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ConstructListOfModelers( self, model, modeler_list ):`

# Calls

- [IssueDeprecationWarning](../../../../functions/python_scripts/kratos_utilities/IssueDeprecationWarning.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [_CreateModelers](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateModelers.md)
- [run_modelers](../../../../functions/tests/test_cad_tessellation_modeler/run_modelers.md)