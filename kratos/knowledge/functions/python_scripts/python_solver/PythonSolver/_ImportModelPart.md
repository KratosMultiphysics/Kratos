---
type: Python Method
title: _ImportModelPart
resource: python_scripts/python_solver.py#L148-L167
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/single_import_model_part/SingleImportModelPart/Import
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_solver/PythonSolver/_GetRestartSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/LoadRestart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _ImportModelPart(self, model_part, model_part_import_settings):`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [Import](../../../../functions/utilities/single_import_model_part/SingleImportModelPart/Import.md)
- [_GetRestartSettings](../../../../functions/python_scripts/python_solver/PythonSolver/_GetRestartSettings.md)
- [LoadRestart](../../../../functions/python_scripts/restart_utility/RestartUtility/LoadRestart.md)