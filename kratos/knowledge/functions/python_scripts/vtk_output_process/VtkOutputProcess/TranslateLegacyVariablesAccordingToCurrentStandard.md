---
type: Python Method
title: TranslateLegacyVariablesAccordingToCurrentStandard
resource: python_scripts/vtk_output_process.py#L41-L61
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/deprecation_management/DeprecationManager/HasDeprecatedVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings: KratosMultiphysics.Parameters) -> None: # Defining a string to help the user understand where the warnings come from (in case any is thrown)`

# Calls

- [HasDeprecatedVariable](../../../../functions/python_scripts/deprecation_management/DeprecationManager/HasDeprecatedVariable.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)