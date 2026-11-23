---
type: Python Method
title: TranslateLegacyVariablesAccordingToCurrentStandard
resource: python_scripts/unv_output_process.py#L54-L62
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

`def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings: KratosMultiphysics.Parameters) -> None:`

# Calls

- [HasDeprecatedVariable](../../../../functions/python_scripts/deprecation_management/DeprecationManager/HasDeprecatedVariable.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)