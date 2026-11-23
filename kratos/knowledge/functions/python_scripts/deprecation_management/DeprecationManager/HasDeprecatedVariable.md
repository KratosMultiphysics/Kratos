---
type: Python Method
title: HasDeprecatedVariable
resource: python_scripts/deprecation_management.py#L20-L34
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/tikz_output_process/TikZOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/unv_output_process/UnvOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtk_output_process/VtkOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def HasDeprecatedVariable(context_string, parameters, old_variable_name, new_variable_name):`

# Called by

- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)
- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/tikz_output_process/TikZOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)
- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/unv_output_process/UnvOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)
- [TranslateLegacyVariablesAccordingToCurrentStandard](../../../../functions/python_scripts/vtk_output_process/VtkOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)