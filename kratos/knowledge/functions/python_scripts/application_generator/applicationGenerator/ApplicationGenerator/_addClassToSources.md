---
type: Python Method
title: _addClassToSources
resource: python_scripts/application_generator/applicationGenerator.py#L228-L234
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/application_generator/utils/templateRule/TemplateRule/GetRule
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _addClassToSources(self, entityType, fromContainer):`

# Calls

- [GetRule](../../../../../functions/python_scripts/application_generator/utils/templateRule/TemplateRule/GetRule.md)

# Called by

- [Generate](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate.md)