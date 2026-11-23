---
type: Python Method
title: _generateVariablesFrom
resource: python_scripts/application_generator/applicationGenerator.py#L272-L335
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_applyTemplateRulesToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _generateVariablesFrom(self, entityType, fromContainer):`

# Calls

- [remove](../../../../../functions/sources/kratos_filesystem/remove.md)
- [_applyTemplateRulesToFile](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_applyTemplateRulesToFile.md)

# Called by

- [Generate](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate.md)