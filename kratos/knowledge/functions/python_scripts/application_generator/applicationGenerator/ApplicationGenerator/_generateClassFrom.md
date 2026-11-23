---
type: Python Method
title: _generateClassFrom
resource: python_scripts/application_generator/applicationGenerator.py#L236-L270
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/exists
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/templateRule/TemplateRule/GetRule
    resolved_by: tree-sitter
    confidence: exact
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

`def _generateClassFrom(self, entityType, fromContainer): # Generate elements`

# Calls

- [exists](../../../../../functions/sources/kratos_filesystem/exists.md)
- [GetRule](../../../../../functions/python_scripts/application_generator/utils/templateRule/TemplateRule/GetRule.md)
- [remove](../../../../../functions/sources/kratos_filesystem/remove.md)
- [_applyTemplateRulesToFile](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_applyTemplateRulesToFile.md)

# Called by

- [Generate](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate.md)