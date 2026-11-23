---
type: Python Method
title: _applyTemplateRulesToFile
resource: python_scripts/application_generator/applicationGenerator.py#L216-L226
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateVariablesFrom
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _applyTemplateRulesToFile(self, src, dst, rules):`

# Called by

- [Generate](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate.md)
- [_generateClassFrom](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom.md)
- [_generateVariablesFrom](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateVariablesFrom.md)