---
type: Python Method
title: Generate
resource: python_scripts/application_generator/applicationGenerator.py#L130-L213
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/rename
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_addClassToSources
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_applyTemplateRulesToFile
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

`def Generate(self):`

# Calls

- [rename](../../../../../functions/sources/kratos_filesystem/rename.md)
- [_addClassToSources](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_addClassToSources.md)
- [remove](../../../../../functions/sources/kratos_filesystem/remove.md)
- [_applyTemplateRulesToFile](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_applyTemplateRulesToFile.md)
- [_generateClassFrom](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom.md)
- [_generateVariablesFrom](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateVariablesFrom.md)