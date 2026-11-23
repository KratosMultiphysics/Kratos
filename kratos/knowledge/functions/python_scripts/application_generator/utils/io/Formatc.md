---
type: Python Function
title: Formatc
resource: python_scripts/application_generator/utils/io.py#L23-L33
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/classMemberCreator/ClassMemberCreator/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/templateRule/TemplateRule/GetRule
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Formatc(stringList, where):`

# Called by

- [__init__](../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/__init__.md)
- [__init__](../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/__init__.md)
- [__init__](../../../../../functions/python_scripts/application_generator/classes/classMemberCreator/ClassMemberCreator/__init__.md)
- [GetRule](../../../../../functions/python_scripts/application_generator/utils/templateRule/TemplateRule/GetRule.md)