---
type: Python Method
title: GetRule
resource: python_scripts/application_generator/utils/templateRule.py#L10-L21
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/application_generator/utils/io/Formatc
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_addClassToSources
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/AddClassMemberVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/conditionCreator/ConditionCreator/AddFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/conditionCreator/ConditionCreator/AddDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/elementCreator/ElementCreator/AddFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/elementCreator/ElementCreator/AddDofs
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetRule(self, ruleName):`

# Calls

- [Formatc](../../../../../../functions/python_scripts/application_generator/utils/io/Formatc.md)

# Called by

- [_addClassToSources](../../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_addClassToSources.md)
- [_generateClassFrom](../../../../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom.md)
- [AddClassMemberVariables](../../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/AddClassMemberVariables.md)
- [AddFlags](../../../../../../functions/python_scripts/application_generator/classes/conditionCreator/ConditionCreator/AddFlags.md)
- [AddDofs](../../../../../../functions/python_scripts/application_generator/classes/conditionCreator/ConditionCreator/AddDofs.md)
- [AddFlags](../../../../../../functions/python_scripts/application_generator/classes/elementCreator/ElementCreator/AddFlags.md)
- [AddDofs](../../../../../../functions/python_scripts/application_generator/classes/elementCreator/ElementCreator/AddDofs.md)