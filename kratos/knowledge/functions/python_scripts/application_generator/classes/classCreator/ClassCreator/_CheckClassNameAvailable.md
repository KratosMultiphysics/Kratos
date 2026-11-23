---
type: Python Method
title: _CheckClassNameAvailable
resource: python_scripts/application_generator/classes/classCreator.py#L205-L258
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_CheckClassNameInDir
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/io/GetApplicationsDirectory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/io/GetKratosDirectory
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _CheckClassNameAvailable( self, name, context="Application", searchSimilar=False ):`

# Calls

- [_CheckClassNameInDir](../../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_CheckClassNameInDir.md)
- [GetApplicationsDirectory](../../../../../../functions/python_scripts/application_generator/utils/io/GetApplicationsDirectory.md)
- [GetKratosDirectory](../../../../../../functions/python_scripts/application_generator/utils/io/GetKratosDirectory.md)

# Called by

- [__init__](../../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/__init__.md)