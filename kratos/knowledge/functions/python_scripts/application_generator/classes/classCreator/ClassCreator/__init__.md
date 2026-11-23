---
type: Python Method
title: __init__
resource: python_scripts/application_generator/classes/classCreator.py#L16-L94
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/application_generator/utils/io/TestCamel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/io/Formatc
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/io/ToUpperFromCamel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/io/ToLowerFromCamel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_CheckClassNameAvailable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__( self, name, base=None, template=None, members=None, procedures=None, author='KratosAppGenerator' ): # Check that name is in camelcase, or at least something parseable`

# Calls

- [TestCamel](../../../../../../functions/python_scripts/application_generator/utils/io/TestCamel.md)
- [Formatc](../../../../../../functions/python_scripts/application_generator/utils/io/Formatc.md)
- [ToUpperFromCamel](../../../../../../functions/python_scripts/application_generator/utils/io/ToUpperFromCamel.md)
- [ToLowerFromCamel](../../../../../../functions/python_scripts/application_generator/utils/io/ToLowerFromCamel.md)
- [_CheckClassNameAvailable](../../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_CheckClassNameAvailable.md)