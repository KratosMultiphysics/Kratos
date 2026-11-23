---
type: Python Method
title: _CheckClassNameInDir
resource: python_scripts/application_generator/classes/classCreator.py#L169-L190
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_GenerateCandidateFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_IsClassDefinition
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_CheckClassNameAvailable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _CheckClassNameInDir( self, where, name, appRep, appSim, repFiles, simFiles ):`

# Calls

- [_GenerateCandidateFiles](../../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_GenerateCandidateFiles.md)
- [_IsClassDefinition](../../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_IsClassDefinition.md)

# Called by

- [_CheckClassNameAvailable](../../../../../../functions/python_scripts/application_generator/classes/classCreator/ClassCreator/_CheckClassNameAvailable.md)