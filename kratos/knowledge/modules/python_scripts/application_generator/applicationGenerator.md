---
type: Python Module
title: applicationGenerator
resource: python_scripts/application_generator/applicationGenerator.py#L1-L344
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sys
    resolved_by: tree-sitter
    confidence: exact
  - target: external/shutil
    resolved_by: tree-sitter
    confidence: exact
  - target: external/difflib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/classes-elementcreator
    resolved_by: tree-sitter
    confidence: exact
  - target: external/classes-conditioncreator
    resolved_by: tree-sitter
    confidence: exact
  - target: external/classes-processcreator
    resolved_by: tree-sitter
    confidence: exact
  - target: external/classes-variablecreator
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utils-io
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utils-constants
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utils-templaterule
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ApplicationGenerator](../../../classes/python_scripts/application_generator/applicationGenerator/ApplicationGenerator.md)
- [__init__](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/__init__.md)
- [GenerateFile](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/GenerateFile.md)
- [AddElements](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/AddElements.md)
- [AddConditions](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/AddConditions.md)
- [AddProcesses](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/AddProcesses.md)
- [AddVariables](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/AddVariables.md)
- [AddCustomStrategies](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/AddCustomStrategies.md)
- [AddCustomVariables](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/AddCustomVariables.md)
- [AddCustomConstitutiveLaws](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/AddCustomConstitutiveLaws.md)
- [Generate](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate.md)
- [_applyTemplateRulesToFile](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_applyTemplateRulesToFile.md)
- [_addClassToSources](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_addClassToSources.md)
- [_generateClassFrom](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom.md)
- [_generateVariablesFrom](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateVariablesFrom.md)
- [_addGenericClass](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_addGenericClass.md)

# Imports

- `os`
- `sys`
- `shutil`
- `difflib`
- `classes.elementCreator`
- `classes.conditionCreator`
- `classes.processCreator`
- `classes.variableCreator`
- `utils.io`
- `utils.constants`
- `utils.templateRule`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)