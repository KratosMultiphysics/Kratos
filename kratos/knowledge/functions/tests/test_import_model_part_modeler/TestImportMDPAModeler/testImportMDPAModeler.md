---
type: Python Method
title: testImportMDPAModeler
resource: tests/test_import_model_part_modeler.py#L23-L189
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfTables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def testImportMDPAModeler(self): # Set up the import model part modeler`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [NumberOfTables](../../../../functions/includes/model_part/NumberOfTables.md)
- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)