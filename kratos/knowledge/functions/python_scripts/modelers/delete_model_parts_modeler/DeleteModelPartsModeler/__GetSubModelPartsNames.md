---
type: Python Method
title: __GetSubModelPartsNames
resource: python_scripts/modelers/delete_model_parts_modeler.py#L63-L68
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetSubModelPartsNames(self, model_part):`

# Calls

- [FullName](../../../../../functions/includes/model_part/FullName.md)

# Called by

- [SetupModelPart](../../../../../functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/SetupModelPart.md)