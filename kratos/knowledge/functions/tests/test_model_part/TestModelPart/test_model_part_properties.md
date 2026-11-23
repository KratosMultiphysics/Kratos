---
type: Python Method
title: test_model_part_properties
resource: tests/test_model_part.py#L294-L380
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/AddSubProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/NumberOfSubproperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/HasSubProperties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_model_part_properties(self):`

# Calls

- [RecursivelyHasProperties](../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [AddSubProperties](../../../../functions/includes/properties/AddSubProperties.md)
- [NumberOfSubproperties](../../../../functions/includes/properties/NumberOfSubproperties.md)
- [HasSubProperties](../../../../functions/includes/properties/HasSubProperties.md)