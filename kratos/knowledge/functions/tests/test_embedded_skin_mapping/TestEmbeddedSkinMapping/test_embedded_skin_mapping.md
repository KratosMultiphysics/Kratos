---
type: Python Method
title: test_embedded_skin_mapping
resource: tests/test_embedded_skin_mapping.py#L10-L79
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_embedded_skin_mapping(self):`

# Calls

- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetInt](../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [GenerateSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin.md)