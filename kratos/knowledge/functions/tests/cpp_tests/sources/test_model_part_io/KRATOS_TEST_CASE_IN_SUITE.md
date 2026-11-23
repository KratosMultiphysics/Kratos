---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_model_part_io.cpp#L26-L173
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/SetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE( ModelPartIOSubModelPartsDivision, KratosCoreFastSuite)`

# Calls

- [SetNumberOfColors](../../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [NumberOfSubModelParts](../../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)
- [remove](../../../../../functions/sources/kratos_filesystem/remove.md)