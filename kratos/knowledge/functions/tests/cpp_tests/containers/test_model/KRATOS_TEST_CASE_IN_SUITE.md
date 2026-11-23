---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/containers/test_model.cpp#L23-L46
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/RenameModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/GetModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ModelGetModelPart, KratosCoreFastSuite)`

# Calls

- [HasModelPart](../../../../../functions/containers/model/Model/HasModelPart.md)
- [DeleteModelPart](../../../../../functions/containers/model/Model/DeleteModelPart.md)
- [RenameModelPart](../../../../../functions/containers/model/Model/RenameModelPart.md)
- [GetModel](../../../../../functions/python_scripts/project/Project/GetModel.md)
- [GetModelPartNames](../../../../../functions/containers/model/Model/GetModelPartNames.md)
- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)