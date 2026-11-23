---
type: C++ Method
title: GetSubModelPartNames
resource: sources/model_part.cpp#L2161-L2171
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/model/Model/GetModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/GetRecursiveSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> ModelPart::GetSubModelPartNames() const`

# Calls

- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)

# Called by

- [GetModelPartNames](../../../../functions/containers/model/Model/GetModelPartNames.md)
- [RemoveSubModelPart](../../../../functions/sources/model_part/ModelPart/RemoveSubModelPart.md)
- [ErrorNonExistingSubModelPart](../../../../functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart.md)
- [WriteSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/WriteSubModelPartBlock.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [GetRecursiveSubModelPartNames](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/GetRecursiveSubModelPartNames.md)