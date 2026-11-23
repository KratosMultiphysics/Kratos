---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility.cpp#L40-L127
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ComputeTags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/GetSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/GetRecursiveSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ReadTagsFromJson
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(AssignUniqueModelPartCollectionTagUtility, KratosCoreFastSuite)`

# Calls

- [Create2DGeometry](../../../../../functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry.md)
- [ComputeTags](../../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ComputeTags.md)
- [GetSubModelPartNames](../../../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)
- [GetRecursiveSubModelPartNames](../../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/GetRecursiveSubModelPartNames.md)
- [WriteTagsToJson](../../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson.md)
- [ReadTagsFromJson](../../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ReadTagsFromJson.md)
- [GetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)