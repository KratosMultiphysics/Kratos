---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/utilities/test_mpi_assign_unique_model_part_collection_tag_utility.cpp#L39-L112
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
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ReadTagsFromJson
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(AssignMPIUniqueModelPartCollectionTagUtility, KratosMPICoreFastSuite)`

# Calls

- [Create2DGeometry](../../../../../../functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry.md)
- [ComputeTags](../../../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ComputeTags.md)
- [WriteTagsToJson](../../../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson.md)
- [ReadTagsFromJson](../../../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ReadTagsFromJson.md)