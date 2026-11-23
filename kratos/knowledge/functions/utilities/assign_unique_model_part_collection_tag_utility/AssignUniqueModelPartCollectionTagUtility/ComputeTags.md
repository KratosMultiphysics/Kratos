---
type: C++ Method
title: ComputeTags
resource: utilities/assign_unique_model_part_collection_tag_utility.cpp#L40-L170
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/GetRecursiveSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/SetParallelModelPartAndSubModelPartCollectionsAndCombinations
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_mpi_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/DebugAssignUniqueModelPartCollectionTag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignUniqueModelPartCollectionTagUtility::ComputeTags( IndexIndexMapType& rNodeTags, IndexIndexMapType& rCondTags, IndexIndexMapType& rElemTags, IndexStringMapType& rCollections )`

# Calls

- [GetRecursiveSubModelPartNames](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/GetRecursiveSubModelPartNames.md)
- [SetParallelModelPartAndSubModelPartCollectionsAndCombinations](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/SetParallelModelPartAndSubModelPartCollectionsAndCombinations.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_mpi_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [DebugAssignUniqueModelPartCollectionTag](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/DebugAssignUniqueModelPartCollectionTag.md)