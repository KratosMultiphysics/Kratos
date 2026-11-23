---
type: C++ Method
title: WriteTagsToJson
resource: utilities/assign_unique_model_part_collection_tag_utility.cpp#L258-L282
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_mpi_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Parameters AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson( const std::string& rFilename, const IndexStringMapType& rCollections )`

# Calls

- [AddEmptyArray](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyArray.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_mpi_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)