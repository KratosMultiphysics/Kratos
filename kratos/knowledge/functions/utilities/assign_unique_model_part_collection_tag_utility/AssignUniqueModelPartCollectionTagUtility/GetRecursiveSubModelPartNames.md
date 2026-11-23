---
type: C++ Method
title: GetRecursiveSubModelPartNames
resource: utilities/assign_unique_model_part_collection_tag_utility.cpp#L288-L311
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/GetSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ComputeTags
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames( ModelPart& rThisModelPart, std::string Prefix )`

# Calls

- [GetSubModelPartNames](../../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [ComputeTags](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ComputeTags.md)