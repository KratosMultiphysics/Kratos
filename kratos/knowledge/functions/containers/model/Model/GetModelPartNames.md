---
type: C++ Method
title: GetModelPartNames
resource: containers/model.cpp#L204-L235
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/GetSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> Model::GetModelPartNames() const`

# Calls

- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [GetSubModelPartNames](../../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)

# Called by

- [Check](../../../../functions/processes/time_averaging_process/TimeAveragingProcess/Check.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE.md)