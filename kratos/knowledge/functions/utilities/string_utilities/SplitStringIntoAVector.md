---
type: C++ Function
title: SplitStringIntoAVector
resource: utilities/string_utilities.cpp#L245-L257
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> SplitStringIntoAVector(const std::string& rText)`

# Called by

- [FillNodalConnectivitiesFromConstraintBlockInList](../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList.md)
- [ReadConstraintsBlock](../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock.md)
- [ReadConstraintsConnectivitiesBlock](../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock.md)
- [FillNodalConnectivitiesFromConstraintBlock](../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlock.md)
- [DivideConstraintsBlock](../../../functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE.md)