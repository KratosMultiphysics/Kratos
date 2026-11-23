---
type: C++ Method
title: RenameModelPart
resource: containers/model.cpp#L87-L104
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Model::RenameModelPart( const std::string& OldName, const std::string& NewName )`

# Calls

- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE.md)