---
type: C++ Method
title: GetInitialPositionsVector
resource: utilities/variable_utils.cpp#L410-L425
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_variable_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_solutionsteps_vector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`TVectorType VariableUtils::GetInitialPositionsVector( const ModelPart::NodesContainerType& rNodes, const unsigned int Dimension)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_variable_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_solutionsteps_vector](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_solutionsteps_vector.md)