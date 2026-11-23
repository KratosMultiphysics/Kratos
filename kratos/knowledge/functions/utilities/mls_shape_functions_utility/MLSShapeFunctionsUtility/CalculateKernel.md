---
type: C++ Method
title: CalculateKernel
resource: utilities/mls_shape_functions_utility.cpp#L33-L40
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_mls_shape_functions_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mls_shape_functions_utility/MLSShapeFunctionsUtility/CalculateShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double MLSShapeFunctionsUtility::CalculateKernel( const array_1d<double,3>& rRadVect, const double h)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_mls_shape_functions_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateShapeFunctions](../../../../functions/utilities/mls_shape_functions_utility/MLSShapeFunctionsUtility/CalculateShapeFunctions.md)