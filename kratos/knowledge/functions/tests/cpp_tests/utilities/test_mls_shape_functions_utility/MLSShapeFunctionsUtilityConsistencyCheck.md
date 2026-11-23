---
type: C++ Function
title: MLSShapeFunctionsUtilityConsistencyCheck
resource: tests/cpp_tests/utilities/test_mls_shape_functions_utility.cpp#L51-L144
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_mls_shape_functions_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MLSShapeFunctionsUtilityConsistencyCheck( std::function<double(array_1d<double,3>)>& rAnalyticalSolutionFunction, std::function<array_1d<double,2>(array_1d<double,3>)>& rAnalyticalSolutionGradientFunction, std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&, Matrix&)>& rMLSShapeFunctionsAndGradientsFunction)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_mls_shape_functions_utility/KRATOS_TEST_CASE_IN_SUITE.md)