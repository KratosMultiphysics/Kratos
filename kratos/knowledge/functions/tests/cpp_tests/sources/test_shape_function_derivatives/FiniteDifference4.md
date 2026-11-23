---
type: C++ Function
title: FiniteDifference4
resource: tests/cpp_tests/sources/test_shape_function_derivatives.cpp#L128-L132
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/check_scalar_base_process/aux_object_cpp_callback/f
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradient
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double FiniteDifference4(const std::function<double(double)>& f, double delta)`

# Calls

- [f](../../../../../functions/python_scripts/check_scalar_base_process/aux_object_cpp_callback/f.md)

# Called by

- [TestShapeFunctionsLocalGradient](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradient.md)