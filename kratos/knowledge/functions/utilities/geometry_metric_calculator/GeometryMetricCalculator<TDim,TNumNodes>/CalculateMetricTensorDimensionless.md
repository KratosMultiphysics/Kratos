---
type: C++ Method
title: CalculateMetricTensorDimensionless
resource: utilities/geometry_metric_calculator.cpp#L117-L134
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_geometry_metric_calculator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometryMetricCalculator<TDim,TNumNodes>::CalculateMetricTensorDimensionless( const GeometryType& rGeometry, BoundedMatrix<double, TDim, TDim>& rMetricTensorDimensionless, double& rReferenceElementSize, double& rMetricInfimum, double& rMetricSupremum)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_geometry_metric_calculator/KRATOS_TEST_CASE_IN_SUITE.md)