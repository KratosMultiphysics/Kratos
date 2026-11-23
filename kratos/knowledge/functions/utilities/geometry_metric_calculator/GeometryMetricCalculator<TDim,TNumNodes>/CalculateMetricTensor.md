---
type: C++ Method
title: CalculateMetricTensor
resource: utilities/geometry_metric_calculator.cpp#L88-L114
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/GaussSeidelEigenSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometryMetricCalculator<TDim,TNumNodes>::CalculateMetricTensor( const GeometryType& rGeometry, BoundedMatrix<double, TDim, TDim>& rMetricTensor, double& rReferenceElementSize, double& rMetricInfimum, double& rMetricSupremum)`

# Calls

- [GaussSeidelEigenSystem](../../../../functions/utilities/math_utils/GaussSeidelEigenSystem.md)