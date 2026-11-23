---
type: C++ Method
title: CalculateSensitivity
resource: utilities/line_sensitivity_utility.cpp#L26-L52
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LineSensitivityUtility::CalculateSensitivity( const ShapeParameter Deriv, double& rIntegrationWeightSensitivity) const`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)