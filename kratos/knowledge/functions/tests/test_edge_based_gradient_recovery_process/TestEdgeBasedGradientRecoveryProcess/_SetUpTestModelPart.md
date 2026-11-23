---
type: Python Method
title: _SetUpTestModelPart
resource: tests/test_edge_based_gradient_recovery_process.py#L15-L36
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/testScalarHistoricalAndHistorical
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/testScalarNonHistoricalAndNonHistorical
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/testVectorHistoricalAndNonHistorical
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _SetUpTestModelPart(self, OriginVar = None, GradientVar = None): # Create the background mesh`

# Calls

- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetInt](../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)

# Called by

- [testScalarHistoricalAndHistorical](../../../../functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/testScalarHistoricalAndHistorical.md)
- [testScalarNonHistoricalAndNonHistorical](../../../../functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/testScalarNonHistoricalAndNonHistorical.md)
- [testVectorHistoricalAndNonHistorical](../../../../functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/testVectorHistoricalAndNonHistorical.md)