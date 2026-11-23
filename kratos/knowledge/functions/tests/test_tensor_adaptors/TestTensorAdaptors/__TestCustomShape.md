---
type: Python Method
title: __TestCustomShape
resource: tests/test_tensor_adaptors.py#L1039-L1055
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomHistoricalVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorProperties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __TestCustomShape(self, container, TensorAdaptorType):`

# Calls

- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)

# Called by

- [test_CustomHistoricalVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomHistoricalVariableTensorAdaptor.md)
- [test_CustomVariableTensorAdaptorNodes](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorNodes.md)
- [test_CustomVariableTensorAdaptorConditions](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorConditions.md)
- [test_CustomVariableTensorAdaptorElements](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorElements.md)
- [test_CustomVariableTensorAdaptorProperties](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CustomVariableTensorAdaptorProperties.md)