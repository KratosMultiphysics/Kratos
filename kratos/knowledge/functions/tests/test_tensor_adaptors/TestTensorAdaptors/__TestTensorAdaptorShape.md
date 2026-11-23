---
type: Python Method
title: __TestTensorAdaptorShape
resource: tests/test_tensor_adaptors.py#L1022-L1037
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeHistoricalVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorProperties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __TestTensorAdaptorShape(self, container, TensorAdaptorType):`

# Calls

- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)

# Called by

- [test_ShapeHistoricalVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeHistoricalVariableTensorAdaptor.md)
- [test_ShapeVariableTensorAdaptorNodes](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorNodes.md)
- [test_ShapeVariableTensorAdaptorConditions](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorConditions.md)
- [test_ShapeVariableTensorAdaptorElements](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorElements.md)
- [test_ShapeVariableTensorAdaptorProperties](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_ShapeVariableTensorAdaptorProperties.md)