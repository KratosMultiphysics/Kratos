---
type: C++ Method
title: pGetStorage
resource: tensor_adaptors/tensor_adaptor.cpp#L150-L153
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtu_output/AddFieldsUsingTensorAdaptorImpl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/for
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddFields
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_tensor_adaptors_to_python/AddBaseTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename TensorAdaptor<TDataType>::Storage::Pointer TensorAdaptor<TDataType>::pGetStorage()`

# Called by

- [AddFieldsUsingTensorAdaptorImpl](../../../../functions/input_output/vtu_output/AddFieldsUsingTensorAdaptorImpl.md)
- [for](../../../../functions/input_output/vtu_output/for.md)
- [AddFields](../../../../functions/input_output/vtu_output/AddFields.md)
- [WriteUnstructuredGridData](../../../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)
- [AddBaseTensorAdaptor](../../../../functions/python/add_tensor_adaptors_to_python/AddBaseTensorAdaptor.md)