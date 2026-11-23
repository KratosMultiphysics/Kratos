---
type: C++ Function
title: GetWritingNDData
resource: input_output/vtu_output.cpp#L173-L219
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtu_output/AddConnectivityData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddFieldsUsingTensorAdaptorImpl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/for
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddFields
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`NDDataPointerType GetWritingNDData( const std::vector<IndexType>& rWritingIndices, NDDataPointerType pNDData)`

# Called by

- [AddConnectivityData](../../../functions/input_output/vtu_output/AddConnectivityData.md)
- [AddFieldsUsingTensorAdaptorImpl](../../../functions/input_output/vtu_output/AddFieldsUsingTensorAdaptorImpl.md)
- [for](../../../functions/input_output/vtu_output/for.md)
- [AddFields](../../../functions/input_output/vtu_output/AddFields.md)