---
type: C++ Function
title: FindUnstructuredGridData
resource: input_output/vtu_output.cpp#L704-L749
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/AddTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/ReplaceTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<Globals::DataLocation, std::vector<VtuOutput::UnstructuredGridData>::iterator> FindUnstructuredGridData( const TContainerType& rContainer, std::vector<VtuOutput::UnstructuredGridData>& rUnstructuredGridDataList)`

# Called by

- [AddTensorAdaptor](../../../functions/input_output/vtu_output/VtuOutput/AddTensorAdaptor.md)
- [ReplaceTensorAdaptor](../../../functions/input_output/vtu_output/VtuOutput/ReplaceTensorAdaptor.md)