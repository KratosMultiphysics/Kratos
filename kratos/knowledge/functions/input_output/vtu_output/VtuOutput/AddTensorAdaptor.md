---
type: C++ Method
title: AddTensorAdaptor
resource: input_output/vtu_output.cpp#L949-L1017
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/FindUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtuOutput::AddTensorAdaptor( const std::string& rTensorAdaptorName, TTensorAdaptorPointerType pTensorAdaptor)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [FindUnstructuredGridData](../../../../functions/input_output/vtu_output/FindUnstructuredGridData.md)