---
type: C++ Function
title: for
resource: input_output/vtu_output.cpp#L141-L155
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/GetWritingNDData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TRY for (const auto& r_model_part_data : rListOfUnstructuredGridData)`

# Calls

- [GetWritingNDData](../../../functions/input_output/vtu_output/GetWritingNDData.md)
- [pGetStorage](../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage.md)
- [visit](../../../functions/python/variable_indexing_python/visit.md)