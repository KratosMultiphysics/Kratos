---
type: C++ Function
title: GetEndianness
resource: input_output/vtu_output.cpp#L88-L91
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/IsBigEndian
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string GetEndianness()`

# Calls

- [IsBigEndian](../../../functions/input_output/vtu_output/IsBigEndian.md)

# Called by

- [WriteUnstructuredGridData](../../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)
- [WriteIntegrationPointData](../../../functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData.md)
- [PrintOutput](../../../functions/input_output/vtu_output/VtuOutput/PrintOutput.md)