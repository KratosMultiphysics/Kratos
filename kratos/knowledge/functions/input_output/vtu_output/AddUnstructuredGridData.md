---
type: C++ Function
title: AddUnstructuredGridData
resource: input_output/vtu_output.cpp#L580-L652
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetNodesContainer
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/VtuOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddUnstructuredGridData( std::vector<VtuOutput::UnstructuredGridData>& rOutput, ModelPart& rModelPart, const IndexType EchoLevel, const bool UseSubModelParts)`

# Calls

- [FullName](../../../functions/includes/model_part/FullName.md)
- [GetNodesContainer](../../../functions/input_output/vtu_output/GetNodesContainer.md)

# Called by

- [VtuOutput](../../../functions/input_output/vtu_output/VtuOutput/VtuOutput.md)