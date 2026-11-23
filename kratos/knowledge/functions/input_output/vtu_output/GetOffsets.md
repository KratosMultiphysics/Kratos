---
type: C++ Function
title: GetOffsets
resource: input_output/vtu_output.cpp#L258-L272
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/AddConnectivityData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`NDData<int>::Pointer GetOffsets( const std::vector<IndexType>& rWritingIndices, const TContainerType& rContainer)`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [AddConnectivityData](../../../functions/input_output/vtu_output/AddConnectivityData.md)