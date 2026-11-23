---
type: C++ Function
title: GetGeometryTypes
resource: input_output/vtu_output.cpp#L222-L255
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

`NDData<unsigned char>::Pointer GetGeometryTypes( std::vector<IndexType>& rWritingIndices, const TContainerType& rContainer, const IndexType EchoLevel)`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [AddConnectivityData](../../../functions/input_output/vtu_output/AddConnectivityData.md)