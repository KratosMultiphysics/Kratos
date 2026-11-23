---
type: C++ Method
title: FillNodalConnectivitiesFromElementBlockInList
resource: sources/model_part_io.cpp#L1187-L1257
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::FillNodalConnectivitiesFromElementBlockInList( ConnectivitiesContainerType& rNodalConnectivities, std::unordered_set<SizeType>& rElementsIds)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadNodalGraphFromEntitiesList](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList.md)