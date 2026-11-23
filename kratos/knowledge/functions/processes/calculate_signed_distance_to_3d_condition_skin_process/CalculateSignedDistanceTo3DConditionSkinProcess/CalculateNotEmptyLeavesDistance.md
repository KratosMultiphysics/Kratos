---
type: C++ Method
title: CalculateNotEmptyLeavesDistance
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1543-L1579
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetKey
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/CalculateCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateNotEmptyLeavesDistance(CellType* pCell)`

# Calls

- [GetKey](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetKey.md)
- [CalculateCoordinates](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalculateCoordinates.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)