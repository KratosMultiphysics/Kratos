---
type: C++ Method
title: SelectLowQualityElementNodes
resource: processes/mesh_local_smoothing_process.cpp#L97-L119
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/Quality
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshLocalSmoothingProcess::SelectLowQualityElementNodes()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Quality](../../../../functions/geometries/geometry/Quality.md)

# Called by

- [Execute](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/Execute.md)