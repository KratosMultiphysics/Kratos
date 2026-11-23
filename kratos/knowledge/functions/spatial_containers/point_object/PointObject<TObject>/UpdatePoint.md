---
type: C++ Method
title: UpdatePoint
resource: spatial_containers/point_object.cpp#L53-L64
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckAndPerformSearch
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/point_object/PointObject<TObject>/PointObject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/point_object/PointObject<TObject>/pSetObject
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PointObject<TObject>::UpdatePoint()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CheckAndPerformSearch](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckAndPerformSearch.md)
- [PointObject](../../../../functions/spatial_containers/point_object/PointObject<TObject>/PointObject.md)
- [pSetObject](../../../../functions/spatial_containers/point_object/PointObject<TObject>/pSetObject.md)