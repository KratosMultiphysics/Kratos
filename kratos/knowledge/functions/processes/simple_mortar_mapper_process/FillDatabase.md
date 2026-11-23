---
type: C++ Function
title: FillDatabase
resource: processes/simple_mortar_mapper_process.h#L751-L788
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/point_object/PointObject<TObject>/pGetObject
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FillDatabase( TEntity& rGeometricalObject, KDTreeType& rTreePoints, const SizeType AllocationSize, const double SearchFactor )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)
- [pGetObject](../../../functions/spatial_containers/point_object/PointObject<TObject>/pGetObject.md)