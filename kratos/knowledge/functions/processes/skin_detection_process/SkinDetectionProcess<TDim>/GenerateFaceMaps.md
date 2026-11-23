---
type: C++ Method
title: GenerateFaceMaps
resource: processes/skin_detection_process.cpp#L68-L155
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GenerateBoundariesEntities
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SkinDetectionProcess<TDim>::GenerateFaceMaps( HashMapVectorIntType& rInverseFaceMap, HashMapVectorIntIdsType& rPropertiesFaceMap ) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GenerateBoundariesEntities](../../../../functions/geometries/geometry/GenerateBoundariesEntities.md)

# Called by

- [Execute](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute.md)