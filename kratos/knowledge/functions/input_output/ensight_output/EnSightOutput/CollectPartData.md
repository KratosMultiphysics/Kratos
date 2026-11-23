---
type: C++ Method
title: CollectPartData
resource: input_output/ensight_output.cpp#L2198-L2263
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/GetEnSightName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/UpdatePartData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnSightOutput::CollectPartData( const ModelPart& rSubModelPart, PartData& rPartData ) const`

# Calls

- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetEnSightName](../../../../functions/input_output/ensight_output/EnSightOutput/GetEnSightName.md)

# Called by

- [UpdatePartData](../../../../functions/input_output/ensight_output/EnSightOutput/UpdatePartData.md)