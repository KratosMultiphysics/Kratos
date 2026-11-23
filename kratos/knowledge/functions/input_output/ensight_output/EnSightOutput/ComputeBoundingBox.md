---
type: C++ Method
title: ComputeBoundingBox
resource: input_output/ensight_output.cpp#L2793-L2841
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`BoundingBox<Point> EnSightOutput::ComputeBoundingBox( const ModelPart& rModelPart, const double Coefficient )`

# Called by

- [WriteGeometryFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometryFile.md)