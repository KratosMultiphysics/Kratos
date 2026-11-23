---
type: C++ Method
title: ComputeGlobalCoordinates
resource: utilities/tessellation_utilities/curve_tessellation.h#L264-L278
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ComputeGlobalCoordinates( CoordinatesArrayType& rGlobalCoordinates, const CoordinatesArrayType& crLocaCoordinates, const GeometryType& rGeometry, bool to_surface_parameter = false )`

# Called by

- [ComputeTessellation](../../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation.md)