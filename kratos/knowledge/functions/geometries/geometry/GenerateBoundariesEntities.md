---
type: C++ Function
title: GenerateBoundariesEntities
resource: geometries/geometry.h#L2102-L2112
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GeneratePoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateFaceMaps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual GeometriesArrayType GenerateBoundariesEntities() const`

# Calls

- [GeneratePoints](../../../functions/geometries/geometry/GeneratePoints.md)

# Called by

- [GenerateFaceMaps](../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateFaceMaps.md)
- [SplitBoundary](../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary.md)