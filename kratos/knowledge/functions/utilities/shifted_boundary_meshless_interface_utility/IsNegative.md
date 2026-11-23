---
type: C++ Function
title: IsNegative
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L70-L81
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsNegative( const GeometryType& rGeometry, const Variable<double>& rLevelSetVariable)`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [SetInterfaceFlags](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags.md)