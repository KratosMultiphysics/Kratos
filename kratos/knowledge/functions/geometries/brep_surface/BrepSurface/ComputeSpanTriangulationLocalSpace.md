---
type: C++ Method
title: ComputeSpanTriangulationLocalSpace
resource: geometries/brep_surface.h#L545-L563
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/ComputeSpanTriangulation
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_geometries_to_python/AddGeometriesToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ComputeSpanTriangulationLocalSpace( const double u0, const double u1, const double v0, const double v1, std::vector<Matrix>& rTrianglesLocalSpace) const`

# Calls

- [ComputeSpanTriangulation](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/ComputeSpanTriangulation.md)

# Called by

- [AddGeometriesToPython](../../../../functions/python/add_geometries_to_python/AddGeometriesToPython.md)