---
type: C++ Function
title: InCone
resource: utilities/geometry_utilities/brep_trimming_utilities.h#L247-L260
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/IsConvex
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/Triangulate_OPT
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static bool InCone(array_1d<double, 2>& p1, array_1d<double, 2>& p2, array_1d<double, 2>& p3, array_1d<double, 2>& p)`

# Calls

- [IsConvex](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/IsConvex.md)

# Called by

- [Triangulate_OPT](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/Triangulate_OPT.md)