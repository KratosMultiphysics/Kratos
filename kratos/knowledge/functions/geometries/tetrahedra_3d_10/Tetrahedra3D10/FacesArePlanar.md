---
type: C++ Method
title: FacesArePlanar
resource: geometries/tetrahedra_3d_10.h#L1181-L1196
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/Norm3
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FacesArePlanar() const`

# Calls

- [Norm3](../../../../functions/utilities/math_utils/Norm3.md)

# Called by

- [HasIntersection](../../../../functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/HasIntersection.md)