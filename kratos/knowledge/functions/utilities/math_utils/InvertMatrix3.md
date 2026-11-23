---
type: C++ Function
title: InvertMatrix3
resource: utilities/math_utils.h#L465-L502
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/ShapeFunctionsIntegrationPointsGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/InverseCayleyTransform
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/InvertMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void InvertMatrix3( const TMatrix1& rInputMatrix, TMatrix2& rInvertedMatrix, double& rInputMatrixDet )`

# Called by

- [ShapeFunctionsIntegrationPointsGradients](../../../functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/ShapeFunctionsIntegrationPointsGradients.md)
- [InverseCayleyTransform](../../../functions/utilities/beam_math_utilities/InverseCayleyTransform.md)
- [InvertMatrix](../../../functions/utilities/math_utils/InvertMatrix.md)