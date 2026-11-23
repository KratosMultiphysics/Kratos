---
type: C++ Method
title: SetDistancesVector
resource: utilities/embedded_skin_utility.cpp#L301-L319
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const Vector EmbeddedSkinUtility<TDim>::SetDistancesVector(const Element &rElement)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [GenerateSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin.md)
- [InterpolateMeshVariableToSkinSpecialization](../../../../functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization.md)