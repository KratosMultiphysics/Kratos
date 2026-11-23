---
type: C++ Method
title: SetupModelPart
resource: modeler/cad_tessellation_modeler.cpp#L36-L126
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeBoundaryTessellation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeSurfaceTriangulation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CadTessellationModeler::SetupModelPart()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [ComputeBoundaryTessellation](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeBoundaryTessellation.md)
- [ComputeSurfaceTriangulation](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeSurfaceTriangulation.md)