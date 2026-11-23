---
type: C++ Method
title: ComputeSurfaceTriangulation
resource: modeler/cad_tessellation_modeler.cpp#L161-L247
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsExactSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsApproxSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeDiscretizationError
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<BoundedMatrix<double,3,3>> CadTessellationModeler::ComputeSurfaceTriangulation( const BrepSurfaceType& rSurfaceGeometry, const std::vector<array_1d<double, 2>>& rBoundaryLoop)`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [InsertGaussPointsExactSurface](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsExactSurface.md)
- [InsertGaussPointsApproxSurface](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsApproxSurface.md)
- [ComputeDiscretizationError](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeDiscretizationError.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/SetupModelPart.md)