---
type: C++ Method
title: InsertGaussPointsApproxSurface
resource: modeler/cad_tessellation_modeler.cpp#L286-L330
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/quadrature/Quadrature/GenerateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeSurfaceTriangulation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<BoundedMatrix<double,3,3>> CadTessellationModeler::InsertGaussPointsApproxSurface( const BrepSurfaceType& rSurfaceGeometry, const std::vector<double>& rPointsCoordinates, const std::vector<IndexType>& rTriangleConnectivities)`

# Calls

- [GenerateIntegrationPoints](../../../../functions/integration/quadrature/Quadrature/GenerateIntegrationPoints.md)

# Called by

- [ComputeSurfaceTriangulation](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeSurfaceTriangulation.md)