---
type: C++ Method
title: ComputeDiscretizationError
resource: modeler/cad_tessellation_modeler.cpp#L332-L355
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeSurfaceTriangulation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CadTessellationModeler::ComputeDiscretizationError( const std::vector<BoundedMatrix<double,3,3>>& rGaussPointsExact, const std::vector<BoundedMatrix<double,3,3>>& rGaussPointsApprox)`

# Called by

- [ComputeSurfaceTriangulation](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeSurfaceTriangulation.md)