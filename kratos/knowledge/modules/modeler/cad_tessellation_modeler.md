---
type: C++ Module
title: cad_tessellation_modeler
resource: modeler/cad_tessellation_modeler.cpp#L1-L358
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-tessellation-utilities-delaunator-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/cad-tessellation-modeler-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Create](../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/Create.md)
- [SetupModelPart](../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/SetupModelPart.md)
- [ComputeBoundaryTessellation](../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeBoundaryTessellation.md)
- [ComputeSurfaceTriangulation](../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeSurfaceTriangulation.md)
- [InsertGaussPointsExactSurface](../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsExactSurface.md)
- [InsertGaussPointsApproxSurface](../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsApproxSurface.md)
- [ComputeDiscretizationError](../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeDiscretizationError.md)

# Imports

- `"includes/define.h"`
- `"utilities/tessellation_utilities/delaunator_utilities.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`
- `"cad_tessellation_modeler.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)