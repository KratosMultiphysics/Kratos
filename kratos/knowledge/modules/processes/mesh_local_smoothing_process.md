---
type: C++ Module
title: mesh_local_smoothing_process
resource: processes/mesh_local_smoothing_process.cpp#L1-L223
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-exception-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-mesh-local-smoothing-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-measure-mesh-quality-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-nodal-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modeler-tetrahedra-ball-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-checks-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [MeshLocalSmoothingProcess](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/MeshLocalSmoothingProcess.md)
- [Execute](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/Execute.md)
- [Info](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/Info.md)
- [PrintInfo](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PrintInfo.md)
- [PrintData](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PrintData.md)
- [SelectLowQualityElementNodes](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/SelectLowQualityElementNodes.md)
- [PerformSmoothing](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PerformSmoothing.md)
- [FindOptimumPositionsAndWeights](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/FindOptimumPositionsAndWeights.md)
- [MoveNodeIfImprovesMinimumQuality](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/MoveNodeIfImprovesMinimumQuality.md)
- [InterpolateNodeOptimumPosition](../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/InterpolateNodeOptimumPosition.md)

# Imports

- `"includes/define.h"`
- `"includes/exception.h"`
- `"includes/kratos_flags.h"`
- `"processes/mesh_local_smoothing_process.h"`
- `"processes/measure_mesh_quality_process.h"`
- `"processes/find_nodal_neighbours_process.h"`
- `"modeler/tetrahedra_ball.h"`
- `"includes/checks.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)