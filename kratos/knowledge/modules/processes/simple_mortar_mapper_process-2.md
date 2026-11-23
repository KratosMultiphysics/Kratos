---
type: C++ Module
title: simple_mortar_mapper_process
resource: processes/simple_mortar_mapper_process.h#L1-L838
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/unordered-map
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-parameters-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-mortar-classes-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-linear-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-atomic-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-exact-mortar-segmentation-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spatial-containers-point-object-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spatial-containers-spatial-containers-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [TLS](../../classes/processes/simple_mortar_mapper_process/TLS.md)
- [Map](../../functions/processes/simple_mortar_mapper_process/Map.md)
- [Info](../../functions/processes/simple_mortar_mapper_process/Info.md)
- [PrintInfo](../../functions/processes/simple_mortar_mapper_process/PrintInfo.md)
- [PrintData](../../functions/processes/simple_mortar_mapper_process/PrintData.md)
- [PerformMortarOperations](../../functions/processes/simple_mortar_mapper_process/PerformMortarOperations.md)
- [ClearIndexes](../../functions/processes/simple_mortar_mapper_process/ClearIndexes.md)
- [FillDatabase](../../functions/processes/simple_mortar_mapper_process/FillDatabase.md)

# Imports

- `<unordered_map>`
- `"processes/process.h"`
- `"includes/kratos_parameters.h"`
- `"includes/model_part.h"`
- `"includes/mortar_classes.h"`
- `"spaces/ublas_space.h"`
- `"linear_solvers/linear_solver.h"`
- `"utilities/atomic_utilities.h"`
- `"utilities/exact_mortar_segmentation_utility.h"`
- `"spatial_containers/point_object.h"`
- `"spatial_containers/spatial_containers.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)