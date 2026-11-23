---
type: C++ Module
title: example_benchmark
resource: benchmarks/example_benchmark.cpp#L1-L46
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/benchmark-benchmark-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-point-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-3d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometry-utilities-nearest-point-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [BM_TriangleNearestPoint](../../functions/benchmarks/example_benchmark/BM_TriangleNearestPoint.md)

# Imports

- `<benchmark/benchmark.h>`
- `"geometries/point.h"`
- `"geometries/triangle_3d_3.h"`
- `"utilities/geometry_utilities/nearest_point_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)