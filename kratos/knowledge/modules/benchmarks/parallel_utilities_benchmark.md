---
type: C++ Module
title: parallel_utilities_benchmark
resource: benchmarks/parallel_utilities_benchmark.cpp#L1-L120
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/utility
    resolved_by: tree-sitter
    confidence: exact
  - target: external/numeric
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/unordered-map
    resolved_by: tree-sitter
    confidence: exact
  - target: external/benchmark-benchmark-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RHSElement](../../classes/benchmarks/parallel_utilities_benchmark/RHSElement.md)
- [RHSElement](../../functions/benchmarks/parallel_utilities_benchmark/RHSElement/RHSElement.md)
- [CalculateRHS](../../functions/benchmarks/parallel_utilities_benchmark/RHSElement/CalculateRHS.md)
- [GetAccumRHSValue](../../functions/benchmarks/parallel_utilities_benchmark/RHSElement/GetAccumRHSValue.md)
- [SetAccumRHSValue](../../functions/benchmarks/parallel_utilities_benchmark/RHSElement/SetAccumRHSValue.md)
- [BM_VectorPower](../../functions/benchmarks/parallel_utilities_benchmark/BM_VectorPower.md)
- [BM_VectorReduction](../../functions/benchmarks/parallel_utilities_benchmark/BM_VectorReduction.md)
- [BM_ThreadLocalStorage](../../functions/benchmarks/parallel_utilities_benchmark/BM_ThreadLocalStorage.md)

# Imports

- `<utility>`
- `<numeric>`
- `<iostream>`
- `<unordered_map>`
- `<benchmark/benchmark.h>`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)