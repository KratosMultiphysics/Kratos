---
type: C++ Module
title: openmp_utils
resource: utilities/openmp_utils.h#L1-L287
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/stdio-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/vector
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/omp-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/ctime
    resolved_by: tree-sitter
    confidence: exact
  - target: external/parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetNumThreads](../../functions/utilities/openmp_utils/GetNumThreads.md)
- [GetCurrentNumberOfThreads](../../functions/utilities/openmp_utils/GetCurrentNumberOfThreads.md)
- [GetNumberOfProcessors](../../functions/utilities/openmp_utils/GetNumberOfProcessors.md)
- [IsDynamic](../../functions/utilities/openmp_utils/IsDynamic.md)
- [ThisThread](../../functions/utilities/openmp_utils/ThisThread.md)
- [IsInParallel](../../functions/utilities/openmp_utils/IsInParallel.md)
- [GetCurrentTime](../../functions/utilities/openmp_utils/GetCurrentTime.md)
- [DivideInPartitions](../../functions/utilities/openmp_utils/DivideInPartitions.md)
- [PartitionedIterators](../../functions/utilities/openmp_utils/PartitionedIterators.md)
- [SetNumThreads](../../functions/utilities/openmp_utils/SetNumThreads.md)
- [PrintOMPInfo](../../functions/utilities/openmp_utils/PrintOMPInfo.md)
- [CreatePartition](../../functions/utilities/openmp_utils/CreatePartition.md)

# Imports

- `<stdio.h>`
- `<vector>`
- `<iostream>`
- `<omp.h>`
- `<ctime>`
- `"parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)