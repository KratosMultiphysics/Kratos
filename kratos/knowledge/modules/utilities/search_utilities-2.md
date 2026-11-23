---
type: C++ Module
title: search_utilities
resource: utilities/search_utilities.h#L1-L1108
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/vector
    resolved_by: tree-sitter
    confidence: exact
  - target: external/numeric
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-bounding-box-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-point-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spatial-containers-point-object-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DistributedSearchInformation](../../classes/utilities/search_utilities/DistributedSearchInformation.md)
- [Reserve](../../functions/utilities/search_utilities/DistributedSearchInformation/Reserve.md)
- [Shrink](../../functions/utilities/search_utilities/DistributedSearchInformation/Shrink.md)
- [Check](../../functions/utilities/search_utilities/DistributedSearchInformation/Check.md)
- [Clear](../../functions/utilities/search_utilities/DistributedSearchInformation/Clear.md)
- [PointIsInsideBoundingBox](../../functions/utilities/search_utilities/PointIsInsideBoundingBox.md)
- [PointIsInsideBoundingBox](../../functions/utilities/search_utilities/PointIsInsideBoundingBox-2.md)
- [PointIsInsideBoundingBox](../../functions/utilities/search_utilities/PointIsInsideBoundingBox-3.md)
- [SynchronousPointSynchronization](../../functions/utilities/search_utilities/SynchronousPointSynchronization.md)
- [SynchronousPointSynchronizationWithBoundingBox](../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithBoundingBox.md)
- [SynchronousPointSynchronizationWithRecvSizes](../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithRecvSizes.md)
- [SynchronousPointSynchronizationWithRadius](../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithRadius.md)
- [PrepareSearch](../../functions/utilities/search_utilities/PrepareSearch.md)
- [PrepareOutputSearch](../../functions/utilities/search_utilities/PrepareOutputSearch.md)
- [PreparePointsSearch](../../functions/utilities/search_utilities/PreparePointsSearch.md)
- [ParallelSearch](../../functions/utilities/search_utilities/ParallelSearch.md)
- [SynchronizePointsWithRanks](../../functions/utilities/search_utilities/SynchronizePointsWithRanks.md)
- [SynchronizePointsWithBoundingBox](../../functions/utilities/search_utilities/SynchronizePointsWithBoundingBox.md)
- [SynchronizeRadius](../../functions/utilities/search_utilities/SynchronizeRadius.md)
- [CalculateNumberOfPoints](../../functions/utilities/search_utilities/CalculateNumberOfPoints.md)

# Imports

- `<vector>`
- `<numeric>`
- `"geometries/bounding_box.h"`
- `"geometries/point.h"`
- `"includes/data_communicator.h"`
- `"includes/variables.h"`
- `"includes/model_part.h"`
- `"spatial_containers/point_object.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)