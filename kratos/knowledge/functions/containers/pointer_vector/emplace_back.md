---
type: C++ Function
title: emplace_back
resource: containers/pointer_vector.h#L285-L288
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/FixedSizeMemoryPool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/AddChunk
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/Tokenize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/MergeMatrices
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/Clear
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/StlBvhTree/BuildNodeInternal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/StlBvhTree/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_pointer_vector/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/profiler/Profiler/Profiler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void emplace_back(Args&&... args)`

# Called by

- [FixedSizeMemoryPool](../../../functions/includes/fixed_size_memory_pool/FixedSizeMemoryPool/FixedSizeMemoryPool.md)
- [AddChunk](../../../functions/includes/thread_fixed_size_memory_pool/ThreadFixedSizeMemoryPool/AddChunk.md)
- [Tokenize](../../../functions/input_output/obj_io/ObjIO/Tokenize.md)
- [ColorConnectedCellsToThisCell](../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell.md)
- [MergeMatrices](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/MergeMatrices.md)
- [WriteIntermediateState](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState.md)
- [Clear](../../../functions/sources/model_part/ModelPart/Clear.md)
- [BuildNodeInternal](../../../functions/spatial_containers/stl_bvh_tree/StlBvhTree/BuildNodeInternal.md)
- [Build](../../../functions/spatial_containers/stl_bvh_tree/StlBvhTree/Build.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/containers/test_pointer_vector/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE.md)
- [FileNameDataCollector](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector.md)
- [Profiler](../../../functions/utilities/profiler/Profiler/Profiler.md)
- [ComputeTessellation](../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation.md)