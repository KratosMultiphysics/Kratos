---
type: C++ Method
title: IndexPartition
resource: utilities/parallel_utilities.h#L462-L481
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_geometrical_utilities_to_python/VectorizedFindHelper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/CheckAndCorrectZeroDiagonalValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_atomic_utilities/KRATOS_TEST_CASE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_parallel_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/particles_utilities/InterpolateValuesAtCoordinates
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IndexPartition(TIndexType Size, int Nchunks = ParallelUtilities::GetNumThreads())`

# Called by

- [VectorizedFindHelper](../../../../functions/python/add_geometrical_utilities_to_python/VectorizedFindHelper.md)
- [CheckAndCorrectZeroDiagonalValues](../../../../functions/spaces/ublas_space/CheckAndCorrectZeroDiagonalValues.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE](../../../../functions/tests/cpp_tests/utilities/test_atomic_utilities/KRATOS_TEST_CASE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_parallel_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [InterpolateValuesAtCoordinates](../../../../functions/utilities/particles_utilities/InterpolateValuesAtCoordinates.md)