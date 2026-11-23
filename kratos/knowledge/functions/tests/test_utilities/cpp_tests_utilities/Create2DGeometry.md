---
type: C++ Function
title: Create2DGeometry
resource: tests/test_utilities/cpp_tests_utilities.cpp#L35-L77
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_mpi_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_calculate_nodal_area_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_compute_nodal_gradient_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_find_nodal_h_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_find_nodal_neighbours_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/GenerateGenericModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_binbased_fast_point_locator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_constraint_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_model_part_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_result_database/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Create2DGeometry( ModelPart& rModelPart, const std::string& rEntityName, const bool Initialize, const bool Elements )`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_mpi_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_calculate_nodal_area_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_compute_nodal_gradient_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_find_nodal_h_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_find_nodal_neighbours_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [GenerateGenericModelPart](../../../../functions/tests/cpp_tests/sources/test_model_part/GenerateGenericModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [Initialize](../../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/Initialize.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_binbased_fast_point_locator/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_constraint_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_model_part_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_result_database/KRATOS_TEST_CASE_IN_SUITE.md)