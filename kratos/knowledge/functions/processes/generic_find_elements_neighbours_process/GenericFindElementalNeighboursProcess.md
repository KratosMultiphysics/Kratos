---
type: C++ Function
title: GenericFindElementalNeighboursProcess
resource: processes/generic_find_elements_neighbours_process.h#L77-L80
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_generic_neighbour_elements/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_shifted_boundary_meshless_interface_utility/AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_generic_find_elemental_neighbours_process/TestGenericFindElementalNeighboursProcess/testGenericFindElementalNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GenericFindElementalNeighboursProcess(ModelPart& model_part) : mrModelPart(model_part)`

# Called by

- [ExecuteInitialize](../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/ExecuteInitialize.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/processes/test_generic_neighbour_elements/KRATOS_TEST_CASE_IN_SUITE.md)
- [AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest](../../../functions/tests/cpp_tests/utilities/test_shifted_boundary_meshless_interface_utility/AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest.md)
- [testGenericFindElementalNeighboursProcess](../../../functions/tests/test_generic_find_elemental_neighbours_process/TestGenericFindElementalNeighboursProcess/testGenericFindElementalNeighboursProcess.md)