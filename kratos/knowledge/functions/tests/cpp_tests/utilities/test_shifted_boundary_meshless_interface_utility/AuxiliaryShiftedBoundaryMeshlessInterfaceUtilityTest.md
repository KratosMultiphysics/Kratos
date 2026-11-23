---
type: C++ Function
title: AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest
resource: tests/cpp_tests/utilities/test_shifted_boundary_meshless_interface_utility.cpp#L31-L69
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_shifted_boundary_meshless_interface_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest( Model& rModel, Parameters TestSettings)`

# Calls

- [GenericFindElementalNeighboursProcess](../../../../../functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess.md)
- [CalculateExtensionOperator](../../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_shifted_boundary_meshless_interface_utility/KRATOS_TEST_CASE_IN_SUITE.md)