---
type: C++ Method
title: CalculateExtensionOperator
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L184-L195
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_shifted_boundary_meshless_interface_utility/AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ShiftedBoundaryMeshlessInterfaceUtility::CalculateExtensionOperator()`

# Calls

- [CalculateGradientBasedConformingExtensionBasis](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedConformingExtensionBasis](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedNonConformingExtensionBasis](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis.md)

# Called by

- [AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest](../../../../functions/tests/cpp_tests/utilities/test_shifted_boundary_meshless_interface_utility/AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest.md)