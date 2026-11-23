---
type: C++ Function
title: GetStandardModifiedShapeFunctionsFactory
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L111-L123
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModifiedShapeFunctionsFactoryType GetStandardModifiedShapeFunctionsFactory(const GeometryType& rGeometry)`

# Called by

- [CalculateGradientBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedNonConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis.md)