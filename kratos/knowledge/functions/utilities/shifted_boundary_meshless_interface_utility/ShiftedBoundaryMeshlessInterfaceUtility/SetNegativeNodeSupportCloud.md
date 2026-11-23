---
type: C++ Method
title: SetNegativeNodeSupportCloud
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L998-L1082
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetRequiredNumberOfPoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ShiftedBoundaryMeshlessInterfaceUtility::SetNegativeNodeSupportCloud( const NodeType& rNegativeNode, PointerVector<NodeType>& rCloudNodes, Matrix& rCloudCoordinates)`

# Calls

- [GetRequiredNumberOfPoints](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetRequiredNumberOfPoints.md)

# Called by

- [CalculateMeshlessBasedConformingExtensionBasis](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)