---
type: C++ Module
title: shifted_boundary_meshless_interface_utility
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L1-L1181
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/containers-pointer-vector-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-global-pointer-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modified-shape-functions-triangle-2d-3-modified-shape-functions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modified-shape-functions-tetrahedra-3d-4-modified-shape-functions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-element-size-calculator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-mls-shape-functions-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-rbf-shape-functions-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-shifted-boundary-meshless-interface-utility-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [IsSplit](../../functions/utilities/shifted_boundary_meshless_interface_utility/IsSplit.md)
- [IsNegative](../../functions/utilities/shifted_boundary_meshless_interface_utility/IsNegative.md)
- [SetNodalDistancesVector](../../functions/utilities/shifted_boundary_meshless_interface_utility/SetNodalDistancesVector.md)
- [GetStandardModifiedShapeFunctionsFactory](../../functions/utilities/shifted_boundary_meshless_interface_utility/GetStandardModifiedShapeFunctionsFactory.md)
- [ShiftedBoundaryMeshlessInterfaceUtility](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/ShiftedBoundaryMeshlessInterfaceUtility.md)
- [CalculateExtensionOperator](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator.md)
- [GetDefaultParameters](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetDefaultParameters.md)
- [CalculateGradientBasedConformingExtensionBasis](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedConformingExtensionBasis](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedNonConformingExtensionBasis](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis.md)
- [SetInterfaceFlags](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags.md)
- [GetMLSShapeFunctionsAndGradientsFunction](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetMLSShapeFunctionsAndGradientsFunction.md)
- [GetMLSShapeFunctionsFunction](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetMLSShapeFunctionsFunction.md)
- [GetRBFShapeFunctionsFunction](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetRBFShapeFunctionsFunction.md)
- [GetElementSizeFunction](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetElementSizeFunction.md)
- [SetSplitElementSupportCloud](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSplitElementSupportCloud.md)
- [SetNegativeNodeSupportCloud](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetNegativeNodeSupportCloud.md)
- [CalculateKernelRadius](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateKernelRadius.md)
- [GetRequiredNumberOfPoints](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetRequiredNumberOfPoints.md)
- [SetSurrogateBoundaryNodalGradientWeights](../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSurrogateBoundaryNodalGradientWeights.md)

# Imports

- `"containers/pointer_vector.h"`
- `"includes/variables.h"`
- `"includes/global_pointer_variables.h"`
- `"modified_shape_functions/triangle_2d_3_modified_shape_functions.h"`
- `"modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"`
- `"utilities/element_size_calculator.h"`
- `"utilities/mls_shape_functions_utility.h"`
- `"utilities/rbf_shape_functions_utility.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`
- `"utilities/shifted_boundary_meshless_interface_utility.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)