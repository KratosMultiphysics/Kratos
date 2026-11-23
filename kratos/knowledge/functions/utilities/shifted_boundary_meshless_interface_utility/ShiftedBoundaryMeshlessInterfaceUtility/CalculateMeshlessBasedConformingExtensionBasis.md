---
type: C++ Method
title: CalculateMeshlessBasedConformingExtensionBasis
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L484-L677
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/GetStandardModifiedShapeFunctionsFactory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetMLSShapeFunctionsFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetRBFShapeFunctionsFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetElementSizeFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetNegativeNodeSupportCloud
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateKernelRadius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/SetNodalDistancesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ShiftedBoundaryMeshlessInterfaceUtility::CalculateMeshlessBasedConformingExtensionBasis()`

# Calls

- [SetInterfaceFlags](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetStandardModifiedShapeFunctionsFactory](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/GetStandardModifiedShapeFunctionsFactory.md)
- [GetMLSShapeFunctionsFunction](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetMLSShapeFunctionsFunction.md)
- [GetRBFShapeFunctionsFunction](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetRBFShapeFunctionsFunction.md)
- [GetElementSizeFunction](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetElementSizeFunction.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [SetNegativeNodeSupportCloud](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetNegativeNodeSupportCloud.md)
- [CalculateKernelRadius](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateKernelRadius.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [SetNodalDistancesVector](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/SetNodalDistancesVector.md)
- [ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputePositiveSideInterfaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CalculateExtensionOperator](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator.md)